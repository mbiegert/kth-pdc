/* 
 * Homework Assignment #1 in the course
 * Parallel Computations For Large Scale Problems
 * taught at KTH Stockholm by Prof. Michael Hanke
 * in Spring 2019
 * 
 * Author: Max Biegert
 * 
 * Mandelbrot.c - generate a mandelbrot set for
 * f(z) = z^2 + d
 * 
 * generates an ascii file per default
 * 
 * use ./mandelbrot --output-mode png --output-file <filename> --palette so
 * to create png output (if compiled with -DENABLE_PNG, requires libpng)
 * available palettes for png are
 * "so": my favourite, found on stack overflow (see below in function write_png_...)
 * "interpolated": a linear interpolation of the colour values
 * "parula": requires colormaps_parula.h on compiling, filled in with the Matlab parula palette
 * 
 */

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<mpi.h>
#include<complex.h>
#include<assert.h>
#include<getopt.h>
#include<libgen.h> // for basename()
#ifdef ENABLE_PNG
#include<png.h>
#endif

// the documentation and declaration of the functions
// is postponed to the bottom in order to not disturb
// the reading flow, starting with main()
unsigned char calc_mandelbrot_pixel(const double _Complex d, const double bound, const unsigned int num_colors);
#ifdef ENABLE_PNG
int write_png_image_to_disk(unsigned char* image, const unsigned int PIXEL_WIDTH,
                            const unsigned int PIXEL_HEIGHT, const unsigned int NCOLORS,
                            const char* const filename, const char* const palette);
#endif
int write_ascii_image_to_disk(unsigned char* image, const unsigned int PIXEL_WIDTH,
                            const unsigned int PIXEL_HEIGHT, const unsigned int NCOLORS, const char* const filename);

enum {
    OUTPUTMODE_ASCII,
    OUTPUTMODE_PNG,
};

int main(int argc, char **argv) {
    int rank, size, rc;
    MPI_Status status;
    int output_mode = OUTPUTMODE_ASCII;
    char* outfile_name = "output";
    char* palette_name = "so";

    // parameters
    const unsigned int NCOLORS = 256; // Do not change this value. Producing image data
                                      // is only supported for 256 colors.
    const double BOUND = 2.0;
    const unsigned int PIXEL_HEIGHT = 8192;
    const unsigned int PIXEL_WIDTH = 8192;
    // define the viewing area
    // specify the bottom left coordinate and height, width
    // default values: -BOUND, -BOUND, 4.0, 4.0
    const double BOTTOM_LEFT_X = -2.0;
    const double BOTTOM_LEFT_Y = -1.0;
    const double IMAGE_HEIGHT = 2.0;
    const double IMAGE_WIDTH = 2.0;

    // steppings
    // We have a viewing area of [-BOUND, BOUND] x [-BOUND, BOUND]
    // and need to calculate a horizontal stepping dx and vertical
    // stepping dy and map this to actual complex numbers in the viewing area
    // (interpret R^2 \cong C)
    const double dx = IMAGE_WIDTH/(double)(PIXEL_WIDTH-1);
    const double dy = IMAGE_HEIGHT/(double)(PIXEL_HEIGHT-1);

    // initialise and setup MPI
    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // parse options, but only do so, if we our rank is 0
    if (rank == 0) {
        static struct option long_options[] =
        {
            {"output-mode", required_argument, NULL, 'm'},
            {"output-name", required_argument, NULL, 'o'},
            {"palette", required_argument, NULL, 'p'},
        };
        char ch;
        while ((ch = getopt_long(argc, argv, "m:o:p:", long_options, NULL)) != -1) {
            switch (ch) {
                case 'm':
                    if (strcmp(optarg, "png") == 0) {
                        #ifdef ENABLE_PNG
                        output_mode = OUTPUTMODE_PNG;
                        #else
                        printf("PNG support not compiled, defaulting to ascii output.");
                        #endif
                    }
                    break;
                case 'o':
                    outfile_name = basename(optarg);
                    break;
                case 'p':
                    palette_name = optarg;
                    break;
            };
        }
    }

    // We use a scatter distribution, where the process
    // of rank r does the work on row r, r+P, r+2P, ...
    //
    // First determine how many rows we will compute according to the
    // scatter distribution, to allocate the correct amount of memory needed.
    unsigned int num_rows = PIXEL_HEIGHT / size;
    const unsigned int rest = PIXEL_HEIGHT % size;
    // if the number of processes does not divide the number of rows we need
    // to also calculate the remaining rows
    if (rank < rest)
        ++num_rows;

    // This array holds the calculated colors and is a contiguous memory area
    // The first PIXEL_WIDTH uns chars are the first row this process calculated and so on
    // we start from the bottom and work our way upwards.
    unsigned char *color = (unsigned char*) malloc(PIXEL_WIDTH * num_rows * sizeof(unsigned char));
    if (!color) {
        printf("malloc failed on process of rank %d.\n", rank);
        return 1;
    }

    // Now we can start to calculate the actual colors.
    for (int y = 0; y < num_rows; ++y) {
        // now account for the correct row
        int local_row_offset = rank + (y * size);
        double dimag = local_row_offset * dy + BOTTOM_LEFT_Y;

        // now iterate over all pixels in the row and calculate the color
        for (int x = 0; x < PIXEL_WIDTH; ++x) {
            double dreal = x*dx + BOTTOM_LEFT_X;
            double _Complex d = dreal + dimag * _Complex_I;
            color[x + y*PIXEL_WIDTH] = calc_mandelbrot_pixel(d, BOUND, NCOLORS);
        }
    }

    // after calculating we have to collect all the data
    if (rank == 0) {
        // This is the collector process and we need to gather all the stuff.

        // First we allocate another array to hold the whole image
        // again the first PIXEL_WIDTH unsigned characters correspond to the
        // row at the very bottom
        unsigned char* image = malloc((unsigned long)PIXEL_WIDTH*PIXEL_HEIGHT * sizeof(unsigned char));
        if (!image) {
            printf("Failed to allocate memory for the output image.");
            return 2;
        }
        
        // now loop over all processes and receive their data, start with this process itself
        int received_rank = 0;
        unsigned int received_length = PIXEL_WIDTH*num_rows;
        for (int i = 0; i<size; ++i) {
            // only receive and overwrite the color array, if we already copied it
            if (i > 0) {
                // we receive from any source, the actual sender can later be found in "status"
                // note also that it is legitimate to use PIXEL_WIDTH*num_rows, since the process 0
                // always has the highest number of rows to calculate
                rc = MPI_Recv(color, PIXEL_WIDTH*num_rows, MPI_UNSIGNED_CHAR,
                        MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, &status);
                if (rc) {
                    printf("MPI_Recv failed. Error code: %d\n", rc);
                    return 3;
                }
                received_rank = status.MPI_SOURCE;
                rc = MPI_Get_count(&status, MPI_UNSIGNED_CHAR, &received_length);
            }

            // now loop over the received data and copy it to our image array
            // note here again, that we start counting rows from the bottom up,
            // count columns from left to right and one continous chunk in the color
            // array corresponds to one *row*
            //
            // loop over all rows received
            for (int received_row = 0; received_row * PIXEL_WIDTH < received_length; ++received_row) {
                // copy one row to the correct offset in the image array
                assert(received_rank + (unsigned long) received_row * size * PIXEL_WIDTH < (unsigned long) PIXEL_WIDTH * PIXEL_HEIGHT);
                memcpy(image + received_rank + (unsigned long) received_row * size * PIXEL_WIDTH,
                        color + (unsigned long)received_row * PIXEL_WIDTH,
                        PIXEL_WIDTH * sizeof(unsigned char));
            }
        }

        // our use for color is over
        free(color);

        // now we need to write the image array to disk
        if (output_mode == OUTPUTMODE_ASCII) {
            rc = write_ascii_image_to_disk(image, PIXEL_WIDTH, PIXEL_HEIGHT, NCOLORS, outfile_name);
        }
        #ifdef ENABLE_PNG
        else if (output_mode == OUTPUTMODE_PNG) {
            rc = write_png_image_to_disk(image, PIXEL_WIDTH, PIXEL_HEIGHT, NCOLORS, outfile_name, palette_name);
        }
        #endif

        free(image);
    }
    else {
        // If we are not the collector, just send it over to the process of rank 0
        rc = MPI_Send(color, PIXEL_WIDTH*num_rows, MPI_UNSIGNED_CHAR, 0, 0, MPI_COMM_WORLD);
        free(color);
    }

    rc = MPI_Finalize();
    return rc;
}

// this function calculates the color value for
// the complex number d, using the bound b and the number of colors N
// for the mandelbrot function f(z) = z^2 + d
unsigned char calc_mandelbrot_pixel(const double _Complex d, const double bound, const unsigned int num_colors) {
    unsigned int count = 0;
    double _Complex z = 0;
    do {
        ++count;
        z = z * z + d;
    } while (cabs(z) < bound && count < num_colors);
    return count-1;
}

// this function takes an image array of length PIXEL_WIDTH*PIXEL_HEIGHT
// and produces a file in ascii format of the name <filename>.txt
int write_ascii_image_to_disk(unsigned char* image, const unsigned int PIXEL_WIDTH, const unsigned int PIXEL_HEIGHT,
                                const unsigned int NCOLORS, const char* const filename) {
    char* full_outfile_name = malloc(strlen(filename) + 4);
    strcpy(full_outfile_name, filename);
    strcat(full_outfile_name, ".txt");
    FILE* fp = fopen(full_outfile_name, "w");
    if (!fp) {
        printf("Failed to open file \"%s\".\n", full_outfile_name);
        return 4;
    }

    // also here reverse rows
    for (unsigned int i = PIXEL_HEIGHT; i > 0; --i) {
        for (unsigned int j = 0; j < PIXEL_WIDTH; ++j)
            fprintf(fp, "%hhu ", image[j + PIXEL_HEIGHT*(i-1)]);
        
        fprintf(fp, "\n");
    }
    fclose(fp);
    return 0;
}

#ifdef ENABLE_PNG
// this function takes an image array of length PIXEL_WIDTH*PIXEL_HEIGHT
// and produces an image file in png format of the name <filename>.txt
int write_png_image_to_disk(unsigned char* image, const unsigned int PIXEL_WIDTH, const unsigned int PIXEL_HEIGHT,
                                const unsigned int NCOLORS, const char* const filename, const char* const palette) {
    int rc = 0;
    // an array containing pointers to every row of the image
    // rows[0] is a pointer to the start of the top most row
    png_bytep rows[PIXEL_HEIGHT];

    // open file and initialise libpng
    FILE *fp;
    char* full_outfile_name = malloc(strlen(filename) + 4);
    strcpy(full_outfile_name, filename);
    strcat(full_outfile_name, ".png");
    fp = fopen(full_outfile_name, "w");
    if (!fp) {
        printf("Failed to open file \"%s\".\n", full_outfile_name);
        return 4;
    }
    png_structp png_write_ptr = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!png_write_ptr) {
        printf("Failed to initialise libpng.\n");
        fclose(fp);
        return 5;
    }
    png_infop png_info_ptr = png_create_info_struct(png_write_ptr);
    if (!png_info_ptr) {
        printf("Failed to initialise libpng.\n");
        png_destroy_write_struct(&png_write_ptr, (png_infopp) NULL);
        fclose(fp);
        return 6;
    }

    if (setjmp(png_jmpbuf(png_write_ptr))) {
        printf("Failed to initialise libpng.\n");
        rc = 7;
        goto finalise;
    }

    png_init_io(png_write_ptr, fp);

    // set up all the image parameters
    png_set_IHDR(png_write_ptr, png_info_ptr, PIXEL_WIDTH, PIXEL_HEIGHT, /*bit_depth*/ 8, PNG_COLOR_TYPE_PALETTE,
        PNG_INTERLACE_NONE, PNG_COMPRESSION_TYPE_DEFAULT, PNG_FILTER_TYPE_DEFAULT);

    // specify the color palette
    #include "colormaps.h"
    if (strcmp(palette, "parula") == 0) {
        #include "colormap_parula.h"
        png_set_PLTE(png_write_ptr, png_info_ptr, color_map_256_parula, color_map_256_parula_size);
    }
    else if (strcmp(palette, "interpolated")) {
        png_set_PLTE(png_write_ptr, png_info_ptr, color_map_256_interpolated, color_map_256_interpolated_size);
    }
    else {
        png_set_PLTE(png_write_ptr, png_info_ptr, color_map_256_so, color_map_256_so_size);
    }

    // write the header and all chunks before the actual image data
    png_write_info(png_write_ptr, png_info_ptr);

    // we need to construct an array of row_pointers of type png_byte which is
    // defined as unsigned char
    // now we reverse the row ordering, as libpng expects the rows from top to bottom
    for (int row_index = 0; row_index < PIXEL_HEIGHT; ++row_index) {
        // the index formula is derived from image + PIXEL_WIDTH*PIXEL_HEIGHT - row_index*PIXEL_WIDTH
        rows[row_index] = image + (unsigned long)(PIXEL_HEIGHT - row_index) * PIXEL_WIDTH;
    }
    png_write_image(png_write_ptr, rows);
    png_write_end(png_write_ptr, NULL);

    finalise:
    png_destroy_write_struct(&png_write_ptr, &png_info_ptr);
    fclose(fp);
    return rc;
}
#endif
