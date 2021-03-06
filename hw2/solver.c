/* 
 * Homework Assignment #2 in the course
 * Parallel Computations For Large Scale Problems
 * taught at KTH Stockholm by Prof. Michael Hanke
 * in Spring 2019
 * 
 * Author: Max Biegert
 * 
 * solver.c - compute a solution u of the 1D differential equation
 * u''(x) + r(x)u(x) = f(x) for 0 < x < 1, r(x) <= 0
 * u(0) = u(1) = 0
 * on a linear grid of stepsize h.
 * 
 * The equation will be solved by a Jacobi iteration and data
 * will be linearly distributed among the nodes.
 * The final result will be written to an ascii text file, where each
 * row is one value. The first line contains the number of iterations.
 * 
 * The python program "generate_plot.py" can be used to create a plot of
 * the computed data.
 * 
 * You can use ./solver --output|-o "filename" to write the output to the
 * specified file.
 */

#include<stdio.h>
#include<stdlib.h>
#include<mpi.h>
#include<assert.h>
#include<getopt.h>
#include<libgen.h> // for basename()
#include<math.h>

static const double PI = 3.14159265358979323846;

// "stiffness" parameter
static const double k = 2.0;

// this function is used to plot the error
double u_exact(double x) {
    return sin(2*PI*x*exp(k*(x-1)));
}

// choose a function r(x) <= 0
double r(double x) {
    //double ret = -1;
    double ret = -x;
    assert(ret <= 0);
    return ret;
}

// choose a right hand side f(x)
double f(double x) {
    // this right hand side produces a nice sine curve
    //return -(2*PI) * (2*PI) * sin(2*PI*x) - x*sin(2*PI*x);

    const double cosine = cos(2*PI*x*exp(k*(x-1)));
    const double sine = sin(2*PI*x*exp(k*(x-1)));
    const double expo = exp(k*(x-1));
    return cosine*(2*PI*k*expo*(2+k*x)) - sine*pow(2*PI*expo*(k*x+1), 2) + r(x)*u_exact(x);
}

double min(double a, double b) {
    if (a < b) {
        return a;
    }
    else {
        return b;
    }
}

int main(int argc, char **argv) {

    // parameters, adjust these to your liking before compiling
    const double lower_limit = 0.0;
    const double upper_limit = 1.0;
    // we expect lower_limit < upper_limit

    // grid_size specifies the number of evaluations of the function
    // the first grid point corresponds to lower_limit and the last
    // grid point to upper_limit
    const size_t grid_size = 1002; // > 1
    const size_t iteration_steps = 10000000;

    // default output file name
    char* outfile_name = "udata.txt";

    // calculate the step size given the above parameters
    const double grid_step_size = (upper_limit - lower_limit) / (grid_size - 1);

    // initialise some counters
    int i,j = 0;

    // parse the output file name from the command line
    static struct option long_options[] =
    {
        {"output", required_argument, NULL, 'o'},
    };
    char ch;
    while ((ch = getopt_long(argc, argv, "m:o:p:", long_options, NULL)) != -1) {
        switch (ch) {
            case 'o':
                outfile_name = basename(optarg);
                break;
        };
    }


    // initialise and setup MPI
    int rank, size, rc;
    MPI_Status status;
    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // First determine how many points we will compute according to the
    // linear distribution, to allocate the correct amount of memory needed.
    // The overlap a=2, so we need one ghost point to the left and one to the right.
    // We use a little trick to ease the computations (and stay in sync with the
    // provided skeleton). We linearly distribute the number of *inner* grid points,
    // that is excluding the end points. Then every the process gets two ghost points
    // and in the case of the process 0 and size-1 these correspond to the end points
    // and are set to 0 according to the boundary condition.
    size_t array_size = (grid_size - 2) / size + 2;

    // if the number of processes does not divide the number of points evenly we need
    // to also distribute the remaining points
    if (rank < (grid_size-2) % size)
        ++array_size;

    // allocate arrays to hold the computed values and  initialise them to zero
    // we need to arrays, since we cannot update the array in place
    double* u_old = calloc(array_size, sizeof(*u_old));
    double* u_new = calloc(array_size, sizeof(*u_new));
    if (!u_old || !u_new) {
        return 1;
    }

    // instead of evaluating the functions r and f every time, we do it
    // once and write it in an array.
    double* ff = malloc((array_size-2) * sizeof(*ff));
    double* rr = malloc((array_size-2) * sizeof(*rr));

    // before we can evaulate f and r we first need to determine the actual
    // x value that corresponds to the first of our inner grid points
    double offset = lower_limit + 
            (rank * (grid_size-2) / size + min(rank, (grid_size-2)%size))*grid_step_size;
    offset = offset + grid_step_size; // this step is necessary since we want to
                                      // evaluate at the inner grid points
    for (i = 0; i < array_size-2; ++i) {
        ff[i] = f(offset + i*grid_step_size);
        rr[i] = r(offset + i*grid_step_size);
    }

    for (i = 0; i < iteration_steps; ++i) {
        // calculation
        // 
        // now we can start to calculate the actual values using the
        // three point stencil
        for (j = 1; j < array_size - 1; ++j) {
            // apply the update function with the stencil once to
            // the whole array
            u_new[j] = u_old[j-1] + u_old[j+1] - grid_step_size*grid_step_size*ff[j];
            u_new[j] = u_new[j] / (2.0 - grid_step_size*grid_step_size*rr[j]);
        }

        // swap the arrays, note that we only reassign addresses here
        double* tmp = u_new;
        u_new = u_old;
        u_old = tmp;

        // communication
        // 
        // next we need to exchange the overlap with the previous and next processor
        // we do this in a red-black pattern, where even ranks are black and
        // odd ranks are red
        // basically even ranks send first and recv afterwards and odds do it the
        // other way round
        if (size > 1 && i < iteration_steps) {
            // first handle the two edge cases
            if (rank == 0) {
                // only send up
                rc = MPI_Send(u_old+array_size-2, 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
                rc = MPI_Recv(u_old+array_size-1, 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
            }
            else if (rank == size-1) {
                if (rank % 2 == 0) {
                    // we need to send first
                    rc = MPI_Send(u_old+1, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
                    rc = MPI_Recv(u_old, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
                }
                else {
                    // we need to receive first
                    rc = MPI_Recv(u_old, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
                    rc = MPI_Send(u_old+1, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
                }
            }
            else {
                if (rank % 2 == 0) {
                    // send up -> recv up -> send down -> recv down
                    rc = MPI_Send(u_old+array_size-2, 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
                    rc = MPI_Recv(u_old+array_size-1, 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
                    rc = MPI_Send(u_old+1, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
                    rc = MPI_Recv(u_old, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
                }
                else {
                    // recv down -> send down -> recv up -> recv down
                    rc = MPI_Recv(u_old, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &status);
                    rc = MPI_Send(u_old+1, 1, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD);
                    rc = MPI_Recv(u_old+array_size-1, 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &status);
                    rc = MPI_Send(u_old+array_size-2, 1, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD);
                }
            }
        }
    }

    // @TODO: This for loop is wrong! It produces too high offsets.
    // lastly calculate the error
    double* error = malloc(array_size * sizeof(*error));
    for (i = 0; i<array_size; ++i) {
        error[i] = fabs(u_exact(offset+(i-1)*grid_step_size) - u_old[i]);
    }

    // we are done and only have to output
    char mode = 'a';
    int temp;
    if (rank == 0) {
        mode = 'w';
    }

    if (rank > 0) {
        // wait for the previous process to be finished
        MPI_Recv(&temp, 1, MPI_INTEGER, rank-1, 0, MPI_COMM_WORLD, &status);
    }

    // open the file
    FILE* fp = fopen(outfile_name, &mode);
    if (!fp) {
        printf("Failed to open file \"%s\" on processor of rank %d.\n", outfile_name, rank);
        return -1;
    }

    if (rank == 0) {
        fprintf(fp, "%d\n", iteration_steps);
        fprintf(fp, "%f %f\n", u_old[0], error[0]);
    }

    // write the array out
    for (i = 1; i<array_size-1; ++i) {
        fprintf(fp, "%f %f\n", u_old[i], error[i]);
    }

    if (rank == size-1) {
        fprintf(fp, "%f %f\n", u_old[array_size-1], error[array_size-1]);
    }

    fclose(fp);

    // notify the next process, if applicable
    if (rank < size-1) {
        MPI_Send(&rank, 1, MPI_INTEGER, rank+1, 0, MPI_COMM_WORLD);
    }

    // clean up
    free(u_new);
    free(u_old);
    free(ff);
    free(rr);

    rc = MPI_Finalize();
    return rc;
}
