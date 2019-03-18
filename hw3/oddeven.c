/* 
 * Homework Assignment #3 in the course
 * Parallel Computations For Large Scale Problems
 * taught at KTH Stockholm by Prof. Michael Hanke
 * in Spring 2019
 * 
 * Author: Max Biegert
 * 
 * oddeven.c - sort a randomly generated array
 * 
 * use ./oddeven 100 to randomly create 100 numbers in [0,1)
 * and sort them. Output will be generated to stdout
 * 
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <assert.h>
#include <time.h>

// expects left and right to be two sorted arrays of given length
// afterwards [left,right] will be two sorted arrays, where every element
// in left is smaller or equal to any element in right
void merge(double* left, int size_left, double* right, int size_right);

// do a full merge sort of the given array
void merge_sort(double* array, int length);

#ifdef DEBUG
#define DPRINTF(fmt, args...) printf(fmt, ## args)
#else
#define DPRINTF(fmt, args...)
#endif

int main(int argc, char **argv) {
    int rank, size, rc, i;
    MPI_Status status;

    // initialise and setup MPI
    rc = MPI_Init(&argc, &argv);
    rc = MPI_Comm_size(MPI_COMM_WORLD, &size);
    rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // initialise  random number generator
    srandom(rank+1);

    // Find problem size N from command line
    if (argc < 2) {
        fprintf(stderr, "No size N given\n");
        return 1;
    }
    int N = atoi(argv[1]);

    // local size, linear data distribution
    const int local_size = N/size + (rank < N%size);
    const int max_neighbour_size = N/size + (rank <= N%size);
    // allocate space for list elements
    double* local_list = malloc(local_size * sizeof(*local_list));
    double* incoming_list = malloc(max_neighbour_size * sizeof(*incoming_list));


    // data generation
    for (i = 0; i < local_size; i++) {
        local_list[i] = ((double) random())/((long)RAND_MAX+1);
    }

    // sort the local list initially
    merge_sort(local_list, local_size);

    //
    // do the actual parallel protocol
    //
    int received_size = 0;
    int evenphase = 0;

    // in even phases even processors communicate up and odd processors down
    // in odd phases even processors communicate down and odd processors up

    // we need P iterations of either even or odd steps
    for (i = 0; i < size; ++i) {
        evenphase = (i%2 == 0);
        if (rank%2 == 0) {
            // evenphase
            // but do not communicate upwards if we are the last
            if (evenphase && rank < size-1) {
                DPRINTF("Rank %d, step %d, communicating up.\n", rank, i);
                MPI_Sendrecv(local_list, local_size, MPI_DOUBLE, rank+1, 0,
                            incoming_list, max_neighbour_size, MPI_DOUBLE,
                            rank+1, 0, MPI_COMM_WORLD, &status);
                // how much did we actually receive
                MPI_Get_count(&status, MPI_DOUBLE, &received_size);
                DPRINTF("Rank %d, step %d, number of received doubles: %d\n", rank, i, received_size);
                // merge such that the smaller numbers are in the local_list
                merge(local_list, local_size, incoming_list, received_size);
            }
            // also do not communicate down
            else if (!evenphase && rank > 0) {
                DPRINTF("Rank %d, step %d, communicating down.\n", rank, i);
                MPI_Sendrecv(local_list, local_size, MPI_DOUBLE, rank-1, 1,
                            incoming_list, max_neighbour_size, MPI_DOUBLE,
                            rank-1, 1, MPI_COMM_WORLD, &status);
                // how much did we actually receive
                MPI_Get_count(&status, MPI_DOUBLE, &received_size);
                DPRINTF("Rank %d, step %d, number of received doubles: %d\n", rank, i, received_size);
                // merge such that the larger numbers are in the local list
                merge(incoming_list, received_size, local_list, local_size);
            }
        }
        // we have an odd rank
        else {
            // odd processors always communicate down, since they cannot have rank==0
            if (evenphase) {
                DPRINTF("Rank %d, step %d, communicating down.\n", rank, i);
                MPI_Sendrecv(local_list, local_size, MPI_DOUBLE, rank-1, 0,
                            incoming_list, max_neighbour_size, MPI_DOUBLE,
                            rank-1, 0, MPI_COMM_WORLD, &status);
                // how much did we actually receive
                MPI_Get_count(&status, MPI_DOUBLE, &received_size);
                DPRINTF("Rank %d, step %d, number of received doubles: %d\n", rank, i, received_size);
                // merge such that the larger numbers are in the local_list
                merge(incoming_list, received_size, local_list, local_size);
            }
            // but do not communicate up if we're the last one
            else if (rank < size-1) {
                DPRINTF("Rank %d, step %d, communicating up.\n", rank, i);
                MPI_Sendrecv(local_list, local_size, MPI_DOUBLE, rank+1, 1,
                            incoming_list, max_neighbour_size, MPI_DOUBLE,
                            rank+1, 1, MPI_COMM_WORLD, &status);
                // how much did we actually receive
                MPI_Get_count(&status, MPI_DOUBLE, &received_size);
                DPRINTF("Rank %d, step %d, number of received doubles: %d\n", rank, i, received_size);
                // merge such that the smaller numbers are in the local list
                merge(local_list, local_size, incoming_list, received_size);
            }
        }
    }

    DPRINTF("Rank %d, awaiting output slot.\n", rank);
    // generate output to fileout
    // create header if root process or wait for previous process
    char filename[15];
    FILE* fp;
    if (rank == 0) {
        // determine filename
        time_t timestamp = time(NULL);
        struct tm* now = gmtime(&timestamp);
        sprintf(filename, "output-%02d%02d", now->tm_hour, now->tm_min);
        fp = fopen(filename, "w");
    }
    else {
        MPI_Recv(filename, 15, MPI_CHAR, rank-1, 0, MPI_COMM_WORLD, &status);
        fp = fopen(filename, "a");
    }

    // print our part of the array
    for (i = 0; i < local_size; ++i) {
        fprintf(fp, "%f ", local_list[i]);
    }
    
    // signal the next process if applicable
    if (rank < size-1) {
        fclose(fp);
        MPI_Send(filename, 15, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD);
    }
    else {
        fprintf(fp, "\n");
        fclose(fp);
    }

    // cleanup and close MPI connections
    free(local_list);
    local_list = NULL;
    free(incoming_list);
    incoming_list = NULL;
    rc = MPI_Finalize();
    return rc;
}

// expects left and right to be two sorted arrays of given length
// afterwards [left,right] will be two sorted arrays, where every element
// in left is smaller or equal to any element in right
// left and right must not overlap
void merge(double* left, int size_left, double* right, int size_right) {
    // first copy the data to two temporary arrays
    double* tmp_left = malloc(size_left * sizeof(*tmp_left));
    double* tmp_right = malloc(size_right * sizeof(*tmp_right));
    memcpy(tmp_left, left, size_left * sizeof(*tmp_left));
    memcpy(tmp_right, right, size_right * sizeof(*tmp_right));

    // Iterate over the left out array and place the numbers
    int left_index = 0, right_index = 0, out_index_left = 0, out_index_right = 0;
    while (left_index < size_left && right_index < size_right) {
        if (tmp_left[left_index] < tmp_right[right_index]) {
            // use the correct output array
            if (out_index_left < size_left) {
                left[out_index_left++] = tmp_left[left_index++];
            }
            else {
                right[out_index_right++] = tmp_left[left_index++];
            }
        }
        else {
            if (out_index_left < size_left) {
                left[out_index_left++] = tmp_right[right_index++];
            }
            else {
                right[out_index_right++] = tmp_right[right_index++];
            }
        }
    }

    // copy the rest of the tmp left array if necessary
    while (left_index < size_left) {
        if (out_index_left < size_left) {
            left[out_index_left++] = tmp_left[left_index++];
        }
        else {
            assert(out_index_right < size_right);
            right[out_index_right++] = tmp_left[left_index++];
        }
    }

    // copy the rest of the tmp right array if necessary
    while (right_index < size_right) {
        if (out_index_left < size_left) {
            left[out_index_left++] = tmp_right[right_index++];
        }
        else {
            assert(out_index_right < size_right);
            right[out_index_right++] = tmp_right[right_index++];
        }
    }

    // free allocated space
    free(tmp_left);
    tmp_left = NULL;
    free(tmp_right);
    tmp_right = NULL;
}

// do a full merge sort (for initial local sorting)
//
void merge_sort(double* array, int length) {
    if (length > 1) {
        int m = length/2;
        merge_sort(array, m);
        merge_sort(array+m, length - m);

        // merge the two sub arrays
        merge(array, m, array+m, length-m);
    }
}
