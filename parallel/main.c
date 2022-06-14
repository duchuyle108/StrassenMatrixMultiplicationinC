#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "matrix.h"

int validate_input(int argc, char const *argv[]);
static double get_wall_seconds();

int main(int argc, char const *argv[])
{
    // Validate input
    if(validate_input(argc, argv))
        return 1;

    int N = pow(2,atoi(argv[1]));
    int n_threads = atoi(argv[2]);
    printf("Running matrix multiplication with matrix size of %d x %d, in parallel with %d threads.\n", N, N, n_threads);
    srand(time(NULL)); //set random seed

    //Initialize matrices
    etype *A, *B, *C;
    init_morton_matrix(&A, N);
    init_morton_matrix(&B, N);
    C = (etype*)malloc(N * N * sizeof(etype));
    
    /* Matrix-matrix multiplication: C = A * B  */
    double time = get_wall_seconds();
    #pragma omp parallel num_threads (n_threads)
    {   
        #pragma omp single
        mat_mul(A, B, C, N);
    }
    time = get_wall_seconds() - time;
    printf("Matrix multiplication time: %lf\n",time);

    /* Verify result */
    //Convert to 2D matrices
    etype **A2, **B2, **C2;
    mem_alloc_grid(&A2, N);
    mem_alloc_grid(&B2, N);
    mem_alloc_grid(&C2, N);
    morton_to_grid(A2, A, N);
    morton_to_grid(B2, B, N);
    morton_to_grid(C2, C, N);
    //Print 2D matrices if needed
    // print_matrix(A2, N);
    // print_matrix(B2, N);
    // print_matrix(C2, N);

    validate_result(A2, B2, C2, N);

    //Free memory allocated for matrices
    del_grid_matrix(A2, N);
    del_grid_matrix(B2, N);
    del_grid_matrix(C2, N);
    free(A);
    free(B);
    free(C);

    return 0;
}

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

int validate_input(int argc, char const *argv[]){
    //Validate number of inputs
    if(argc != 3){
        printf("Program requires 2 input args: n(matrix size) and n_threads\n");
        return 1;
    }
    //Validate matrix size
    if(!isdigit(*argv[1]) || atoi(argv[1]) <= 0){
        printf("Matrix size must be a non-zero number\n"); 
        return 1;
    }
    //Validate number of threads
    if(!isdigit(*argv[2]) || atoi(argv[2]) <= 0){
        printf("Number of threads must be a non-zero number\n"); 
        return 1;
    }
    return 0;
}