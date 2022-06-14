#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include "matrix.h"

#define CUT_OFF 7 //cut-off point at matrix size = 2^CUT_OFF

extern etype ***M, ***T; //Strassen auxiliary and temporary matrices
extern etype **tempA, **tempB, **tempC; //Temporary matrix used inside cut-off condition

int validate_input(int argc, char const *argv[]);
static double get_wall_seconds();

int main(int argc, char const *argv[])
{
    // Validate input
    if(validate_input(argc, argv))
        return 1;

    int power = atoi(argv[1]);
    int N = pow(2,power); //matrix size
    printf("Matrix size: %d\n", N);
    srand(time(NULL)); //set random seed

    //Initialize matrices
    etype *A, *B, *C;
    init_morton_matrix(&A, N);
    init_morton_matrix(&B, N);
    C = (etype*)malloc(N * N * sizeof(etype));
    
    /* Main calculation part for C = A * B */
    double time = get_wall_seconds();
    int i, j, sub_block, cutoff_matrix_size;
    
    //If input matrix size > cut-off point
    if(power > CUT_OFF){
        //Define auxiliary arrays of matricies: M1->M7 and T1->T3 as in theory
        M = (etype***)malloc(7 * sizeof(etype**));
        T = (etype***)malloc(3 * sizeof(etype**));

        /*This part calculate and allocate memory for each auxiliary variable in each layer
        Eg: - For input N = 4; layer 0: each of variables M or T have 2^(4-1) * 2^(4-1) memory blocks
        which is equivalent to the size of 1/4 of the original matrices.
            - M[2][1] equivalent to matrix M2 used in recursive layer 1;
        */
        for(i = 0; i < 7; i++){
            M[i] = (etype**)malloc((power - 1) * sizeof(etype*));
            for(j = 0; j < power - 1; j++){
                sub_block = pow(2, 2 * (power - 1 - j));
                M[i][j] = (etype*)malloc(sub_block* sizeof(etype));
            }
        }
        for(i = 0; i < 3; i++){
            T[i]= (etype**)malloc((power - 1) * sizeof(etype*));
            for(j = 0; j < power - 1; j++){
                sub_block = pow(2, 2 * (power - 1 - j));
                T[i][j] = (etype*)malloc(sub_block* sizeof(etype));
            }
        }
        //Allocate memory for temporary matrices inside cut-off condition
        cutoff_matrix_size = pow(2, CUT_OFF);
        mem_alloc_grid(&tempA, cutoff_matrix_size);
        mem_alloc_grid(&tempB, cutoff_matrix_size);
        mem_alloc_grid(&tempC, cutoff_matrix_size);
    } else{
        cutoff_matrix_size = pow(2, power);
        mem_alloc_grid(&tempA, cutoff_matrix_size);
        mem_alloc_grid(&tempB, cutoff_matrix_size);
        mem_alloc_grid(&tempC, cutoff_matrix_size);
    }

    //Multiplication function call starting with layer 0
    mat_mul(A, B, C, N, 0, pow(2, CUT_OFF));

    time = get_wall_seconds() - time;
    printf("Matrix multiplication time = %lf\n",time);

    /* VERIFY RESULT */
    //Create 2D matrices
    etype **A2, **B2, **C2;
    mem_alloc_grid(&A2, N);
    mem_alloc_grid(&B2, N);
    mem_alloc_grid(&C2, N);

    //Convert Morton-ordering to 2D
    morton_to_grid(A2, A, N);
    morton_to_grid(B2, B, N);
    morton_to_grid(C2, C, N);

    //Visualize matrices if needed
    // print_grid_matrix(A2, N);
    // print_grid_matrix(B2, N);
    // print_grid_matrix(C2, N);

    //Validate result
    validate_result(A2, B2, C2, N);

    /* Free allocated memory */
    //Free memory allocated for matrices
    del_grid_matrix(A2, N);
    del_grid_matrix(B2, N);
    del_grid_matrix(C2, N);
    free(A);
    free(B);
    free(C);
    //Free M, T
    if(power > CUT_OFF){
        for(i = 0; i < 7; i++){
            for(j = 0; j < power - 1; j++){
                free(M[i][j]);
            }
            free(M[i]);
        }
        for(i = 0; i < 3; i++){
            for(j = 0; j < power - 1; j++){
                free(T[i][j]);
            }
            free(T[i]);
        }
        free(M);
        free(T);
    }
    //Free temporary matrix used inside cut-off condition
    del_grid_matrix(tempA, cutoff_matrix_size);
    del_grid_matrix(tempB, cutoff_matrix_size);
    del_grid_matrix(tempC, cutoff_matrix_size);

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
    if(argc != 2){
        printf("Program requires 1 input arg: n(matrix size)\n");
        return 1;
    }
    //Validate matrix size
    if(!isdigit(*argv[1]) || atoi(argv[1]) <= 0){
        printf("Matrix size must be a non-zero number\n"); 
        return 1;
    }
    return 0;
}