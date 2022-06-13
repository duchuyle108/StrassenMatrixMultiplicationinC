#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

typedef double etype;

static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

int validate_input(int argc, char const *argv[]);
void validate_result(int N, etype **A, etype **B, etype **C);
void init_matrix(etype ***arr, int N);
void del_matrix(etype **mat, int N);
void mem_alloc(etype ***mat, int N);
etype rand_etype(etype min, etype max);
void print_matrix(etype **mat, int N);
void mat_mul(etype **A, etype **B, etype **res, int N);
void mat_add(etype **A, etype **B, etype **res, int N);
void mat_sub(etype **A, etype **B, etype **res, int N);

int main(int argc, char const *argv[])
{
    // Validate input
    if(validate_input(argc, argv))
        return 1;

    int N = pow(2,atoi(argv[1]));
    printf("Matrix size: %d\n", N);
    srand(time(NULL));
    //Initialize matrices
    etype **A, **B, **C;
    init_matrix(&A, N);
    init_matrix(&B, N);
    mem_alloc(&C, N);
    
    //Multiply
    double start_time, exec_time;
    start_time = get_wall_seconds();
    mat_mul(A, B, C, N);
    exec_time = get_wall_seconds() - start_time;
    printf("Matrix multiplication time = %lf\n",exec_time);
    validate_result(N, A, B, C);

    //Free memory allocated for matrices
    del_matrix(A, N);
    del_matrix(B, N);
    del_matrix(C, N);

    return 0;
}

void mat_mul(etype **A, etype **B, etype **res, int N){
    if(N == 1)
        **res = **A * **B;
    else {
        etype **A11, **A12, **A21, **A22;
        etype **B11, **B12, **B21, **B22;
        etype **C11, **C12, **C21, **C22;
        etype **M1, **M2, **M3, **M4, **M5, **M6, **M7;
        etype **T1, **T2, **T3, **T4, **T5, **T6, **T7, **T8, **T9, **T10;
        
        //Allocate memory
        mem_alloc(&A11, N/2); mem_alloc(&A12, N/2); mem_alloc(&A21, N/2); mem_alloc(&A22, N/2);
        mem_alloc(&B11, N/2); mem_alloc(&B12, N/2); mem_alloc(&B21, N/2); mem_alloc(&B22, N/2);
        mem_alloc(&C11, N/2); mem_alloc(&C12, N/2); mem_alloc(&C21, N/2); mem_alloc(&C22, N/2);
        mem_alloc(&M1, N/2); mem_alloc(&M2, N/2); mem_alloc(&M3, N/2); mem_alloc(&M4, N/2); 
        mem_alloc(&M5, N/2); mem_alloc(&M6, N/2); mem_alloc(&M7, N/2);
        mem_alloc(&T1, N/2); mem_alloc(&T2, N/2); mem_alloc(&T3, N/2); mem_alloc(&T4, N/2); 
        mem_alloc(&T5, N/2); mem_alloc(&T6, N/2); mem_alloc(&T7, N/2); mem_alloc(&T8, N/2);
        mem_alloc(&T9, N/2); mem_alloc(&T10, N/2);
        
        //Assign sub-matrix
        int i, j;
        for (i = 0; i < N/2; i++)
            for (j = 0; j < N/2; j++){
                A11[i][j] = A[i][j];
                A12[i][j] = A[i][j+N/2];
                A21[i][j] = A[i+N/2][j];
                A22[i][j] = A[i+N/2][j+N/2];
                B11[i][j] = B[i][j];
                B12[i][j] = B[i][j+N/2];
                B21[i][j] = B[i+N/2][j];
                B22[i][j] = B[i+N/2][j+N/2];
            }

        //Calculate M1
        mat_add(A11, A22, T1, N/2);
        mat_add(B11, B22, T2, N/2);
        mat_mul(T1, T2, M1, N/2);

        //Calculate M2
        mat_add(A21, A22, T3, N/2);
        mat_mul(T3, B11, M2, N/2);

        //Calculate M3
        mat_sub(B12, B22, T4, N/2);
        mat_mul(A11, T4, M3, N/2);

        //Calculate M4
        mat_sub(B21, B11, T5, N/2);
        mat_mul(A22, T5, M4, N/2);

        //Calculate M5
        mat_add(A11, A12, T6, N/2);
        mat_mul(T6, B22, M5, N/2);

        //Calculate M6
        mat_sub(A21, A11, T7, N/2);
        mat_add(B11, B12, T8, N/2);
        mat_mul(T7, T8, M6, N/2);

        //Calculate M7
        mat_sub(A12, A22, T9, N/2);
        mat_add(B21, B22, T10, N/2);
        mat_mul(T9, T10, M7, N/2);

        //Calculate C11
        mat_add(M1, M4, T1, N/2);
        mat_sub(M5, M7, T2, N/2);
        mat_sub(T1, T2, C11, N/2);

        //Calculate C12
        mat_add(M3, M5, C12, N/2);

        //Calculate C21
        mat_add(M2, M4, C21, N/2);

        //Calculate C22
        mat_sub(M1, M2, T1, N/2);
        mat_add(M3, M6, T2, N/2);
        mat_add(T1, T2, C22, N/2);

        for (i = 0; i < N/2; i++)
            for (j = 0; j < N/2; j++){
                res[i][j] = C11[i][j];
                res[i][j+N/2] = C12[i][j];
                res[i+N/2][j] = C21[i][j];
                res[i+N/2][j+N/2] = C22[i][j];
            }

        //Delete matrices
        del_matrix(A11, N/2); del_matrix(A12, N/2); del_matrix(A21, N/2); del_matrix(A22, N/2);
        del_matrix(B11, N/2); del_matrix(B12, N/2); del_matrix(B21, N/2); del_matrix(B22, N/2);
        del_matrix(C11, N/2); del_matrix(C12, N/2); del_matrix(C21, N/2); del_matrix(C22, N/2);
        del_matrix(M1, N/2); del_matrix(M2, N/2); del_matrix(M3, N/2); del_matrix(M4, N/2); 
        del_matrix(M5, N/2); del_matrix(M6, N/2); del_matrix(M7, N/2);
        del_matrix(T1, N/2); del_matrix(T2, N/2); del_matrix(T3, N/2); del_matrix(T4, N/2); 
        del_matrix(T5, N/2); del_matrix(T6, N/2); del_matrix(T7, N/2); del_matrix(T8, N/2); 
        del_matrix(T9, N/2); del_matrix(T10, N/2);
    }

}

void mat_add(etype **A, etype **B, etype **res, int N){
    int i, j;
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
            res[i][j] = A[i][j] + B[i][j];
}

void mat_sub(etype **A, etype **B, etype **res, int N){
    int i, j;
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
            res[i][j] = A[i][j] - B[i][j];
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

void validate_result(int N, etype **A, etype **B, etype **C){
    etype **C_test = (etype**)malloc(N * sizeof(etype*));
    int i, j, k;
    for(i = 0; i < N; i++){
        C_test[i] = (etype*)malloc(N * sizeof(etype));
        memset(C_test[i], 0, N * sizeof(etype));
    }
    
    //Naive matrix multiplication
    for(i = 0; i < N; i++)
        for(k = 0; k < N; k++)
            for(j = 0; j < N; j++)
                C_test[i][j] += A[i][k] * B[k][j];
    etype max_diff = 0;
    //Compare results
    for(i = 0; i < N; i++)
        for(j = 0; j < N; j++)
            if(C[i][j] - C_test[i][j] > max_diff)
                max_diff = C[i][j] - C_test[i][j];
    
    del_matrix(C_test, N);
    printf("Max diff in an element is : %3.15lf\n", max_diff);
}

void mem_alloc(etype ***mat, int N){
    int i;
    *mat = (etype**)malloc(N * sizeof(etype*));
    for(i = 0; i < N; i++){
        (*mat)[i] = (etype*)malloc(N * sizeof(etype));
    }
}

void init_matrix(etype ***mat, int N){
    int i, j;
    //Allocate memory
    mem_alloc(mat, N);

    // Initialize matrix values
    for(i =0; i < N; i++)
        for(j=0; j < N; j++)
            (*mat)[i][j] = rand_etype(-10, 10);

}

void del_matrix(etype **mat, int N){
    int i;
    for(i = 0; i < N; i++){
        free(mat[i]);
    }
    free(mat);
}

etype rand_etype(etype min, etype max){
    etype range = max - min;
    return min + (rand() * range / RAND_MAX);

}

void print_matrix(etype **mat, int N){
    int i, j;
    for(i =0; i < N; i++){
        for(j=0; j < N; j++)
            printf("%lf ",mat[i][j]);
        printf("\n");
    }
}


