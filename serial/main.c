#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <time.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

typedef double etype;

etype ***M;
etype ***T;

etype **tempA, **tempB, **tempC;


static double get_wall_seconds() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
  return seconds;
}

void morton_to_grid(etype **grid, etype *morton, int N);
void grid_to_morton(etype **grid, etype *morton, int N);
int validate_input(int argc, char const *argv[]);
void validate_result(etype **A, etype **B, etype **C, int N);
void init_matrix(etype **arr, int N);
void del_matrix(etype **mat, int N);
void mem_alloc(etype **mat, int N);
void mem_alloc_2d(etype ***mat, int N);
etype rand_etype(etype min, etype max);
void print_matrix(etype **mat, int N);
void mat_mul(etype *A, etype *B, etype *res, int N, int step);
void mat_add(etype *A, etype *B, etype *res, int N);
void mat_sub(etype *A, etype *B, etype *res, int N);

int main(int argc, char const *argv[])
{
    // Validate input
    if(validate_input(argc, argv))
        return 1;

    int power = atoi(argv[1]);
    int cut_off = 7;
    int N = pow(2,power);
    int temp_mat_size;
    printf("Matrix size: %d\n", N);
    srand(time(NULL));
    //Initialize matrices
    etype *A, *B, *C;
    init_matrix(&A, N);
    init_matrix(&B, N);
    C = (etype*)malloc(N * N * sizeof(etype));
    
    //Multiply
    int i, j, sub_block;
    if(power > cut_off){
        M = (etype***)malloc(7 * sizeof(etype**));
        T = (etype***)malloc(3 * sizeof(etype**));

        for(i = 0; i < 7; i++){
            M[i] = (etype**)malloc((power - cut_off) * sizeof(etype*));
            for(j = 0; j < power - cut_off; j++){
                sub_block = pow(2, 2 * (power - 1 - j));
                M[i][j] = (etype*)malloc(sub_block* sizeof(etype));
            }
        }

        for(i = 0; i < 3; i++){
            T[i]= (etype**)malloc((power - cut_off) * sizeof(etype*));
            for(j = 0; j < power - cut_off; j++){
                sub_block = pow(2, 2 * (power - 1 - j));
                T[i][j] = (etype*)malloc(sub_block* sizeof(etype));
            }
        }
        temp_mat_size = pow(2, cut_off);
        mem_alloc_2d(&tempA, temp_mat_size);
        mem_alloc_2d(&tempB, temp_mat_size);
        mem_alloc_2d(&tempC, temp_mat_size);
    }
    
    double start_time, exec_time;
    start_time = get_wall_seconds();
    mat_mul(A, B, C, N, 0);
    exec_time = get_wall_seconds() - start_time;
    printf("Matrix multiplication time = %lf\n",exec_time);

    //Convert to 2D matrix
    etype **A2, **B2, **C2;
    mem_alloc_2d(&A2, N);
    mem_alloc_2d(&B2, N);
    mem_alloc_2d(&C2, N);

    morton_to_grid(A2, A, N);
    morton_to_grid(B2, B, N);
    morton_to_grid(C2, C, N);

    // // print_matrix(A2, N);
    // // printf("\n");
    // // print_matrix(B2, N);
    // // printf("\n");
    // // print_matrix(C2, N);
    // // printf("\n");

    validate_result(A2, B2, C2, N);

    del_matrix(A2, N);
    del_matrix(B2, N);
    del_matrix(C2, N);

    //Free memory allocated for matrices
    free(A);
    free(B);
    free(C);

    //Free M, T
    if(power > cut_off){
        for(i = 0; i < 7; i++){
            for(j = 0; j < power - cut_off; j++){
                free(M[i][j]);
            }
            free(M[i]);
        }
        for(i = 0; i < 3; i++){
            for(j = 0; j < power - cut_off; j++){
                free(T[i][j]);
            }
            free(T[i]);
        }
        free(M);
        free(T);
        del_matrix(tempA, temp_mat_size);
        del_matrix(tempB, temp_mat_size);
        del_matrix(tempC, temp_mat_size);
    }

    return 0;
}

void mat_mul(etype *A, etype *B, etype *res, int N, int step){
    if(N <= 512){
        int i,j,k;
        for(i = 0; i < N; i++)
            for(k = 0; k < N; k++)
                tempC[i][k] = 0;

        morton_to_grid(tempA, A, N);
        morton_to_grid(tempB, B, N);

        for(i = 0; i < N; i++)
            for(k = 0; k < N; k++)
                for(j = 0; j < N; j++)
                    tempC[i][j] += tempA[i][k] * tempB[k][j];
        
        grid_to_morton(tempC, res, N);
    }
    else {
        etype *A11, *A12, *A21, *A22;
        etype *B11, *B12, *B21, *B22;
        etype *C11, *C12, *C21, *C22;
        etype *M1, *M2, *M3, *M4, *M5, *M6, *M7;
        etype *T1, *T2, *T3;

        int n = N/2;
        int N1 = N*N/4;
        int N2 = N*N/2;
        int N3 = 3*N*N/4;
        
        // //Allocate memory
        // mem_alloc(&M1, N1); mem_alloc(&M2, N1); mem_alloc(&M3, N1); mem_alloc(&M4, N1); 
        // mem_alloc(&M5, N1); mem_alloc(&M6, N1); mem_alloc(&M7, N1);
        // mem_alloc(&T1, N1); mem_alloc(&T2, N1); mem_alloc(&T3, N1); mem_alloc(&T4, N1); 
        
        int i;
        //Assign sub-matrix
        A11 = A; A12 = &A[N1]; A21 = &A[N2]; A22 = &A[N3];

        B11 = B; B12 = &B[N1]; B21 = &B[N2]; B22 = &B[N3];

        C11 = res; C12 = &res[N1]; C21 = &res[N2]; C22 = &res[N3];

        M1 = M[0][step]; M2 = M[1][step]; M3 = M[2][step]; M4 = M[3][step]; M5 = M[4][step];
        M6 = M[5][step]; M7 = M[6][step]; T1 = T[0][step]; T2 = T[1][step]; T3 = T[2][step];

        //Calculate M1
        mat_add(A11, A22, T1, n);
        mat_add(B11, B22, M1, n);
        mat_mul(T1, M1, M1, n, step + 1);
        
        //Calculate M2
        mat_add(A21, A22, M2, n);
        mat_mul(M2, B11, M2, n, step + 1);

        //Calculate M3
        mat_sub(B12, B22, M3, n);
        mat_mul(A11, M3, M3, n, step + 1);

        //Calculate M4
        mat_sub(B21, B11, M4, n);
        mat_mul(A22, M4, M4, n, step + 1);

        //Calculate M5
        mat_add(A11, A12, M5, n);
        mat_mul(M5, B22, M5, n, step + 1);

        //Calculate M6
        mat_sub(A21, A11, T2, n);
        mat_add(B11, B12, M6, n);
        mat_mul(T2, M6, M6, n, step + 1);

        //Calculate M7
        mat_sub(A12, A22, T3, n);
        mat_add(B21, B22, M7, n);
        mat_mul(T3, M7, M7, n, step + 1);

        for(i = 0; i < N1; i++){
            C11[i] = M1[i] + M4[i] - M5[i] + M7[i];
            C12[i] = M3[i] + M5[i];
            C21[i] = M2[i] + M4[i];
            C22[i] = M1[i] - M2[i] + M3[i] + M6[i];
        }

        //Delete matrices
        // free(M1); free(M2); free(M3); free(M4); free(M5); free(M6); free(M7);
        // free(T1); free(T2); free(T3); free(T4);
    }

}

void mat_add(etype *A, etype *B, etype *res, int N){
    int i;
    int n = N * N;
    for(i = 0; i < n; i++)
            res[i] = A[i] + B[i];
}

void mat_sub(etype *A, etype *B, etype *res, int N){
    int i;
    int n = N * N;
    for(i = 0; i < n; i++)
            res[i] = A[i] - B[i];
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

void validate_result(etype **A, etype **B, etype **C, int N){
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

    printf("Max diff in an element is : %3.15lf\n", max_diff);
    del_matrix(C_test, N);
}

void init_matrix(etype **mat, int N){
    int i;
    //Allocate memory
    int n = N * N;
    *mat = (etype*)malloc(n * sizeof(etype));

    // Initialize matrix values
    for(i =0; i < n; i++)
            (*mat)[i] = rand_etype(-10, 10);

}

void mem_alloc(etype **mat, int N){
    *mat = (etype*)malloc(N * N * sizeof(etype));
}

void mem_alloc_2d(etype ***mat, int N){
    int i;
    *mat = (etype**)malloc(N * sizeof(etype*));
    for(i = 0; i < N; i++)
        (*mat)[i] = (etype*)malloc(N * sizeof(etype));
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

// Convert a Morton-order matrix to a 2D matrix
void morton_to_grid(etype **grid, etype *morton, int N){
    if(N == 2){
        grid[0][0] = morton[0];
        grid[0][1] = morton[1];
        grid[1][0] = morton[2];
        grid[1][1] = morton[3];
    }
    else{
        int n = N/2;
        //Divide morton-order matrix into 4 parts
        etype *morton2 = &morton[N*N/4];
        etype *morton3 = &morton[N*N/2];
        etype *morton4 = &morton[3*N*N/4];
        int i;

        //Divide 2D matrix into 4 parts
        etype **grid11, **grid12, **grid21, **grid22;
        grid11 = (etype**)malloc(n*sizeof(etype*));
        grid12 = (etype**)malloc(n*sizeof(etype*));
        grid21 = (etype**)malloc(n*sizeof(etype*));
        grid22 = (etype**)malloc(n*sizeof(etype*));
        
        for(i = 0; i < n; i++){
            grid11[i] = &grid[i][0];
            grid12[i] = &grid[i][n];
            grid21[i] = &grid[i+n][0];
            grid22[i] = &grid[i+n][n];
        }
        
        //Recursive call in each part
        morton_to_grid(grid11, morton, n);
        morton_to_grid(grid12, morton2, n);
        morton_to_grid(grid21, morton3, n);
        morton_to_grid(grid22, morton4, n);

        free(grid11); free(grid12); free(grid21); free(grid22);
    }
}

// Convert a 2D matrix to a Morton-order matrix
void grid_to_morton(etype **grid, etype *morton, int N){
    if(N == 2){
        morton[0] = grid[0][0];
        morton[1] = grid[0][1];
        morton[2] = grid[1][0];
        morton[3] = grid[1][1];
    }
    else {
        int n = N/2;
        //Divide morton-order matrix into 4 parts
        etype *morton2 = &morton[N*N/4];
        etype *morton3 = &morton[N*N/2];
        etype *morton4 = &morton[3*N*N/4];
        int i;

        //Divide 2D matrix into 4 parts
        etype **grid11, **grid12, **grid21, **grid22;
        grid11 = (etype**)malloc(n*sizeof(etype*));
        grid12 = (etype**)malloc(n*sizeof(etype*));
        grid21 = (etype**)malloc(n*sizeof(etype*));
        grid22 = (etype**)malloc(n*sizeof(etype*));
        
        for(i = 0; i < n; i++){
            grid11[i] = &grid[i][0];
            grid12[i] = &grid[i][n];
            grid21[i] = &grid[i+n][0];
            grid22[i] = &grid[i+n][n];
        }
        
        //Recursive call in each part
        grid_to_morton(grid11, morton, n);
        grid_to_morton(grid12, morton2, n);
        grid_to_morton(grid21, morton3, n);
        grid_to_morton(grid22, morton4, n);

        free(grid11); free(grid12); free(grid21); free(grid22);
    }
}
