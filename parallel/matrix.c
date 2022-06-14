#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "matrix.h"

#define CUT_OFF 256 //Early termination point, varying from system to system

void mat_mul(etype *A, etype *B, etype *res, int N){
    if(N <= CUT_OFF){
        etype **gridA, **gridB, **gridC;
        mem_alloc_grid(&gridA, N);
        mem_alloc_grid(&gridB, N);
        mem_alloc_grid(&gridC, N);
        int i,j,k;
        for(i = 0; i < N; i++)
            for(k = 0; k < N; k++)
                gridC[i][k] = 0;

        morton_to_grid(gridA, A, N);
        morton_to_grid(gridB, B, N);
        
        // #pragma omp taskloop num_tasks(2)
        for(i = 0; i < N; i++)
            for(k = 0; k < N; k++)
                for(j = 0; j < N; j++)
                    gridC[i][j] += gridA[i][k] * gridB[k][j];
        
        grid_to_morton(gridC, res, N);
        del_grid_matrix(gridA, N);
        del_grid_matrix(gridB, N);
        del_grid_matrix(gridC, N);
    }
    else {
        etype *A11, *A12, *A21, *A22;
        etype *B11, *B12, *B21, *B22;
        etype *C11, *C12, *C21, *C22;
        etype *M1, *M2, *M3, *M4, *M5, *M6, *M7;
        etype *T1, *T2, *T3;

        //Auxiliary variables
        int n = N/2;
        int N1 = N*N/4;
        int N2 = N*N/2;
        int N3 = 3*N*N/4;
        
        //Allocate memory to auxiliary matrices
        mem_alloc(&M1, n); mem_alloc(&M2, n); mem_alloc(&M3, n); mem_alloc(&M4, n); 
        mem_alloc(&M5, n); mem_alloc(&M6, n); mem_alloc(&M7, n);
        mem_alloc(&T1, n); mem_alloc(&T2, n); mem_alloc(&T3, n);
        
         /* STAGE 1: Assign sub matrices */
        A11 = A; A12 = &A[N1]; A21 = &A[N2]; A22 = &A[N3];
        B11 = B; B12 = &B[N1]; B21 = &B[N2]; B22 = &B[N3];
        C11 = res; C12 = &res[N1]; C21 = &res[N2]; C22 = &res[N3];
        
        /* STAGE 2: Calculate auxilirary matrices */

        #pragma omp task untied
        {
            //Calculate M1 = (A11 + A22) * (B11 + B22)
            mat_add(A11, A22, T1, n);
            mat_add(B11, B22, M1, n);
            mat_mul(T1, M1, M1, n);
        }
        #pragma omp task untied
        {
            //Calculate M2 = (A21 + A22) * B11
            mat_add(A21, A22, M2, n);
            mat_mul(M2, B11, M2, n);
        }
        #pragma omp task untied
        {
            //Calculate M3 = A11 * (B12 + B22)
            mat_sub(B12, B22, M3, n);
            mat_mul(A11, M3, M3, n);
        }
        #pragma omp task untied
        {
            //Calculate M4 = A22 * (B21 - B11)
            mat_sub(B21, B11, M4, n);
            mat_mul(A22, M4, M4, n);
        }
        #pragma omp task untied
        {
            //Calculate M5 = (A11 + A12) * B22
            mat_add(A11, A12, M5, n);
            mat_mul(M5, B22, M5, n);
        }
        #pragma omp task untied
        {
            //Calculate M6 = (A21 - A11) * (B11 + B12)
            mat_sub(A21, A11, T2, n);
            mat_add(B11, B12, M6, n);
            mat_mul(T2, M6, M6, n);
        }
        #pragma omp task untied
        {
            //Calculate M7 = (A12 - A22) * (B21 + B22)
            mat_sub(A12, A22, T3, n);
            mat_add(B21, B22, M7, n);
            mat_mul(T3, M7, M7, n);
        }
        #pragma omp taskwait

        /* STAGE 3: Aggregate final result */
        int i;
        for(i = 0; i < N1; i++){
            //C11 = M1 + M4 - M5 + M7
            C11[i] = M1[i] + M4[i] - M5[i] + M7[i];
            //C12 = M3 + M5
            C12[i] = M3[i] + M5[i];
            //C21 = M2 + M4
            C21[i] = M2[i] + M4[i];
            //C22 = M1 - M2 + M3 + M6
            C22[i] = M1[i] - M2[i] + M3[i] + M6[i];
        }

        //Delete auxiliary matrices
        free(M1); free(M2); free(M3); free(M4); free(M5); free(M6); free(M7);
        free(T1); free(T2); free(T3);
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

    printf("Max diff in an element: %3.15lf\n", max_diff);
    del_grid_matrix(C_test, N);
}

void init_morton_matrix(etype **mat, int N){
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

void mem_alloc_grid(etype ***mat, int N){
    int i;
    *mat = (etype**)malloc(N * sizeof(etype*));
    for(i = 0; i < N; i++)
        (*mat)[i] = (etype*)malloc(N * sizeof(etype));
}

void del_grid_matrix(etype **mat, int N){
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
            printf("%3.5lf ",mat[i][j]);
        printf("\n");
    }
    printf("\n");
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