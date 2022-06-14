typedef double etype;

//Convert matrix from Morton-ordering to 2D and vice versa
void morton_to_grid(etype **grid, etype *morton, int N);
void grid_to_morton(etype **grid, etype *morton, int N);

//Validate multiplication result
void validate_result(etype **A, etype **B, etype **C, int N);

// For initialize matricies
void mem_alloc(etype **mat, int N);
void mem_alloc_grid(etype ***mat, int N);
void init_morton_matrix(etype **arr, int N);

//Free memory for 2D matrices
void del_grid_matrix(etype **mat, int N);

//random number
etype rand_etype(etype min, etype max);

//Print 2D matricies
void print_grid_matrix(etype **mat, int N);

//Matrix-matric operations, take parameters as Morton-order matricies
void mat_mul(etype *A, etype *B, etype *res, int N, int step,int cut_off); //res = A * B
void mat_add(etype *A, etype *B, etype *res, int N); //res = A + B
void mat_sub(etype *A, etype *B, etype *res, int N); //res = A - B