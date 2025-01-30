#ifndef VECTOR_H
#define VECTOR_H
    #include <stdlib.h> 
    #include "vector_types.h"
    int checkEIDS(int *elems);
    int compare_int_min(void *a, void *b);
    int compare_int_max(void *a, void *b);
    int compare_double_min(void *a, void *b);
    int compare_double_max(void *a, void *b);
    void *find_extreme(void *array, size_t element_size, size_t num_elements, compare_func comp);
    void crossProduct(double v1[3], double v2[3], double result[3]);
    void normalize(double vector[3]);
    void compute_projection_matrix(double n[3], double P[3][3]);
    void matrix_multiply_3x3(double A[3][3], double B[3][3], double result[3][3]);
    void compute_in_plane_stress(double stress[3][3], double P[3][3], double in_plane[3][3]);
    void extract_2x2_tensor(double in_plane[3][3], double t1[3], double t2[3], double tensor_2x2[2][2]);
    int find_eigenvalues_2x2(double tensor[2][2], double *lambda1, double *lambda2);
    void find_eigenvector_2x2(double tensor[2][2], double lambda, double v[2]);
#endif
