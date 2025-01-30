#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "vector_types.h"


int compare_int_min(void *a, void *b)
{
	return (*(int *)a < *(int *)b);
}
int compare_int_max(void *a, void *b)
{
	return (*(int *)a > *(int *)b);
}
int compare_double_min(void *a, void *b)
{
	return (*(double *)a < *(double *)b);
}
int compare_double_max(void *a, void *b)
{
	return (*(double *)a > *(double *)b);
}
void *find_extreme(void *array, size_t element_size, size_t num_elements, compare_func comp)
{
	void *extreme = array;

	for (size_t i = 1; i < num_elements; ++i)
	{
		void *current_element = (char *)array + i * element_size;
		if (comp(current_element, extreme))
		{
			extreme = current_element;
		}
	}

	return extreme;
}
// check the start ID of elements in the mesh file
int checkEIDS(int *elems)
{
	size_t int_size = sizeof(*elems) / sizeof(elems[0]);
	// Find min and max for int array
	int *int_min = (int *)find_extreme(elems, sizeof(int), int_size, compare_int_min);
	// printf("--> ID of elements start from %d!\n",*int_min);
	return *int_min;
}
// Function to calculate the cross product of two vectors
void crossProduct(double v1[3], double v2[3], double result[3]) {
    result[0] = v1[1] * v2[2] - v1[2] * v2[1];
    result[1] = v1[2] * v2[0] - v1[0] * v2[2];
    result[2] = v1[0] * v2[1] - v1[1] * v2[0];
}
// Function to normalize a vector
void normalize(double vector[3]) {
    double magnitude = sqrt(vector[0] * vector[0] + vector[1] * vector[1] + vector[2] * vector[2]);
    if (magnitude > 0.0) {
        vector[0] /= magnitude;
        vector[1] /= magnitude;
        vector[2] /= magnitude;
    }
}
// Function to compute the projection matrix for the tangential plane
void compute_projection_matrix(double n[3], double P[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            P[i][j] = (i == j ? 1.0 : 0.0) - n[i] * n[j];
        }
    }
}
// Function to multiply two 3x3 matrices
void matrix_multiply_3x3(double A[3][3], double B[3][3], double result[3][3])
{
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            result[i][j] = 0.0;
            for (int k = 0; k < 3; k++)
            {
                result[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}
// Function to compute the in-plane stress tensor
void compute_in_plane_stress(double stress[3][3], double P[3][3], double in_plane[3][3])
{
    double temp[3][3];
    matrix_multiply_3x3(P, stress, temp);
    matrix_multiply_3x3(temp, P, in_plane);
}
// Function to extract the 2x2 tensor in the tangential plane
void extract_2x2_tensor(double in_plane[3][3], double t1[3], double t2[3], double tensor_2x2[2][2])
{
    double t1s[3], t2s[3];
    for (int i = 0; i < 3; i++)
    {
        t1s[i] = 0;
        for (int j = 0; j < 3; j++)
            t1s[i] += t1[j] * in_plane[j][i];
    }
    for (int i = 0; i < 3; i++)
    {
        t2s[i] = 0;
        for (int j = 0; j < 3; j++)
            t2s[i] += t2[j] * in_plane[j][i];
    }
    tensor_2x2[0][0] = t1[0] * t1s[0] + t1[1] * t1s[1] + t1[2] * t1s[2];
    tensor_2x2[0][1] = t2[0] * t1s[0] + t2[1] * t1s[1] + t2[2] * t1s[2];
    tensor_2x2[1][0] = t1[0] * t2s[0] + t1[1] * t2s[1] + t1[2] * t2s[2];
    tensor_2x2[1][1] = t2[0] * t2s[0] + t2[1] * t2s[1] + t2[2] * t2s[2];
}
// Function to find eigenvalues of a 2x2 matrix
int find_eigenvalues_2x2(double tensor[2][2], double *lambda1, double *lambda2)
{
    double trace = tensor[0][0] + tensor[1][1];
    double determinant = tensor[0][0] * tensor[1][1] - tensor[0][1] * tensor[1][0];
    if ((trace * trace - 4 * determinant) < 0)
    {
        fprintf(stderr, "ERROR: eigenvalues of the extracted matrix 2*2 is complex number!\n ");
        exit(EXIT_FAILURE);
    }
    double discriminant = sqrt(trace * trace - 4 * determinant);
    *lambda1 = (trace + discriminant) / 2.0;
    *lambda2 = (trace - discriminant) / 2.0;
    return 0; //    successsful signal
}
// Function to find an eigenvector of a 2x2 matrix for a given eigenvalue
void find_eigenvector_2x2(double tensor[2][2], double lambda, double v[2])
{
    double TOL = 1e-8;
	if (fabs(tensor[0][1]) > TOL)
    {
        v[0] = lambda - tensor[1][1];
        v[1] = tensor[0][1];
    }
    else if (fabs(tensor[1][0]) > TOL)
    {
        v[0] = tensor[1][0];
        v[1] = lambda - tensor[0][0];
    }
    else
    {
        // Diagonal matrix case: Eigenvector can be [1, 0] or [0, 1]
        if (fabs(lambda - tensor[0][0]) < TOL)
        {
            v[0] = 1.0;
            v[1] = 0.0;
        }
        else
        {
            v[0] = 0.0;
            v[1] = 1.0;
        }
    }
    double magnitude = sqrt(v[0] * v[0] + v[1] * v[1]);
    if (magnitude < TOL)
    {
        fprintf(stderr, "The magnitude of eigenvector is zero!\n ");
        exit(EXIT_FAILURE);
    }
    for (int i = 0; i < 2; i++)
    {
        v[i] /= magnitude;
    }
}
