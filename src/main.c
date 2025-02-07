// mesh funcs and mesh types
#include "mesh.h"
#include "vector.h"
// CRS lib and CG solver
#include "CRSMat_types.h"
#include "CRSmatfuncs.h"
#include "CGSolver.h"
#include "CGSolver_types.h"
// General funcs
#include "common.h"
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "string.h"
// Operator
#include "gradient.h"
// febio
#include "febiofunc.h"
#include "febio_types.h"
// ploting
#include "gnuplot.h"

// read region from a zfem file
int read_regionmask(char *path, int nelem, int npoin, int *elems, int **region_id2, int **region_idp2)
{
    int e = 0;

    // Allocate memory for region_id
    int *region_id = malloc((size_t)npoin * sizeof(*region_id));
    if (region_id == NULL)
    {
        fprintf(stderr, "Memory allocation failed for region_id.\n");
        return -1;
    }

    // Open the file
    FILE *fptr = fopen(path, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "ERROR: Cannot open file - %s.\n", path);
        free(region_id); // Free allocated memory before returning
        return -1;
    }

    /* Read all lines of the file */
    int buffer = 100;
    char *str;
    char line[buffer];
    int endcount = 0;
    int nscan, iline;
    char test[20];

    // Reading region labels
    while (1)
    {
        str = edit_endline_character(line, buffer, fptr);
        nscan = sscanf(str, "%s", test);

        if (!strcmp(test, "regions"))
        {
            str = edit_endline_character(line, buffer, fptr);
            nscan = sscanf(str, "%s", test);

            if (!strcmp(test, "1"))
            {
                for (iline = 0; iline < npoin; iline++)
                {
                    str = edit_endline_character(line, buffer, fptr);
                    nscan = sscanf(str, "%d", &region_id[iline]);

                    // Error check
                    if (nscan != 1)
                    {
                        fprintf(stderr, "ERROR: Incorrect data on line %d.\n", iline + 1);
                        free(region_id);
                        fclose(fptr);
                        return -1;
                    }
                }
                endcount += 1;
            }
        }

        if (endcount == 1)
        {
            printf("* Done Reading region mask file!\n");
            break;
        }
    }

    // Close the file
    if (fclose(fptr) == EOF)
    {
        printf("Error closing %s\n", path);
        free(region_id);
        return -1;
    }

    // Allocate memory for region_id_ele
    int *region_id_ele = malloc((size_t)nelem * sizeof(*region_id_ele));
    if (region_id_ele == NULL)
    {
        fprintf(stderr, "Memory allocation failed for region_id_ele.\n");
        free(region_id);
        return -1;
    }

    // Fill region_id_ele based on region_id
    // int points[3] = {0, 0, 0};
    for (int ele = 0; ele < nelem; ele++)
    {
        // points[0] = elems[3 * ele];
        // points[1] = elems[3 * ele + 1];
        // points[2] = elems[3 * ele + 2];

        // Assign region ID for the element based on one point (e.g., points[0] - 1)
        region_id_ele[ele] = region_id[elems[3 * ele] - 1];
        if (region_id_ele[ele] > 16)
        {
            fprintf(stderr, "there is error in reading or calculating regions mask\n");
            exit(EXIT_FAILURE);
        }
    }

    /* Return results */
    *region_id2 = region_id_ele;
    *region_idp2 = region_id;

    return e;
}

int main(void)
{
#pragma region mesh_and_data_structures
    int npoin, nelem, *elems;
    double *ptxyz;
    char path[] = {"temp/a06161.1.flds.zfem"};
    CHECK_ERROR(read_zfem(path, &npoin, &nelem, &ptxyz, &elems));
    if (ptxyz == NULL || elems == NULL)
    {
        fprintf(stderr, "Memory allocation (elems/ptxyz) failed.\n");
        return 1;
    }
    /* created required data structure for mesh */
    int Nredge = 3;
    int *esurp, *esurp_ptr, *esure, *open;
    save_esurp(npoin, nelem, elems, &esurp, &esurp_ptr, Nredge);
    if (ptxyz == NULL || elems == NULL)
    {
        fprintf(stderr, "Memory allocation (esurp/esurp_ptr) failed.\n");
        return 1;
    }
    save_esure(nelem, elems, esurp_ptr, esurp, &esure, &open, Nredge);
    if (open == NULL || esure == NULL)
    {
        fprintf(stderr, "Memory allocation (esure/open) failed.\n");
        return 1;
    }
    int opencount = 0;
    for (int i = 0; i < nelem; i++)
    {
        if (open[i] == 0)
            opencount++;
    }
    (opencount == 0) ? printf("! this is open mesh.\n") : printf("* this is close mesh.\n");

    /*calc norm of ele*/
    double *normele;
    CHECK_ERROR(save_normele(nelem, elems, ptxyz, &normele));
    // flip the normal vector to be outward:
    for (int ele = 0; ele < (3 * nelem); ele++)
        normele[ele] = -1 * normele[ele];
    /*calculate the centeroid of each element*/
    double *cen;
    CHECK_ERROR(save_centri3(nelem, elems, ptxyz, &cen));
    if (cen == NULL)
    {
        fprintf(stderr, "Memory allocation (cen) failed.\n");
        return 1;
    }
    /*calculate the center of each edge for each element*/
    double *cenedge;
    CHECK_ERROR(save_cenedgetri3(nelem, elems, ptxyz, &cenedge));
    /*calculate the normal vectors of edges*/
    double *normedge;
    CHECK_ERROR(save_normedge(nelem, ptxyz, elems, normele, &normedge));
    if (normedge == NULL)
    {
        fprintf(stderr, "Memory allocation (normedge) failed.\n");
        return 1;
    }
#pragma endregion mesh_and_data_structures
#pragma region required_mask
    // read reginal mask
    int *region_ele, *region_p;
    CHECK_ERROR(read_regionmask("temp/labels_srf.zfem", nelem, npoin, elems, &region_ele, &region_p));
    if (region_ele == NULL || region_p == NULL)
    {
        fprintf(stderr, "! ERROR in allocation memory or reading read_regionmask.\n");
        exit(EXIT_FAILURE);
    }

    // find a boundary condition
    int *cell_stat = (int *)calloc((size_t)nelem, sizeof(int));
    char bcpath[] = {"/dagon1/achitsaz/poisson/temp/a06161.1.flds.zfem.labels"};
    FILE *fptr = fopen(bcpath, "r");
    if (fptr == NULL)
    {
        fprintf(stderr, "BC file does not open properly\n");
        return 1;
    }
    int buffer = 50, nscan;
    char line[buffer];
    char *str;
    for (int ele = 0; ele < nelem; ele++)
    {
        str = edit_endline_character(line, buffer, fptr);
        nscan = sscanf(str, "%d", &cell_stat[ele]);
        if (nscan != 1)
        {
            fprintf(stderr, "! there is error in the line %d of BC file.\n", ele + 1);
            return 1;
        }
    }
    fclose(fptr);
#pragma endregion required_mask
    // creat mesh struct
    mesh *M1 = (mesh *)malloc(sizeof(mesh));
    if (M1)
    {
        *M1 = (mesh){0}; // Set all integer and pointer fields to 0 or NULL
    }
    if (M1 == NULL)
    {
        fprintf(stderr, "Memory allocation failed for M1 pointer\n");
        exit(EXIT_FAILURE);
    }
    M1->elems = elems;
    M1->npoin = npoin;
    M1->nelem = nelem;
    M1->esure = esure;
    M1->esurp = esurp;
    M1->esurp_ptr = esurp_ptr;
    M1->nredge = Nredge;
    M1->open = open;
    M1->ptxyz = ptxyz;
#ifdef DEBUG
    // coordinate
    M1->numExtraPoints = nelem + 3 * nelem;
    double *extra_ptxyz = calloc((size_t)M1->numExtraPoints * 3, sizeof(double));
    // points for normal of elements
    for (int i = 0; i < (nelem * 3); i++)
        extra_ptxyz[i] = cen[i];
    // // point for edge p1->p2
    // for (int i = 0; i < (nelem* 3); i++)
    //     extra_ptxyz[3*nelem+i] = cen[i];
    // // point for edge p2->p3
    // for (int i = 0; i < (nelem* 3); i++)
    //     extra_ptxyz[6*nelem+i] = cen[i];
    // // point for edge p3->p1
    // for (int i = 0; i < (nelem* 3); i++)
    //     extra_ptxyz[9*nelem+i] = cen[i];

    // point for edge p1->p2
    for (int ele = 0; ele < nelem; ele++)
    {
        for (int i = 0; i < 3; i++)
            extra_ptxyz[3 * nelem + 3 * ele + i] = cenedge[9 * ele + i];
    }

    // point for edge p2->p3
    for (int ele = 0; ele < nelem; ele++)
    {
        for (int i = 0; i < 3; i++)
            extra_ptxyz[6 * nelem + 3 * ele + i] = cenedge[9 * ele + 3 + i];
    }
    // point for edge p3->p1
    for (int ele = 0; ele < nelem; ele++)
    {
        for (int i = 0; i < 3; i++)
            extra_ptxyz[9 * nelem + 3 * ele + i] = cenedge[9 * ele + 6 + i];
    }

    M1->extra_ptxyz = extra_ptxyz;

    double *new_normele = calloc(((size_t)M1->npoin + (size_t)M1->numExtraPoints) * 3, sizeof(double));
    // normal vector of each element
    for (int i = 0; i < nelem; i++)
    {
        for (int j = 0; j < 3; j++)
            new_normele[3 * M1->npoin + 3 * i + j] = normele[3 * i + j];
    }
    // normal vector of edge1 p1 -> p2
    for (int i = 0; i < nelem; i++)
    {
        for (int j = 0; j < 3; j++)
            new_normele[3 * M1->npoin + 3 * nelem + 3 * i + j] = normedge[9 * i + j];
    }
    // normal vector of edge1 p2 -> p3
    for (int i = 0; i < nelem; i++)
    {
        for (int j = 0; j < 3; j++)
            new_normele[3 * M1->npoin + 6 * nelem + 3 * i + j] = normedge[9 * i + 3 + j];
    }
    // normal vector of edge1 p3 -> p1
    for (int i = 0; i < nelem; i++)
    {
        for (int j = 0; j < 3; j++)
            new_normele[3 * M1->npoin + 9 * nelem + 3 * i + j] = normedge[9 * i + 6 + j];
    }
    FunctionWithArgs prtelefield[] =
        {
            {"open", 1, nelem, open, SCA_int_VTK},
            {"BC", 1, nelem, cell_stat, SCA_int_VTK},
            {"regions", 1, nelem, region_ele, SCA_int_VTK},
        };
    size_t countele = sizeof(prtelefield) / sizeof(prtelefield[0]);
    FunctionWithArgs prtpntfield[] = {
        {"normal_vec", 3, (M1->npoin + M1->numExtraPoints), new_normele, VEC_double_VTK}};
    size_t countpnt = 1;
    CHECK_ERROR(SaveVTK("./", "checkmesh", 0, M1, tri3funcVTK, prtelefield, countele, prtpntfield, countpnt));
    free(new_normele);
    free(extra_ptxyz);
#endif
#pragma region solve_Poisson
    // define sparse matrix of coefficient
    int *row_ptr = (int *)calloc((size_t)(nelem + 1), sizeof(int));
    int max_nnz = nelem * (Nredge + 1); // At most non-zero entries per row (element+neighbours)
    double *val = (double *)malloc((size_t)max_nnz * sizeof(double));
    int *col_ind = (int *)malloc((size_t)max_nnz * sizeof(int));
    int nnz = 0;
    double *RHS = (double *)calloc((size_t)nelem, sizeof(double));
    if (row_ptr == NULL || val == NULL || col_ind == NULL || RHS == NULL)
    {
        printf("Memory allocation failed.\n");
        return 1;
    }

    // calculate the symmetrical Positive Definite Coefficient martix and Right-Hand-Sided vector
    int order_in_cells[] = {0, 1, 1, 2, 2, 0};
    for (int ele = 0; ele < nelem; ele++)
    {
        nnz++;
        int IDele = nnz - 1;
        // Known cells
        if (cell_stat[ele] > 0)
        {
            // equation for this element is u=Tb
            val[IDele] = 1;
            col_ind[IDele] = ele;
            RHS[ele] = (cell_stat[ele] == 1) ? 1000 : 0;
        }
        else
        {
            // Unknown cells
            for (int nei = 0; nei < Nredge; nei++)
            {
                int neighbor = esure[Nredge * ele + nei];
                int lp1 = elems[Nredge * ele + order_in_cells[2 * nei]] - 1;
                int lp2 = elems[Nredge * ele + order_in_cells[2 * nei + 1]] - 1;
                double dA = 0;
                dA = SQUARE((ptxyz[3 * lp1] - ptxyz[3 * lp2]));          // x-coordinate
                dA += SQUARE((ptxyz[3 * lp1 + 1] - ptxyz[3 * lp2 + 1])); // y-coordinate
                dA += SQUARE((ptxyz[3 * lp1 + 2] - ptxyz[3 * lp2 + 2])); // z-coordinate
                dA = sqrt(dA);
                double dl = 0;
                dl = SQUARE(cen[Nredge * ele] - cen[Nredge * neighbor]);          // x-coordinate
                dl += SQUARE(cen[Nredge * ele + 1] - cen[Nredge * neighbor + 1]); // y-coordinate
                dl += SQUARE(cen[Nredge * ele + 2] - cen[Nredge * neighbor + 2]); // z-coordinate
                dl = sqrt(dl);
                double coef = dA / dl;
                // Internal flux
                // contribution for diagonal element of coefficient matrix
                val[IDele] += coef;
                col_ind[IDele] = ele;
                // contribution for off-diagonal element of coefficient matrix
                if (cell_stat[neighbor] > 0)
                {
                    // between known and unknown cells
                    RHS[ele] += (cell_stat[neighbor] == 1) ? 1000 : 0;
                }
                else
                {
                    // between 2known cells
                    val[nnz] = -1 * coef;
                    col_ind[nnz] = neighbor;
                    nnz++; // cause creat new element in the coefficient matrix
                }
            }
        }
        row_ptr[ele + 1] = nnz;
    }
    row_ptr[nelem] = nnz;

    // SOLVER SECTION
    CRSMatrix A;
    A.n = nelem;
    A.nnz = nnz;
    A.row_ptr = row_ptr;
    A.col_index = col_ind;
    A.values = val;

    // unknown vector
    double *u = (double *)calloc((size_t)A.n, sizeof(double));

    // Solve using CG
    clock_t start_time, end_time;
    double cpu_time_used;
    SolverConfig config = {100000, 1e-8, true};
    solver_set_config(config);
    // precond_conjugate_gradient(&A, RHS, u);
    start_time = clock();
    conjugate_gradient(&A, RHS, u);
    end_time = clock();
    cpu_time_used = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("* CG Solver execution time : %.2f seconds\n", cpu_time_used);
    // check the result
    FunctionWithArgs prtelefield2[] =
        {
            {"BC_poisson", 1, nelem, cell_stat, SCA_int_VTK},
            {"poisson", 1, nelem, u, SCA_double_VTK},
        };
    size_t countele2 = sizeof(prtelefield2) / sizeof(prtelefield2[0]);
    FunctionWithArgs prtpntfield2[] = {0};
    size_t countpnt2 = 0;
    CHECK_ERROR(SaveVTK("./", "checkpoisson", 0, M1, tri3funcVTK, prtelefield2, countele2, prtpntfield2, countpnt2));
#pragma endregion solve_Poisson
#pragma region calc_local_basis
    // calculate the gradient scaler on each element
    double *norm_grad;
    CHECK_ERROR(gradient_ele_tri3(nelem, elems, ptxyz, esure, normedge, u, &norm_grad));
    if (norm_grad == NULL)
    {
        fprintf(stderr, "! ERROR: grad array is empty\n");
        return -1;
    }
    for (int ele = 0; ele < nelem; ele++)
    {
        double sum = 0;
        sum += SQUARE(norm_grad[3 * ele]);
        sum += SQUARE(norm_grad[3 * ele + 1]);
        sum += SQUARE(norm_grad[3 * ele + 2]);
        sum = sqrt(sum);
        if (sum != 0)
        {
            norm_grad[3 * ele] /= sum;
            norm_grad[3 * ele + 1] /= sum;
            norm_grad[3 * ele + 2] /= sum;
        }
    }

    // calculate the third vector for local coordinate system for each cells
    double *local_coord = (double *)malloc(9 * (size_t)nelem * sizeof(double));
    for (int ele = 0; ele < nelem; ele++)
    {
        double thrid_vec[3];
        double first_vec[3] = {normele[3 * ele], normele[3 * ele + 1], normele[3 * ele + 2]};
        double second_vec[3] = {norm_grad[3 * ele], norm_grad[3 * ele + 1], norm_grad[3 * ele + 2]};
        crossProduct(first_vec, second_vec, thrid_vec);
        local_coord[9 * ele] = first_vec[0];
        local_coord[9 * ele + 1] = first_vec[1];
        local_coord[9 * ele + 2] = first_vec[2];
        local_coord[9 * ele + 3] = second_vec[0];
        local_coord[9 * ele + 4] = second_vec[1];
        local_coord[9 * ele + 5] = second_vec[2];
        local_coord[9 * ele + 6] = thrid_vec[0];
        local_coord[9 * ele + 7] = thrid_vec[1];
        local_coord[9 * ele + 8] = thrid_vec[2];
    }
    // save local system in VTK format
    // coordinate
    M1->numExtraPoints = 3 * nelem;
    double *extra_ptxyz2 = calloc((size_t)M1->numExtraPoints * 3, sizeof(double));
    // points for first vec
    for (int i = 0; i < (nelem * 3); i++)
        extra_ptxyz2[i] = cen[i];
    // points for second vec
    for (int i = 0; i < (nelem * 3); i++)
        extra_ptxyz2[3 * nelem + i] = cen[i];
    // points for third vec
    for (int i = 0; i < (nelem * 3); i++)
        extra_ptxyz2[6 * nelem + i] = cen[i];

    M1->extra_ptxyz = extra_ptxyz2;

    double *new_normele2 = calloc(((size_t)M1->npoin + (size_t)M1->numExtraPoints) * 3, sizeof(double));
    // First vector of each element
    for (int i = 0; i < nelem; i++)
    {
        for (int j = 0; j < 3; j++)
            new_normele2[3 * M1->npoin + 3 * i + j] = local_coord[9 * i + j];
    }
    // Second vector
    for (int i = 0; i < nelem; i++)
    {
        for (int j = 0; j < 3; j++)
            new_normele2[3 * M1->npoin + 3 * nelem + 3 * i + j] = local_coord[9 * i + 3 + j];
    }
    // Third vector
    for (int i = 0; i < nelem; i++)
    {
        for (int j = 0; j < 3; j++)
            new_normele2[3 * M1->npoin + 6 * nelem + 3 * i + j] = local_coord[9 * i + 6 + j];
    }
    // save local sys in the VTK format
    FunctionWithArgs prtelefield3[] =
        {
            {"norm_grad_u", 3, nelem, norm_grad, VEC_double_VTK},
        };
    size_t countele3 = sizeof(prtelefield3) / sizeof(prtelefield3[0]);
    FunctionWithArgs prtpntfield3[] = {
        {"local_normal_vec", 3, (M1->npoin + M1->numExtraPoints), new_normele2, VEC_double_VTK}};
    size_t countpnt3 = 1;
    CHECK_ERROR(SaveVTK("./", "checklocalsys", 0, M1, tri3funcVTK, prtelefield3, countele3, prtpntfield3, countpnt3));
    free(extra_ptxyz2);
    free(new_normele2);
#pragma endregion local_basis

    // read stress tensor from log file
    double *st;
    readfebiolog("temp/pres_0.log", nelem, &st, end_first_step);
    if (st == NULL)
    {
        fprintf(stderr, "there is problem in reading stress tensor\n");
        exit(EXIT_FAILURE);
    }

#pragma region extract_in-plane_tensor
    // analysis tensor in tangential plane
    double *shear_st = (double *)malloc((size_t)nelem * 4 * sizeof(double));
    double *shear_evals_max = (double *)malloc((size_t)nelem * sizeof(double));
    double *shear_evals_min = (double *)malloc((size_t)nelem * sizeof(double));
    double *shear_evects_3D_min = (double *)malloc((size_t)nelem * 3 * sizeof(double));
    double *shear_evects_3D_max = (double *)malloc((size_t)nelem * 3 * sizeof(double));
    for (int ele = 0; ele < nelem; ele++)
    {
        // extract stress and normal, t1 and t2 tangantial vectors from arraies
        double stress[3][3], normal[3], t1[3], t2[3];
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
                stress[i][j] = st[9 * ele + 3 * i + j];
        }
        for (int i = 0; i < 3; i++)
            normal[i] = local_coord[9 * ele + i];
        for (int i = 0; i < 3; i++)
            t1[i] = local_coord[9 * ele + 3 + i];
        for (int i = 0; i < 3; i++)
            t2[i] = local_coord[9 * ele + 6 + i];

        // Compute the projection matrix
        double P[3][3];
        compute_projection_matrix(normal, P);

        // Compute the in-plane stress tensor
        double in_plane[3][3];
        compute_in_plane_stress(stress, P, in_plane);

        // Extract the 2x2 tangential stress tensor
        double tensor_2x2[2][2];
        extract_2x2_tensor(in_plane, t1, t2, tensor_2x2);

        // save in-plane stress tensor
        for (int i = 0; i < 2; i++)
        {
            for (int j = 0; j < 2; j++)
                shear_st[4 * ele + 2 * i + j] = tensor_2x2[i][j];
        }

        // Find eigenvalues
        double lambda1, lambda2;
        find_eigenvalues_2x2(tensor_2x2, &lambda1, &lambda2);

        // Find eigenvectors
        double v2D1[2], v2D2[2];
        find_eigenvector_2x2(tensor_2x2, lambda1, v2D1);
        find_eigenvector_2x2(tensor_2x2, lambda2, v2D2);

        // Reconstruct eigenvectors in 3D space
        double v3D1[3], v3D2[3];
        for (int i = 0; i < 3; i++)
            v3D1[i] = t1[i] * v2D1[0] + t2[i] * v2D1[1];
        for (int i = 0; i < 3; i++)
            v3D2[i] = t1[i] * v2D2[0] + t2[i] * v2D2[1];

        // save eigen values and reconstructed eigen vectors in 3D
        if (fabs(lambda1) > fabs(lambda2))
        {
            shear_evals_max[ele] = lambda1;
            shear_evals_min[ele] = lambda2;
            for (int i = 0; i < 3; i++)
                shear_evects_3D_max[3 * ele + i] = v3D1[i];
            for (int i = 0; i < 3; i++)
                shear_evects_3D_min[3 * ele + i] = v3D2[i];
        }
        else
        {
            shear_evals_max[ele] = lambda2;
            shear_evals_min[ele] = lambda1;
            for (int i = 0; i < 3; i++)
                shear_evects_3D_max[3 * ele + i] = v3D2[i];
            for (int i = 0; i < 3; i++)
                shear_evects_3D_min[3 * ele + i] = v3D1[i];
        }
    }

    // // print in-plane tensor
    // int num_prt = 5;
    // for (int ele =0;ele<nelem;ele++){
    //     if (num_prt<0) break;
    //     if (region_ele[ele]==16 ){
    //         printf("in-plane stress ele : %d \n",ele);
    //         for (int i=0;i<2;i++) {
    //             for (int j=0;j<2;j++)
    //                 printf("%lf ",shear_st[4*ele+2*i+j]);
    //             printf("\n");
    //         }
    //         num_prt--;
    //     }
    // }

    // check the symmetricity of in-plan tensor
    int num_sys = 0;
    int num_non_sys = 0;
    for (int ele = 0; ele < nelem; ele++)
    {
        if (region_ele[ele] == 16 || region_ele[ele] == 8 || region_ele[ele] == 4)
        {
            if (shear_st[4 * ele + 1] - shear_st[4 * ele + 2] > 1)
            {
                // fprintf(stderr,"! in-plan tensor of ele = %d (region: %d) is not symmetric.\n",ele,region_ele[ele]);
                // exit(EXIT_FAILURE);
                num_non_sys++;
            }
            num_sys++;
        }
    }
    (num_non_sys == 0) ? printf("* in-plane stress tensor is symmetric.\n") : printf("! nonsymetric : %d symetric: %d element.\n", num_non_sys, num_sys);
#pragma endregion extract_in - plane_tensor

    // disturbution of eval_min/eval_max
    double *eval_ratio_aneu = (double *)malloc((size_t)nelem * sizeof(double));
    double *eval_ratio = (double *)malloc((size_t)nelem * sizeof(double));
    int num_aneu = 0;
    for (int ele = 0; ele < nelem; ele++)
    {
        eval_ratio[ele] = (shear_evals_max[ele] > 1) ? fabs(shear_evals_min[ele] / shear_evals_max[ele]) : 0;
        if (region_ele[ele] == 16 || region_ele[ele] == 8 || region_ele[ele] == 4)
        {
            eval_ratio_aneu[num_aneu] = fabs(shear_evals_min[ele] / shear_evals_max[ele]);
            num_aneu++;
        }
    }
    // creat Histogram to show disturbution of eval_ratio
    int num_values = num_aneu; // Number of values to analyze
    double max_value = 1;      // Maximum range value
    double min_value = 0;      // Min range value
    int num_bins = 10;         // Number of bins
    int *bins;

    // calculate the disturbution in each bin
    int max_bin_count = compute_distribution_double(eval_ratio_aneu, num_values, &bins, num_bins, max_value, min_value);
    // save histogeram
    plot_histogram(bins, num_bins, max_value, num_values, "Eval_ratio.dat", "Histogram of Eval ratio on aneurysm region", "ratio", "frequency", "Eval_ratio_histogram.png", "min/max", max_bin_count);
    free(bins);
    // find uni or bi-directional region
    int *sdir;
    unibimask(nelem, shear_evals_max, shear_evals_min, 0.2, 10000, &sdir);
    if (sdir == NULL)
    {
        fprintf(stderr, "! there is problem in allocate memory for sdir.\n");
        exit(EXIT_FAILURE);
    }
    // Calculate the Von Mises Stress
    double *von_mises = (double *)malloc((size_t)nelem * sizeof(double));
    for (int ele = 0; ele < nelem; ele++)
        von_mises[ele] = sqrt(SQUARE(shear_evals_max[ele]) + SQUARE(shear_evals_min[ele]));

    // calculate mask according sign of eigen values
    int *eigen_class = (int *)calloc((size_t)nelem, sizeof(int));
    for (int ele = 0; ele < nelem; ele++)
    {

        if (sdir[ele] == 2)
        {
            if (shear_evals_max[ele] > 0 && shear_evals_min[ele] > 0)
                eigen_class[ele] = 1;
            if (shear_evals_max[ele] > 0 && shear_evals_min[ele] < 0)
                eigen_class[ele] = 2;
            if (shear_evals_max[ele] < 0 && shear_evals_min[ele] > 0)
                eigen_class[ele] = 3;
            if (shear_evals_max[ele] < 0 && shear_evals_min[ele] < 0)
                eigen_class[ele] = 4;
        }
        if (sdir[ele] == 1)
        {
            if (shear_evals_max[ele] > 0)
                eigen_class[ele] = 5;
            if (shear_evals_max[ele] < 0)
                eigen_class[ele] = 6;
        }
    }
    // find disturbution in each class just for aneurysm
    int *eigen_class_disturb = (int *)calloc((size_t)nelem, sizeof(int));
    num_aneu = 0;
    for (int ele = 0; ele < nelem; ele++)
    {
        if (region_ele[ele] == 16 || region_ele[ele] == 8 || region_ele[ele] == 4)
        {
            eigen_class_disturb[num_aneu] = eigen_class[ele];
            num_aneu++;
        }
    }

    // creat Histogram to show disturbution of eval_ratio
    int num_values2 = num_aneu; // Number of values to analyze
    int max_value2 = 6;         // Maximum range value
    int min_value2 = 1;         // Min range value
    int num_bins2 = 6;          // Number of bins
    int *bins2;

    // calculate the disturbution in each bin
    int max_bin_count2 = compute_distribution_int(eigen_class_disturb, num_values2, &bins2, num_bins2, max_value2, min_value2);
    // save histogeram
    plot_histogram(bins2, num_bins2, max_value2, num_values2, "eigen_class_disturb.dat", "Histogram of eigen class on aneurysm region", "classes", "frequency", "Eigen_classes_histogram.png", "classes", max_bin_count2);
    free(bins2);

    // calculate the rotation angle for max principal stress 
    double *rotation_max_t1 = (double *)calloc ((size_t)nelem,sizeof(double));
    if (!rotation_max_t1){
        fprintf(stderr,"! Failure in Memory Allication rotation_max_t1 .\n");
        exit(EXIT_FAILURE);
    }
    for (int ele =0;ele<nelem;ele++){
        if (shear_evals_max[ele] < 1) continue;
        double evec_c[3] = {shear_evects_3D_max[3 * ele], shear_evects_3D_max[3 * ele + 1], shear_evects_3D_max[3 * ele + 2]};
        double t1[3]={local_coord[9*ele+3],local_coord[9*ele+4],local_coord[9*ele+5]};
        double dot_product = evec_c[0] * t1[0] + evec_c[1] * t1[1] + evec_c[2] * t1[2];
        double degree = acos(fabs(dot_product)) * 180.0 / PI;
        // double degree_min = 90;
        // if (degree > 45)
        // {
        //     double cen_landa_min = shear_evals_min[ele];
        //     double cen_landa_max = shear_evals_max[ele];
        //     double ratio = cen_landa_min / cen_landa_max;
        //     if (ratio > 0.1)
        //     {
        //         double evec_c_min[3] = {shear_evects_3D_min[3 * ele], shear_evects_3D_min[3 * ele + 1], shear_evects_3D_min[3 * ele + 2]};
        //         double dot_product_min = evec_c_min[0] * t1[0] + evec_c_min[1] * t1[1] + evec_c_min[2] * t1[2];
        //         degree_min = acos(fabs(dot_product_min)) * 180.0 / PI;
        //     }
        // }
        // degree = MIN(degree, degree_min);
        rotation_max_t1[ele]= MAX(degree, rotation_max_t1[ele]);
    }
    // calculate relative angels of eigen vectors
    double *evec_max_angles = (double *)calloc((size_t)nelem * 3, sizeof(double));
    if (!evec_max_angles)
    {
        fprintf(stderr, "! Failure in Memory Allocation evec_max_angles.\n");
        exit(EXIT_FAILURE);
    }
    for (int ele = 0; ele < nelem; ele++)
    {
        double evec_c[3] = {shear_evects_3D_max[3 * ele], shear_evects_3D_max[3 * ele + 1], shear_evects_3D_max[3 * ele + 2]};
        if (shear_evals_max[ele] < 1)
            continue;
        for (int nei = 0; nei < 3; nei++)
        {
            int ele_nei = esure[3 * ele + nei];
            double evec_nei[3] = {shear_evects_3D_max[3 * ele_nei], shear_evects_3D_max[3 * ele_nei + 1], shear_evects_3D_max[3 * ele_nei + 2]};
            double dot_product = evec_c[0] * evec_nei[0] + evec_c[1] * evec_nei[1] + evec_c[2] * evec_nei[2];
            double degree = acos(fabs(dot_product)) * 180.0 / PI;
            double degree_min = 90;
            if (degree > 45)
            {
                double nei_landa_min = shear_evals_min[ele_nei];
                double cen_landa_max = shear_evals_max[ele];
                double ratio = nei_landa_min / cen_landa_max;
                if (ratio > 0.5)
                {
                    double evec_nei_min[3] = {shear_evects_3D_min[3 * ele_nei], shear_evects_3D_min[3 * ele_nei + 1], shear_evects_3D_min[3 * ele_nei + 2]};
                    double dot_product_min = evec_c[0] * evec_nei_min[0] + evec_c[1] * evec_nei_min[1] + evec_c[2] * evec_nei_min[2];
                    degree_min = acos(fabs(dot_product_min)) * 180.0 / PI;
                }
            }
            degree = MIN(degree, degree_min);
            evec_max_angles[ele] = MAX(degree, evec_max_angles[ele]);
        }
    }
    // save tensor analysis in VTK format
    M1->numExtraPoints = (4) * nelem;
    double *extra_ptxyz3 = calloc((size_t)M1->numExtraPoints * 3, sizeof(double));
    // points for first eigvec
    for (int i = 0; i < (nelem * 3); i++)
        extra_ptxyz3[i] = cen[i];
    // points for second eigvec
    for (int i = 0; i < (nelem * 3); i++)
        extra_ptxyz3[3 * nelem + i] = cen[i];
    // points for first orivec
    for (int i = 0; i < (nelem * 3); i++)
        extra_ptxyz3[2 * 3 * nelem + i] = cen[i];
    // points for second orivec
    for (int i = 0; i < (nelem * 3); i++)
        extra_ptxyz3[3 * 3 * nelem + i] = cen[i];

    M1->extra_ptxyz = extra_ptxyz3;

    double *new_normele3 = calloc(((size_t)M1->npoin + (size_t)M1->numExtraPoints) * 3, sizeof(double));
    // First vector of each element
    for (int i = 0; i < nelem; i++)
    {
        for (int j = 0; j < 3; j++)
            new_normele3[3 * M1->npoin + 3 * i + j] = shear_evals_max[i] * shear_evects_3D_max[3 * i + j];
    }
    // Second vector
    for (int i = 0; i < nelem; i++)
    {
        for (int j = 0; j < 3; j++)
            new_normele3[3 * M1->npoin + 3 * nelem + 3 * i + j] = shear_evals_min[i] * shear_evects_3D_min[3 * i + j];
    }

    double *ori = calloc(((size_t)M1->npoin + (size_t)M1->numExtraPoints) * 3, sizeof(double));
    double *ori_sign = calloc((size_t)M1->npoin + (size_t)M1->numExtraPoints, sizeof(double));
    // First vector
    for (int i = 0; i < nelem; i++)
    {
        ori_sign [M1->npoin + i] = (shear_evals_max[i]>0) ? 1 : -1 ;
        for (int j = 0; j < 3; j++)
            ori[3 * M1->npoin + 3 * i + j] = shear_evals_max[i] * shear_evects_3D_max[3 * i + j];
    }
    // Second vector
    for (int i = 0; i < nelem; i++)
    {
        ori_sign [M1->npoin + nelem + i] = (shear_evals_min[i]>0) ? 1 : -1 ;
        for (int j = 0; j < 3; j++)
            ori[3 * M1->npoin + 3 * nelem + 3 * i + j] = shear_evals_min[i] * shear_evects_3D_min[3 * i + j];
    }
    // Reversed first vector
    for (int i = 0; i < nelem; i++)
    {
        ori_sign [M1->npoin + 2 * nelem + i] = (shear_evals_max[i]>0) ? 1 : -1 ;
        for (int j = 0; j < 3; j++)
            ori[3 * M1->npoin + 6 * nelem + 3 * i + j] = -1 * shear_evals_max[i] * shear_evects_3D_max[3 * i + j];
    }
    // Reversed second vector
    for (int i = 0; i < nelem; i++)
    {
        ori_sign [M1->npoin + 3 * nelem + i] = (shear_evals_min[i]>0) ? 1 : -1 ;
        for (int j = 0; j < 3; j++)
            ori[3 * M1->npoin + 9 * nelem + 3 * i + j] = -1 * shear_evals_min[i] * shear_evects_3D_min[3 * i + j];
    }

    FunctionWithArgs prtelefield4[] =
        {
            {"uni/bi_region", 1, nelem, sdir, SCA_int_VTK},
            {"eigen_class", 1, nelem, eigen_class, SCA_int_VTK},
            {"EValue_max", 1, nelem, shear_evals_max, SCA_double_VTK},
            {"Shear_Von_Mises_Stress", 1, nelem, von_mises, SCA_double_VTK},
            {"Eval_ratio", 1, nelem, eval_ratio, SCA_double_VTK},
            {"evec_max_angles", 1, nelem, evec_max_angles, SCA_double_VTK},
            {"rotation_max_t1", 1, nelem, rotation_max_t1, SCA_double_VTK}};            
    size_t countele4 = sizeof(prtelefield4) / sizeof(prtelefield4[0]);
    FunctionWithArgs prtpntfield4[] = {
        {"ori_sign", 1, (M1->npoin + M1->numExtraPoints), ori_sign, SCA_double_VTK},
        {"eval_dot_evec", 3, (M1->npoin + M1->numExtraPoints), new_normele3, VEC_double_VTK},
        {"ori", 3, (M1->npoin + M1->numExtraPoints), ori, VEC_double_VTK}};
    size_t countpnt4 = 3;
    CHECK_ERROR(SaveVTK("./", "tensor_analysis", 0, M1, tri3funcVTK, prtelefield4, countele4, prtpntfield4, countpnt4));
    free(extra_ptxyz3);
    free(new_normele3);
#pragma region free_dynamics_alloc
    // free dynamics arraies
    free(elems);
    free(ptxyz);
    free(open);
    free(esure);
    free(esurp);
    free(esurp_ptr);
    free(row_ptr);
    free(val);
    free(col_ind);
    free(RHS);
    free(cen);
    free(M1);
    free(cell_stat);
    free(u);
    free(cenedge);
    free(normele);
    free(normedge);
    free(norm_grad);
    free(local_coord);
    free(st);
    free(shear_evals_max);
    free(shear_evals_min);
    free(shear_evects_3D_max);
    free(shear_evects_3D_min);
    free(sdir);
    free(eval_ratio);
    free(eval_ratio_aneu);
    free(shear_st);
    free(region_ele);
    free(region_p);
    free(von_mises);
    free(eigen_class);
    free(eigen_class_disturb);
    free(evec_max_angles);
    free(ori);
    free(ori_sign);
    free(rotation_max_t1);
#pragma endregion free_dynamics_alloc
    return 0; // success signal
}