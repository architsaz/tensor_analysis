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
// Operator
#include "gradient.h"
// febio
#include "febiofunc.h"
#include "febio_types.h"
// ploting 
#include "gnuplot.h"

int main (void){
    int npoin, nelem, *elems;
    double *ptxyz;
    char path [] = {"temp/a06161.1.flds.zfem"};
    CHECK_ERROR(read_zfem(path,&npoin,&nelem,&ptxyz,&elems));
    if (ptxyz == NULL || elems == NULL)
    {
        fprintf(stderr,"Memory allocation (elems/ptxyz) failed.\n");
        return 1;
    }
    // created required data structure for mesh 
    int Nredge = 3;
    int *esurp,*esurp_ptr,*esure,*open;
    save_esurp(npoin,nelem,elems,&esurp,&esurp_ptr,Nredge);
    if (ptxyz == NULL || elems == NULL)
    {
        fprintf(stderr,"Memory allocation (esurp/esurp_ptr) failed.\n");
        return 1;
    }
    save_esure(nelem,elems,esurp_ptr,esurp,&esure,&open,Nredge);
    if (open == NULL || esure == NULL)
    {
        fprintf(stderr,"Memory allocation (esure/open) failed.\n");
        return 1;
    }
    int opencount = 0;
    for (int i=0;i<nelem;i++){
        if (open[i]==0) opencount++;
    }
    (opencount==0) ? printf("! this is open mesh.\n") : printf("* this is close mesh.\n");
    // calc norm of ele
    double *normele;
    CHECK_ERROR(save_normele(nelem,elems,ptxyz,&normele));
    // flip the normal vector to be outward:
    for (int ele = 0; ele < (3 * nelem); ele++)
        normele[ele] = -1 * normele[ele];
    // calculate the centeroid of each element
    double *cen;
    CHECK_ERROR(save_centri3(nelem,elems,ptxyz,&cen));
    if (cen == NULL){
        fprintf(stderr,"Memory allocation (cen) failed.\n");
        return 1;        
    }
    // calculate the center of each edge for each element
    double *cenedge;
    CHECK_ERROR(save_cenedgetri3(nelem,elems,ptxyz,&cenedge));
    // calculate the normal vectors of edges
    double *normedge;
    CHECK_ERROR(save_normedge(nelem,ptxyz,elems,normele,&normedge)); 
    if (normedge == NULL){
        fprintf(stderr,"Memory allocation (normedge) failed.\n");
        return 1;        
    }
    // find a boundary condition
    int *cell_stat = (int *)calloc((size_t)nelem,sizeof(int));
    char bcpath [] = {"/dagon1/achitsaz/poisson/temp/a06161.1.flds.zfem.labels"};
    FILE *fptr = fopen(bcpath,"r");
    if (fptr == NULL){
        fprintf(stderr,"BC file does not open properly\n");
        return 1;
    }
    int buffer = 50, nscan;
    char line[buffer];
    char *str;
    for (int ele =0;ele<nelem;ele++){
        str=edit_endline_character(line,buffer,fptr);
        nscan = sscanf(str, "%d", &cell_stat[ele]);
        if (nscan!=1){
            fprintf(stderr,"! there is error in the line %d of BC file.\n",ele+1);
            return 1;
        }
    }
    fclose(fptr);
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
    M1->elems =elems;M1->npoin =npoin;M1->nelem=nelem;M1->esure=esure;M1->esurp=esurp;M1->esurp_ptr=esurp_ptr;M1->nredge=Nredge;M1->open=open;M1->ptxyz=ptxyz;
    #ifdef DEBUG
        // coordinate 
            M1->numExtraPoints = nelem+3*nelem;
            double *extra_ptxyz = calloc((size_t)M1->numExtraPoints * 3, sizeof(double));
            // points for normal of elements
            for (int i = 0; i < (nelem* 3); i++)
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
            for (int ele = 0; ele < nelem; ele++){
                for (int i=0;i<3;i++)
                extra_ptxyz[3*nelem+3*ele+i] = cenedge[9*ele+i];
            }
                
            // point for edge p2->p3    
            for (int ele = 0; ele < nelem; ele++){
                for (int i=0;i<3;i++)
                extra_ptxyz[6*nelem+3*ele+i] = cenedge[9*ele+3+i];
            }
            // point for edge p3->p1    
            for (int ele = 0; ele < nelem; ele++){
                for (int i=0;i<3;i++)
                extra_ptxyz[9*nelem+3*ele+i] = cenedge[9*ele+6+i];
            }

            M1->extra_ptxyz = extra_ptxyz;

            double *new_normele = calloc(((size_t)M1->npoin + (size_t)M1->numExtraPoints) * 3, sizeof(double));
            // normal vector of each element
            for (int i = 0; i < nelem; i++)
            {
                for (int j=0;j<3;j++)
                new_normele[3*M1->npoin + 3*i+j] = normele[3*i+j];
            }
            // normal vector of edge1 p1 -> p2
            for (int i = 0; i < nelem; i++)
            {
                for (int j=0;j<3;j++)
                new_normele[3*M1->npoin +3*nelem+ 3*i+j] = normedge[9*i+j];
            }
            // normal vector of edge1 p2 -> p3
            for (int i = 0; i < nelem; i++)
            {
                for (int j=0;j<3;j++)
                new_normele[3*M1->npoin +6*nelem+ 3*i+j] = normedge[9*i+3+j];
            }
            // normal vector of edge1 p3 -> p1
            for (int i = 0; i < nelem; i++)
            {
                for (int j=0;j<3;j++)
                new_normele[3*M1->npoin +9*nelem+ 3*i+j] = normedge[9*i+6+j];
            }
           FunctionWithArgs prtelefield[] =
            {
                {"open", 1, nelem, open, SCA_int_VTK},
                {"BC", 1, nelem, cell_stat, SCA_int_VTK},
            };
        size_t countele = sizeof(prtelefield) / sizeof(prtelefield[0]);
        FunctionWithArgs prtpntfield[] = {
            {"normal_vec", 3, (M1->npoin + M1->numExtraPoints), new_normele, VEC_double_VTK}
        };
        size_t countpnt = 1;
        CHECK_ERROR(SaveVTK("./", "checkmesh", 0, M1, tri3funcVTK, prtelefield, countele, prtpntfield, countpnt));
        free (new_normele);
        free (extra_ptxyz);
    #endif
    // define sparse matrix of coefficient
    int *row_ptr = (int *)calloc((size_t)(nelem + 1), sizeof(int));
    int max_nnz = nelem * (Nredge+1); // At most non-zero entries per row (element+neighbours)
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
    int order_in_cells [] = {0,1,1,2,2,0};
    for (int ele = 0; ele < nelem; ele++)
    {
        nnz++;
        int IDele = nnz - 1;
        // Known cells 
        if (cell_stat[ele]>0){
            // equation for this element is u=Tb
            val[IDele] = 1;
            col_ind[IDele] = ele;
            RHS[ele] = (cell_stat[ele] == 1) ? 1000 : 0; 
        }else{
        // Unknown cells 
            for (int nei = 0; nei < Nredge; nei++)
            {
                int neighbor = esure[Nredge * ele + nei]; 
                int lp1 = elems [Nredge * ele + order_in_cells[2*nei]]-1; 
                int lp2 = elems [Nredge * ele + order_in_cells[2*nei+1]]-1; 
                double dA = 0;
                dA = SQUARE((ptxyz [3 * lp1 ] - ptxyz [3 * lp2 ] ));          // x-coordinate  
                dA += SQUARE((ptxyz [3 * lp1 + 1] - ptxyz [3 * lp2 +1 ] ));   // y-coordinate
                dA += SQUARE((ptxyz [3 * lp1 + 2] - ptxyz [3 * lp2 +2 ] ));   // z-coordinate
                dA = sqrt (dA); 
                double dl = 0;
                dl = SQUARE(cen [Nredge * ele] - cen [Nredge *neighbor]);     // x-coordinate
                dl += SQUARE(cen [Nredge * ele + 1] - cen [Nredge *neighbor + 1]);     // y-coordinate
                dl += SQUARE(cen [Nredge * ele + 2] - cen [Nredge *neighbor + 2]);     // z-coordinate
                dl = sqrt (dl);
                double coef = dA/dl;
                // Internal flux  
                // contribution for diagonal element of coefficient matrix 
                val[IDele] += coef;
                col_ind[IDele] = ele;
                // contribution for off-diagonal element of coefficient matrix 
                if (cell_stat[neighbor] >0 ){
                    // between known and unknown cells
                    RHS [ele] += (cell_stat[neighbor] == 1) ? 1000 : 0;
                }else{
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
    clock_t start_time,end_time;
    double cpu_time_used;
    SolverConfig config = {100000, 1e-8, true};
    solver_set_config(config);
    // precond_conjugate_gradient(&A, RHS, u);
    start_time = clock();
    conjugate_gradient(&A, RHS, u);
    end_time = clock();
    cpu_time_used = (double)(end_time - start_time) / CLOCKS_PER_SEC;
    printf("* CG Solver execution time : %.2f seconds\n", cpu_time_used);

    // calculate the gradient scaler on each element
    double *norm_grad;
    CHECK_ERROR(gradient_ele_tri3(nelem,elems,ptxyz,esure,normedge,u,&norm_grad));
    if (norm_grad == NULL ){
        fprintf(stderr,"! ERROR: grad array is empty\n");
        return -1;
    }
    for (int ele=0;ele<nelem;ele++){
        double sum = 0;
        sum += SQUARE (norm_grad[3*ele]);
        sum += SQUARE (norm_grad[3*ele+1]);
        sum += SQUARE (norm_grad[3*ele+2]);
        sum = sqrt (sum);
        if (sum != 0){
        norm_grad [3*ele] /=  sum;
        norm_grad [3*ele+1] /=  sum;
        norm_grad [3*ele+2] /=  sum;
        }
    }

    // calculate the third vector for local coordinate system for each cells 
    double *local_coord = (double *)malloc (9*(size_t)nelem*sizeof(double));
    for (int ele=0;ele<nelem;ele++){
        double thrid_vec [3];
        double first_vec [3] = {norm_grad[3*ele],norm_grad[3*ele+1],norm_grad[3*ele+2]};
        double second_vec [3] = {normele[3*ele],normele[3*ele+1],normele[3*ele+2]};
        crossProduct(first_vec,second_vec,thrid_vec);
        local_coord [9*ele]=first_vec[0];
        local_coord [9*ele+1]=first_vec[1];
        local_coord [9*ele+2]=first_vec[2];
        local_coord [9*ele+3]=second_vec[0];
        local_coord [9*ele+4]=second_vec[1];
        local_coord [9*ele+5]=second_vec[2];
        local_coord [9*ele+6]=thrid_vec[0];
        local_coord [9*ele+7]=thrid_vec[1];
        local_coord [9*ele+8]=thrid_vec[2];
    }
    // read stress tensor from log file 
    double *st;
    readfebiolog("temp/pres_0.log",nelem,&st,end_first_step);
    if (st==NULL){
        fprintf(stderr,"there is problem in reading stress tensor\n");
        exit(EXIT_FAILURE);
    }
    // analysis tensor in tangential plane 
    double *shear_st = (double *)malloc((size_t)nelem*4*sizeof(double));
    double *shear_evals_max = (double *)malloc((size_t)nelem*sizeof(double));
    double *shear_evals_min = (double *)malloc((size_t)nelem*sizeof(double));
    double *shear_evects_3D_min = (double *)malloc((size_t)nelem*3*sizeof(double));
    double *shear_evects_3D_max = (double *)malloc((size_t)nelem*3*sizeof(double));
    for (int ele=0;ele<nelem;ele++){
        // extract stress and normal, t1 and t2 tangantial vectors from arraies
        double stress [3][3],normal [3],t1 [3],t2 [3];
        for (int i =0;i<3;i++){
            for (int j =0;j<3;j++)
            stress [i][j]=st[9*ele+3*i+j];
        }
        for (int i =0;i<3;i++)
            normal [i] = local_coord[9*ele+i];
        for (int i =0;i<3;i++)
            t1[i]= local_coord[9*ele+3+i];
        for (int i =0;i<3;i++)
            t2[i]= local_coord[9*ele+6+i];
        
        // Compute the projection matrix
        double P[3][3];
        compute_projection_matrix(normal, P);

        // Compute the in-plane stress tensor
        double in_plane[3][3];
        compute_in_plane_stress(stress, P, in_plane);

        // Extract the 2x2 tangential stress tensor
        double tensor_2x2[2][2];
        extract_2x2_tensor(in_plane, t1, t2, tensor_2x2);

        //save in-plane stress tensor
        for (int i=0;i<2;i++) {
            for (int j=0;j<2;j++)
            shear_st[4*ele+i]=tensor_2x2[i][j];
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
        if (fabs(lambda1)>fabs(lambda2)){
            shear_evals_max [ele] = lambda1;
            shear_evals_min [ele] = lambda2;
            for (int i=0;i<3;i++)
                shear_evects_3D_max[3*ele+i]=v3D1[i];
            for (int i=0;i<3;i++)
                shear_evects_3D_min[3*ele+i]=v3D2[i];
        }else{
            shear_evals_max [ele] = lambda2;
            shear_evals_min [ele] = lambda1;
            for (int i=0;i<3;i++)
                shear_evects_3D_max[3*ele+i]=v3D2[i];
            for (int i=0;i<3;i++)
                shear_evects_3D_min[3*ele+i]=v3D1[i];
        }    
    }
    // disturbution of eval_min/eval_max 
    double *eval_ratio = (double *)malloc((size_t)nelem*sizeof(double));
    for (int ele=0;ele<nelem;ele++){
        eval_ratio[ele] = fabs(shear_evals_min[ele]/shear_evals_max[ele]);
    }
    // creat Histogram to show disturbution of eval_ratio
        int num_values = nelem;  // Number of values to analyze
        double max_value = 1;  // Maximum range value
        int num_bins = 20 ;     // Number of bins
        int *bins;

        // calculate the disturbution in each bin
        int max_bin_count=compute_disturbution_histogram(eval_ratio,num_values,&bins,num_bins,max_value);
        // save histogeram
        plot_histogram(bins,num_bins,max_value,num_values,"Eval_ratio.dat","Histogram of Eval min/max","ratio","frequency" ,"Eval_ratio_histogram.png","Evalratio",max_bin_count);
        free(bins);
    // find uni or bi-directional region
    int *sdir;
    unibimask(nelem,shear_evals_max,shear_evals_min,10000,&sdir);
    if (sdir == NULL){
        fprintf (stderr,"! there is problem in allocate memory for sdir.\n");
        exit(EXIT_FAILURE);
    }
    // save in VTK format
    // coordinate 
        M1->numExtraPoints = 3*nelem;
        double *extra_ptxyz2 = calloc((size_t)M1->numExtraPoints * 3, sizeof(double));
        // points for first vec 
        for (int i = 0; i < (nelem* 3); i++)
            extra_ptxyz2[i] = cen[i];
        // points for second vec 
        for (int i = 0; i < (nelem* 3); i++)
            extra_ptxyz2[3*nelem+i] = cen[i];
        // points for third vec    
        for (int i = 0; i < (nelem* 3); i++)
            extra_ptxyz2[6*nelem+i] = cen[i];

        M1->extra_ptxyz = extra_ptxyz2;

        double *new_normele2 = calloc(((size_t)M1->npoin + (size_t)M1->numExtraPoints) * 3, sizeof(double));
        // First vector of each element
        for (int i = 0; i < nelem; i++)
        {
            for (int j=0;j<3;j++)
            new_normele2[3*M1->npoin + 3*i+j] = local_coord[9*i+j];
        }
        // Second vector 
        for (int i = 0; i < nelem; i++)
        {
            for (int j=0;j<3;j++)
            new_normele2[3*M1->npoin +3*nelem+ 3*i+j] = local_coord[9*i+3+j];
        }
        // Third vector 
        for (int i = 0; i < nelem; i++)
        {
            for (int j=0;j<3;j++)
            new_normele2[3*M1->npoin +6*nelem+ 3*i+j] = local_coord[9*i+6+j];
        }

    FunctionWithArgs prtelefield2[] =
        {
            {"BC_poisson", 1, nelem, cell_stat, SCA_int_VTK},
            {"poisson", 1, nelem, u, SCA_double_VTK},
            {"uni/bi_region", 1, nelem, sdir, SCA_int_VTK},
            {"EValue_max", 1, nelem, shear_evals_max, SCA_double_VTK},
            {"EValue_min", 1, nelem, shear_evals_min, SCA_double_VTK},
            {"norm_grad_u", 3, nelem, norm_grad, VEC_double_VTK},
            {"Evect_max", 3, nelem, shear_evects_3D_max, VEC_double_VTK},
            {"Evect_min", 3, nelem, shear_evects_3D_min, VEC_double_VTK},
        };
    size_t countele2 = sizeof(prtelefield2) / sizeof(prtelefield2[0]);
    FunctionWithArgs prtpntfield2[] = {
        {"local_normal_vec", 3, (M1->npoin + M1->numExtraPoints), new_normele2, VEC_double_VTK}
    };
    size_t countpnt2 = 1;
    CHECK_ERROR(SaveVTK("./", "checkmesh", 1, M1, tri3funcVTK, prtelefield2, countele2, prtpntfield2, countpnt2));
    free (extra_ptxyz2);free(new_normele2);


    // free dynamics arraies
    free (elems);
    free (ptxyz);
    free (open);
    free (esure);
    free (esurp);
    free (esurp_ptr);
    free (row_ptr);
    free (val);
    free (col_ind);
    free (RHS);
    free (cen);
    free (M1);
    free (cell_stat);
    free (u);
    free (cenedge);
    free (normele);
    free (normedge);
    free (norm_grad);
    free (local_coord);
    free (st);
    free (shear_evals_max);
    free (shear_evals_min);
    free (shear_evects_3D_max);
    free (shear_evects_3D_min);
    free (sdir);
    free (eval_ratio);
    free (shear_st);
    return 0; // success signal
}