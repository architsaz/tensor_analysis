#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "febio_types.h"
#include "common.h"

// read stress tensor in the log file of febio
int readfebiolog(char *path, int mesh_nelem, double **st2, read_time logtime)
{

    int e = 0;
    double *st;
    // Open log file :
    FILE *fptr = fopen(path, "r");
    char str[256];
    if (fptr == NULL)
    {
        fprintf(stderr, "* ERROR: there is no file in this path : %s.\n", path);
        exit(EXIT_FAILURE);
    }
    // find the number of shell element
    int nelem = 0;
    while (fgets(str, 256, fptr) != NULL)
    {
        if (sscanf(str, "	Number of shell elements ....................... : %d", &nelem) == 1)
            break;
    }
    fclose(fptr);
    printf("* number of shell element is %d\n", nelem);
    if (mesh_nelem != nelem)
    {
        fprintf(stderr, "ERROR: Number of shell element in the file : %s does not match with what zfem file on data directory", path);
        exit(EXIT_FAILURE);
    }
    // find the maximum time in the file
    double time_value, max_time;
    max_time = 0.0;
    fptr = fopen(path, "r");
    while (fgets(str, 256, fptr) != NULL)
    {
        if (sscanf(str, "Time = %lf", &time_value) == 1)
            max_time = MAX(time_value, max_time);
    }
    fclose(fptr);
    printf("* %s study execute till time = %lf\n", path, max_time);
    // save the stress tensor
    st = calloc(9 * (size_t)nelem, sizeof(*st));
    double logtime_value = 0.0;
    if (logtime == end_first_step)
        logtime_value = 1.0;
    if (logtime == end_second_step)
        logtime_value = 2.0;
    if (logtime == time_max)
        logtime_value = max_time;
    // st [9] = [sxx,sxy,sxz;syx,syy,syz;szx,szy,szz]
    fptr = fopen(path, "r");
    int junk, nscan = 0;
    int find_time = 0;
    while (fgets(str, 256, fptr) != NULL)
    {
        if (sscanf(str, "Time = %lf", &time_value) == 1)
        {
            if (time_value == logtime_value)
            {
                printf("* start to read stress tensor at time: %lf\n", time_value);
                find_time++;
                fgets(str, 256, fptr);
                for (int ele = 0; ele < nelem; ele++)
                {
                    fgets(str, 256, fptr);
                    nscan = 0;
                    // log_st[6]= [sxx(0),syy(4),szz(8),sxy(1)(3),syz(5)(7),sxz(2)(6)]
                    nscan = sscanf(str, "%d %lf %lf %lf %lf %lf %lf", &junk, &st[9 * ele], &st[9 * ele + 4], &st[9 * ele + 8],
                                   &st[9 * ele + 1], &st[9 * ele + 5], &st[9 * ele + 2]);
                    if (nscan != 7)
                    {
                        fprintf(stderr, "there is error on number of element in line %d", ele);
                        exit(EXIT_FAILURE);
                    }
                    st[9 * ele + 3] = st[9 * ele + 1];
                    st[9 * ele + 7] = st[9 * ele + 5];
                    st[9 * ele + 6] = st[9 * ele + 2];
                }
            }
        }
    }
    if (find_time == 0)
    {
        fprintf(stderr, "ERROR: can not find time %lf in the log file!\n", max_time);
        exit(EXIT_FAILURE);
    }
    // for(int ele=0;ele<10;ele++){
    //     printf("ele: %d ",ele);
    //     for(int i =0;i<9;i++) printf("%lf ",st[9*ele+i]);
    //     printf("\n");
    // }
    fclose(fptr);
    printf("* the stress tensor saved !\n");

    *st2 = st;
    return e;
}
// find unidirectional or bidirectional stress region mask:
void unibimask(int nelem, double *smax, double *smin, double critical_ratio, double threshold, int **sdir2)
{
    int *sdir = (int *)calloc((size_t)nelem, sizeof(int));
    for (int ele = 0; ele < nelem; ele++)
    {
        if (fabs(smin[ele]) < 1)
            continue; // avoid unpressurized region
        double ratio = fabs(smin[ele]) / fabs(smax[ele]);
        sdir[ele] = (ratio < critical_ratio) ? 1 : 2; // unidiretional : 1 bidirectional : 2
        if (fabs(smax[ele]) <= threshold)
        {
            // sdir[ele] = (sdir[ele] == 4) ? 2 : 1;
            sdir[ele] = 0;
        }
    }
    *sdir2 = sdir;
}