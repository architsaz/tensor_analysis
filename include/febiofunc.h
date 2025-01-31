#ifndef FEBIOFUNCS_H
#define FEBIOFUNCS_H
    #include "febio_types.h"
    int readfebiolog(char *path,int nelem, double **st2, read_time logtime);
    void unibimask(int nelem, double *smax, double *smin,double threshold, int **sdir2);
#endif