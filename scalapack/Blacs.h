#ifndef OPTIMET_BLACS_H
#define OPTIMET_BLACS_H

#include <complex>
#include "OptimetFC.h"

extern "C" {
void OPTIMET_FC_GLOBAL_(sl_init, SL_INIT)(int *ictxt, int *nprow, int *npcol);
void OPTIMET_FC_GLOBAL_(blacs_setup, BLACS_SETUP)(int *nprow, int *npcol);
void OPTIMET_FC_GLOBAL_(blacs_gridinfo, BLACS_GRIDINFO)(int *context, int *nrows, int *ncols,
                                                        int *this_row, int *this_col);
void OPTIMET_FC_GLOBAL_(blacs_gridexit, BLACS_GRIDEXIT)(int *ictxt);
void OPTIMET_FC_GLOBAL_(blacs_exit, BLACS_EXIT)(int *status);
void OPTIMET_FC_GLOBAL_(blacs_pinfo, BLACS_PINFO)(int *rank, int *size);
void OPTIMET_FC_GLOBAL_(blacs_get, BLACS_GET)(int *context, int *what, int *output);
void OPTIMET_FC_GLOBAL_(blacs_gridinit, BLACS_GRIDINIT)(int *context, char *order, int *rows,
                                                        int *cols);
int OPTIMET_FC_GLOBAL(numroc, NUMROC)(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
void Cpsgemr2d(int m, int n, float *ptrmyblock, int ia, int ja, int *ma, float *ptrmynewblock,
               int ib, int jb, int *mb, int globcontext);
void Cpdgemr2d(int m, int n, double *ptrmyblock, int ia, int ja, int *ma, double *ptrmynewblock,
               int ib, int jb, int *mb, int globcontext);
void Cpcgemr2d(int m, int n, std::complex<float> *ptrmyblock, int ia, int ja, int *ma,
               std::complex<float> *ptrmynewblock, int ib, int jb, int *mb, int globcontext);
void Cpzgemr2d(int m, int n, std::complex<double> *ptrmyblock, int ia, int ja, int *ma,
               std::complex<double> *ptrmynewblock, int ib, int jb, int *mb, int globcontext);
void OPTIMET_FC_GLOBAL(psgesv, PSGESV)(int *n, int *nrhs, float *a, int *ia, int *ja, int *desca,
                                       int *ipiv, float *b, int *ib, int *jb, int *descb,
                                       int *info);
void OPTIMET_FC_GLOBAL(pdgesv, PDGESV)(int *n, int *nrhs, double *a, int *ia, int *ja, int *desca,
                                       int *ipiv, double *b, int *ib, int *jb, int *descb,
                                       int *info);
void OPTIMET_FC_GLOBAL(pcgesv, PCGESV)(int *n, int *nrhs, std::complex<float> *a, int *ia, int *ja,
                                       int *desca, int *ipiv, std::complex<float> *b, int *ib,
                                       int *jb, int *descb, int *info);
void OPTIMET_FC_GLOBAL(pzgesv, PZGESV)(int *n, int *nrhs, std::complex<double> *a, int *ia, int *ja,
                                       int *desca, int *ipiv, std::complex<double> *b, int *ib,
                                       int *jb, int *descb, int *info);
} /* optimet */
#endif
