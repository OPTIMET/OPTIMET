#ifndef OPTIMET_BLACS_H
#define OPTIMET_BLACS_H

#include "Types.h"
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
void OPTIMET_FC_GLOBAL(pdgemr2d, PDGEMR2D)(int *m, int *n, double *A, int *IA, int *JA, int *descA,
                                           double *B, int *IB, int *JB, int *descB, int *gcontext);
int OPTIMET_FC_GLOBAL(numroc, NUMROC)(int *n, int *nb, int *iproc, int *isrcproc, int *nprocs);
void Cpdgemr2d(int m, int n, double *ptrmyblock, int ia, int ja, int *ma, double *ptrmynewblock,
               int ib, int jb, int *mb, int globcontext);

} /* optimet */
#endif
