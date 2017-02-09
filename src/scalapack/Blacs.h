// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

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
void OPTIMET_FC_GLOBAL(blacs_gridmap, BLACS_GRIDMAP)(
    int *context, int *usermap, int *ldau, int *nprow, int *npcol);
int OPTIMET_FC_GLOBAL_(blacs_pnum, BLACS_PNUM)(int *context, int *i, int *j);
void OPTIMET_FC_GLOBAL_(blacs_pcoord, BLACS_PCOORD)(int *context, int *pnum, int *row, int *col);

#define OPTIMET_MACRO(letter, LETTER, TYPE)                                                   \
  void OPTIMET_FC_GLOBAL(letter ## gebs2d, LETTER ## GEBS2D)(                                 \
      int *context, char const *, char const *, int *m, int *n, TYPE const*, int *lda);       \
  void OPTIMET_FC_GLOBAL(letter ## gebr2d, LETTER ## GEBR2D)(                                 \
      int *context, char const *, char const *, int *m, int *n, TYPE*,                        \
      int *lda, int *rsrc, int *csrc);                                                        \
  void Cp ## letter ## gemr2d(int m, int n, TYPE *ptrmyblock, int ia, int ja, int *ma, TYPE   \
      *ptrmynewblock, int ib, int jb, int *mb, int globcontext);                              \
  void OPTIMET_FC_GLOBAL(p ## letter ## gesv, P ## LETTER ## GESV)(int *n, int *nrhs,         \
      TYPE *a, int *ia, int *ja, int *desca, int *ipiv, TYPE *b, int *ib, int *jb,            \
      int *descb, int *info);

OPTIMET_MACRO(i, I, int);
OPTIMET_MACRO(s, S, float);
OPTIMET_MACRO(d, D, double);
OPTIMET_MACRO(c, C, std::complex<float>);
OPTIMET_MACRO(z, Z, std::complex<double>);
#undef OPTIMET_MACRO

#define OPTIMET_MACRO(func, FUNC)                                                             \
  int OPTIMET_FC_GLOBAL(indx ## func, INDX ## FUNC)(int*, int*, int*, int*, int*)
OPTIMET_MACRO(g2l, G2L);
OPTIMET_MACRO(l2g, L2G);
OPTIMET_MACRO(g2p, G2P);
#undef OPTIMET_MACRO

#define OPTIMET_MACRO(func, FUNC, TYPE)                                                           \
  int OPTIMET_FC_GLOBAL(p ## func ##gemm, p ## FUNC ## GEMM)(                                     \
      char* TRANSA, char* TRANSB, int * M, int * N, int * K,                                      \
      TYPE * ALPHA, TYPE * A, int * IA, int * JA, int * DESCA,                                    \
      TYPE * B, int * IB, int * JB, int * DESCB,                                                  \
      TYPE * BETA, TYPE * C, int * IC, int * JC, int * DESCC);
OPTIMET_MACRO(d, D, double);
OPTIMET_MACRO(z, Z, std::complex<double>);
#undef OPTIMET_MACRO
} /* optimet */
#endif
