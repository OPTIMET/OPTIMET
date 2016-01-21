#ifndef OPTIMET_BLACS_H
#define OPTIMET_BLACS_H

#include "Types.h"
#include "OptimetFC.h"

extern "C" {
  void OPTIMET_FC_GLOBAL(sl_init, SL_INIT)(int *ictxt, int *nprow, int *npcol);
  void OPTIMET_FC_GLOBAL(blacs_gridinfo, BLACS_GRIDINFO)(
      int *context, int *nrows, int *ncols, int *this_row, int *this_col);
  void OPTIMET_FC_GLOBAL(blacs_gridexit, BLACS_GRIDEXIT)(int *ictxt);
  void OPTIMET_FC_GLOBAL(blacs_exit, BLACS_EXIT)(int *status);
} /* optimet */
#endif
