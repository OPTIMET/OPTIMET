#include <iostream>
#include "scalapack/Matrix.h"
#include "scalapack/InitExit.h"
#include "scalapack/blacs.h"

namespace optimet {
namespace scalapack {
// Matrix Matrix::transfer(Context other) const {
// }

#define OPTIMET_MACRO(NAME)                                                                      \
  Matrix::EigenMatrix::Index                                                                     \
  Matrix::NAME ## s(Context const &context, Sizes size, Sizes blocks, Index index) {             \
    if(not context.is_valid())                                                                   \
      return 0;                                                                                  \
    int n = static_cast<int>(size.NAME ## s);                                                    \
    int nb = static_cast<int>(blocks.NAME ## s);                                                 \
    int iproc = static_cast<int>(context.NAME());                                                \
    int isrc = static_cast<int>(index.NAME) + 1;                                                 \
    int nprocs = static_cast<int>(context.NAME ## s());                                          \
    auto const result = OPTIMET_FC_GLOBAL(numroc, NUMROC)(&n, &nb, &iproc, &isrc, &nprocs);      \
    return static_cast<Matrix::EigenMatrix::Index>(result);                                      \
  }
  OPTIMET_MACRO(row);
  OPTIMET_MACRO(col);
#undef OPTIMET_MACRO
} // scalapack
} // optimet
