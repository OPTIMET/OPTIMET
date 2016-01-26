#include <iostream>
#include "scalapack/Matrix.h"
#include "scalapack/InitExit.h"
#include "scalapack/Blacs.h"
#include "AttachDebug.h"

namespace optimet {
namespace scalapack {
// Matrix Matrix::transfer(Context other) const {
// }

#define OPTIMET_MACRO(NAME)                                                                        \
  Matrix::EigenMatrix::Index Matrix::NAME##s(Context const &context, Sizes size, Sizes blocks,     \
                                             Index index) {                                        \
    if(not context.is_valid())                                                                     \
      return 0;                                                                                    \
    int n = static_cast<int>(size.NAME##s);                                                        \
    int nb = static_cast<int>(blocks.NAME##s);                                                     \
    int iproc = static_cast<int>(context.NAME());                                                  \
    int isrc = static_cast<int>(index.NAME);                                                       \
    int nprocs = static_cast<int>(context.NAME##s());                                              \
    auto const result = OPTIMET_FC_GLOBAL(numroc, NUMROC)(&n, &nb, &iproc, &isrc, &nprocs);        \
    return static_cast<Matrix::EigenMatrix::Index>(result);                                        \
  }
OPTIMET_MACRO(row);
OPTIMET_MACRO(col);
#undef OPTIMET_MACRO

void Matrix::transfer_to(Context const &un, Matrix &other) const {
  int m = rows();
  int n = cols();
  int one = 1;
  int c = *un;
  // wait_for_debugger();
  OPTIMET_FC_GLOBAL(pdgemr2d, PDGEMR2D)(&n, &m, const_cast<t_real *>(eigen().data()), &one, &one,
            const_cast<int*>(blacs().data()), other.eigen().data(), &one, &one,
            const_cast<int*>(other.blacs().data()), &c);
}
Matrix Matrix::transfer_to(Context const &un, Context const &other, Sizes const &_blocks,
                           Index const &_index) const {
  Matrix result(other, sizes(),
                {_blocks.rows == std::numeric_limits<t_uint>::max() ? blocks().rows : _blocks.rows,
                 _blocks.cols == std::numeric_limits<t_uint>::max() ? blocks().cols : _blocks.cols},
                {_index.row == std::numeric_limits<t_uint>::max() ? index().row : _index.row,
                 _index.col == std::numeric_limits<t_uint>::max() ? index().col : _index.col});
  transfer_to(un, result);
  return result;
}
} // scalapack
} // optimet
