#include "scalapack/Blacs.h"
#include "scalapack/InitExit.h"
#include "scalapack/Matrix.h"
#include <iostream>

namespace optimet {
namespace scalapack {

#define OPTIMET_MACRO(NAME)                                                                        \
  template <class SCALAR>                                                                          \
  typename Matrix<SCALAR>::EigenMatrix::Index Matrix<SCALAR>::NAME##s(                             \
      Context const &context, Sizes size, Sizes blocks, Index index) {                             \
    if(not context.is_valid())                                                                     \
      return 0;                                                                                    \
    int n = static_cast<int>(size.NAME##s);                                                        \
    int nb = static_cast<int>(blocks.NAME##s);                                                     \
    int iproc = static_cast<int>(context.NAME());                                                  \
    int isrc = static_cast<int>(index.NAME);                                                       \
    int nprocs = static_cast<int>(context.NAME##s());                                              \
    auto const result = OPTIMET_FC_GLOBAL(numroc, NUMROC)(&n, &nb, &iproc, &isrc, &nprocs);        \
    return static_cast<typename Matrix<SCALAR>::EigenMatrix::Index>(result);                       \
  }
OPTIMET_MACRO(row);
OPTIMET_MACRO(col);
#undef OPTIMET_MACRO

namespace {
#define OPTIMET_MACRO(LETTER, TYPE)                                                                \
  inline void gemr2d(int m, int n, TYPE *ptrmyblock, int ia, int ja, int *ma, TYPE *ptrmynewblock, \
                     int ib, int jb, int *mb, int globcontext) {                                   \
    Cp##LETTER##gemr2d(m, n, ptrmyblock, ia, ja, ma, ptrmynewblock, ib, jb, mb, globcontext);      \
  }
OPTIMET_MACRO(s, float);
OPTIMET_MACRO(d, double);
OPTIMET_MACRO(c, std::complex<float>);
OPTIMET_MACRO(z, std::complex<double>);
#undef OPTIMET_MACRO
}

template <class SCALAR> void Matrix<SCALAR>::transfer_to(Context const &un, Matrix &other) const {
  if(context().is_valid() or un.is_valid() or other.context().is_valid())
    gemr2d(rows(), cols(), const_cast<SCALAR *>(local().data()), 1, 1,
           const_cast<int *>(blacs().data()), other.local().data(), 1, 1,
           const_cast<int *>(other.blacs().data()), *un);
}
template <class SCALAR>
Matrix<SCALAR> Matrix<SCALAR>::transfer_to(Context const &un, Context const &other,
                                           Sizes const &_blocks, Index const &_index) const {
  Matrix<SCALAR> result(
      other, sizes(),
      {_blocks.rows == std::numeric_limits<t_uint>::max() ? blocks().rows : _blocks.rows,
       _blocks.cols == std::numeric_limits<t_uint>::max() ? blocks().cols : _blocks.cols},
      {_index.row == std::numeric_limits<t_uint>::max() ? index().row : _index.row,
       _index.col == std::numeric_limits<t_uint>::max() ? index().col : _index.col});
  transfer_to(un, result);
  return result;
}
} // scalapack
} // optimet
