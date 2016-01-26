#include <vector>
#include "scalapack/LinearSystemSolver.h"
#include "scalapack/Matrix.h"
#include "scalapack/InitExit.h"
#include "scalapack/Blacs.h"

namespace optimet {
namespace scalapack {

template <class SCALAR>
std::tuple<Matrix<SCALAR>, int>
general_linear_system(Matrix<SCALAR> const &A, Matrix<SCALAR> const &b) {
  if(not(A.context().is_valid() and b.context().is_valid()))
    return {b, 0};
  Matrix<SCALAR> result = b;
  Matrix<SCALAR> Acopy = A;
  auto info = general_linear_system_inplace(Acopy, result);
  return {std::move(result), std::move(info)};
}

namespace {
#define OPTIMET_MACRO(letter, LETTER, TYPE)                                                        \
  inline void gesv(int *n, int *nrhs, TYPE *a, int *ia, int *ja, int *desca, int *ipiv, TYPE *b,   \
                   int *ib, int *jb, int *descb, int *info) {                                      \
    OPTIMET_FC_GLOBAL(p##letter##gesv, P##LETTER##GESV)(                                           \
        n, nrhs, a, ia, ja, desca, ipiv, b, ib, jb, descb, info); \
  }
OPTIMET_MACRO(s, S, float);
OPTIMET_MACRO(d, D, double);
OPTIMET_MACRO(c, C, std::complex<float>);
OPTIMET_MACRO(z, Z, std::complex<double>);
#undef OPTIMET_MACRO
}

template <class SCALAR> int general_linear_system_inplace(Matrix<SCALAR> &A, Matrix<SCALAR> &b) {
  if(A.rows() != A.cols())
    throw std::runtime_error("Matrix should be square");
  if(A.cols() != b.rows())
    throw std::runtime_error("A.cols() != b.rows()");
  if(not(A.context().is_valid() and b.context().is_valid()))
    return 0;
  if(A.context() != b.context())
    throw std::runtime_error("Contexts of A and b must be identical");
  if(A.blocks().rows != A.blocks().cols)
    throw std::runtime_error("Blocs must be square");
  if(b.blocks().rows != b.blocks().cols)
    throw std::runtime_error("Blocs must be square");

  int n = A.rows(), nrhs = b.cols(), one = 1, info;
  std::vector<int> ipiv(A.local().rows() + A.blocks().rows);
  gesv(&n, &nrhs, A.local().data(), &one, &one, const_cast<int *>(A.blacs().data()), ipiv.data(),
       b.local().data(), &one, &one, const_cast<int *>(b.blacs().data()), &info);
  return info;
}

} // scalapack
} // optimet
