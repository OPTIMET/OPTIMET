#include <vector>
#include "scalapack/LinearSystemSolver.h"
#include "scalapack/Matrix.h"
#include "scalapack/InitExit.h"
#include "scalapack/Blacs.h"

namespace optimet {
namespace scalapack {

std::tuple<Matrix, int> general_linear_system(Matrix const &A, Matrix const &b) {
  if(not (A.context().is_valid() and b.context().is_valid()))
    return {b, 0};
  Matrix result = b;
  Matrix Acopy = A;
  auto info = general_linear_system_inplace(Acopy, result);
  return {std::move(result), std::move(info)};
}

int general_linear_system_inplace(Matrix &A, Matrix &b) {
  if(A.rows() != A.cols())
    throw std::runtime_error("Matrix should be square");
  if(A.cols() != b.rows())
    throw std::runtime_error("A.cols() != b.rows()");
  if(not (A.context().is_valid() and b.context().is_valid()))
    return 0;
  if(A.context() != b.context())
    throw std::runtime_error("Contexts of A and b must be identical");
  if(A.blocks().rows != A.blocks().cols)
    throw std::runtime_error("Blocs must be square");
  if(b.blocks().rows != b.blocks().cols)
    throw std::runtime_error("Blocs must be square");

  int n = A.rows(), nrhs = b.cols(), one = 1, info;
  std::vector<int> ipiv(A.local().rows() + A.blocks().rows);
  OPTIMET_FC_GLOBAL(pdgesv, PDGESV)(&n, &nrhs, A.local().data(), &one, &one,
      const_cast<int *>(A.blacs().data()), ipiv.data(), b.local().data(), &one, &one,
      const_cast<int *>(b.blacs().data()), &info);
  return info;
}

} // scalapack
} // optimet
