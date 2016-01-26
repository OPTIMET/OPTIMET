#ifndef OPTIMET_SCALAPACK_LINEAR_SYSTEM_SOLVER_H_
#define OPTIMET_SCALAPACK_LINEAR_SYSTEM_SOLVER_H_

#include <tuple>
#include "Types.h"
#include "scalapack/Matrix.h"

namespace optimet {
namespace scalapack {

  //! Solves a system of linear equations, overwriting the input
  //! \param A: matrix to diagnolize on input, factors L and U on output
  //! \param b: right-hand size on input, X on output
  int general_linear_system_inplace(Matrix &A, Matrix &b);
  //! Solves a system of linear equations
  std::tuple<Matrix, int> general_linear_system(Matrix const &A, Matrix const &b);
}
}
#endif
