#ifndef OPTIMET_SCALAPACK_LINEAR_SYSTEM_SOLVER_H_
#define OPTIMET_SCALAPACK_LINEAR_SYSTEM_SOLVER_H_

#include "Types.h"
#ifdef OPTIMET_MPI

#include <tuple>
#include "scalapack/Matrix.h"

namespace optimet {
namespace scalapack {

//! Solves a system of linear equations, overwriting the input
//! \param A: matrix to diagnolize on input, factors L and U on output
//! \param b: right-hand size on input, X on output
template <class SCALAR> int general_linear_system_inplace(Matrix<SCALAR> &A, Matrix<SCALAR> &b);
//! Solves a system of linear equations
template <class SCALAR>
std::tuple<Matrix<SCALAR>, int>
general_linear_system(Matrix<SCALAR> const &A, Matrix<SCALAR> const &b);

//! Solve a system of linear equations using Belos
template <class SCALAR>
std::tuple<Matrix<SCALAR>, int>
gmres_linear_system(Matrix<SCALAR> const &A, Matrix<SCALAR> const &b);
}
}

# include "scalapack/LinearSystemSolver.hpp"
#endif
#endif
