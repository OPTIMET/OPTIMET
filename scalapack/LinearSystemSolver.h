#ifndef OPTIMET_SCALAPACK_LINEAR_SYSTEM_SOLVER_H_
#define OPTIMET_SCALAPACK_LINEAR_SYSTEM_SOLVER_H_

#include "Types.h"
#ifdef OPTIMET_MPI

#include "mpi/Communicator.h"
#include "scalapack/Belos.h"
#include "scalapack/Matrix.h"

#ifdef OPTIMET_BELOS
#include <Teuchos_ParameterList.hpp>
#endif

#include <tuple>

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

#ifdef OPTIMET_BELOS
//! Solve a system of linear equations using Belos
template <class SCALAR>
std::tuple<Matrix<SCALAR>, int>
gmres_linear_system(Matrix<SCALAR> const &A, Matrix<SCALAR> const &b,
                    Teuchos::RCP<Teuchos::ParameterList> const &parameters,
                    mpi::Communicator const &comm = mpi::Communicator());
template <class SCALAR>
std::tuple<Matrix<SCALAR>, int>
gmres_linear_system(Matrix<SCALAR> const &A, Matrix<SCALAR> const &b,
                    mpi::Communicator const &comm = mpi::Communicator());
#endif
}
}

#include "scalapack/LinearSystemSolver.hpp"
#endif
#endif
