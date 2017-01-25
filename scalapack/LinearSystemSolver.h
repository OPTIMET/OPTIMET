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

#ifndef OPTIMET_SCALAPACK_LINEAR_SYSTEM_SOLVER_H_
#define OPTIMET_SCALAPACK_LINEAR_SYSTEM_SOLVER_H_

#include "Types.h"
#ifdef OPTIMET_SCALAPACK

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
std::tuple<typename Matrix<SCALAR>::ConcreteMatrix, int>
general_linear_system(Matrix<SCALAR> const &A, Matrix<SCALAR> const &b);

#ifdef OPTIMET_BELOS
//! Solve a system of linear equations using Belos
template <class SCALARA, class SCALARB>
std::tuple<typename Matrix<SCALARA>::ConcreteMatrix, int>
gmres_linear_system(Matrix<SCALARA> const &A, Matrix<SCALARB> const &b,
                    Teuchos::RCP<Teuchos::ParameterList> const &parameters,
                    mpi::Communicator const &comm = mpi::Communicator());
template <class SCALARA, class SCALARB>
std::tuple<typename Matrix<SCALARA>::ConcreteMatrix, int>
gmres_linear_system(Matrix<SCALARA> const &A, Matrix<SCALARB> const &b,
                    mpi::Communicator const &comm = mpi::Communicator());
#endif
}
}

#include "scalapack/LinearSystemSolver.hpp"
#endif
#endif
