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

#ifndef OPTIMET_SCALAPACK_LINEAR_SYSTEM_SOLVER_HPP_
#define OPTIMET_SCALAPACK_LINEAR_SYSTEM_SOLVER_HPP_

#include "scalapack/Blacs.h"
#include "scalapack/InitExit.h"
#include "scalapack/LinearSystemSolver.h"
#include "scalapack/Matrix.h"

#include <vector>

#ifdef OPTIMET_BELOS
#include <BelosBlockGmresSolMgr.hpp>
#include <BelosConfigDefs.hpp>
#include <BelosLinearProblem.hpp>
#include <BelosPseudoBlockGmresSolMgr.hpp>
#include <BelosSolverFactory.hpp>
#include <BelosTpetraAdapter.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Tpetra_MultiVector.hpp>
#endif

namespace optimet {
namespace scalapack {

template <class SCALAR>
std::tuple<typename Matrix<SCALAR>::ConcreteMatrix, int>
general_linear_system(Matrix<SCALAR> const &A, Matrix<SCALAR> const &b) {
  typedef typename Matrix<SCALAR>::ConcreteMatrix ConcreteMatrix;
  if(not(A.context().is_valid() and b.context().is_valid()))
    return std::tuple<ConcreteMatrix, int>{ConcreteMatrix(b.context(), b.sizes(), b.blocks()), 0};
  ConcreteMatrix result(b.local(), b.context(), b.sizes(), b.blocks());
  ConcreteMatrix Acopy(A.local(), A.context(), A.sizes(), A.blocks());
  auto info = general_linear_system_inplace(Acopy, result);
  return std::tuple<ConcreteMatrix, int>{std::move(result), std::move(info)};
}

namespace {
#define OPTIMET_MACRO(letter, LETTER, TYPE)                                                        \
  inline void gesv(int *n, int *nrhs, TYPE *a, int *ia, int *ja, int *desca, int *ipiv, TYPE *b,   \
                   int *ib, int *jb, int *descb, int *info) {                                      \
    OPTIMET_FC_GLOBAL(p##letter##gesv, P##LETTER##GESV)                                            \
    (n, nrhs, a, ia, ja, desca, ipiv, b, ib, jb, descb, info);                                     \
  }
OPTIMET_MACRO(s, S, float);
OPTIMET_MACRO(d, D, double);
OPTIMET_MACRO(c, C, std::complex<float>);
OPTIMET_MACRO(z, Z, std::complex<double>);
#undef OPTIMET_MACRO

template <class SCALARA, class SCALARB>
void sane_input(Matrix<SCALARA> const &A, Matrix<SCALARB> const &b) {
  if(A.rows() != A.cols())
    throw std::runtime_error("Matrix should be square");
  if(A.cols() != b.rows())
    throw std::runtime_error("A.cols() != b.rows()");
  if(A.context().is_valid() and b.context().is_valid()) {
    if(A.context() != b.context())
      throw std::runtime_error("Contexts of A and b must be identical");
    if(A.blocks().rows != A.blocks().cols)
      throw std::runtime_error("Blocs must be square");
    if(b.blocks().rows != b.blocks().cols)
      throw std::runtime_error("Blocs must be square");
  }
}
}

template <class SCALAR> int general_linear_system_inplace(Matrix<SCALAR> &A, Matrix<SCALAR> &b) {
  sane_input(A, b);
  int n = A.rows(), nrhs = b.cols(), one = 1, info;
  std::vector<int> ipiv(A.local().rows() + A.blocks().rows);
  gesv(&n, &nrhs, A.local().data(), &one, &one, const_cast<int *>(A.blacs().data()), ipiv.data(),
       b.local().data(), &one, &one, const_cast<int *>(b.blacs().data()), &info);
  return info;
}

#ifdef OPTIMET_BELOS
template <class SCALARA, class SCALARB>
std::tuple<typename Matrix<SCALARA>::ConcreteMatrix, int>
gmres_linear_system(Matrix<SCALARA> const &A, Matrix<SCALARB> const &b,
                    Teuchos::RCP<Teuchos::ParameterList> const &parameters,
                    mpi::Communicator const &comm) {
  typedef typename Matrix<SCALARA>::ConcreteMatrix ConcreteMatrix;
  sane_input(A, b);
  if(b.context() != A.context()) {
    auto const b_in_A = b.transfer_to(A.context());
    return gmres_linear_system(A, b_in_A, parameters, comm);
  }
  if(A.size() == 0)
    return std::tuple<ConcreteMatrix, int>{ConcreteMatrix(b.context(), b.sizes(), b.blocks()), 0};
  auto const splitcomm = comm.split(A.context().is_valid());
  if(not A.context().is_valid())
    return std::tuple<ConcreteMatrix, int>{ConcreteMatrix(b.context(), b.sizes(), b.blocks()), 0};
  // Create the GMRES solver.
  BelosSolverFactory<SCALARA> factory;
  auto solver = factory.create(parameters->get("Solver", "GMRES"), parameters);
  // Create a LinearProblem struct with the problem to solve.
  // A, X, B, and M are passed by (smart) pointer, not copied.
  Teuchos::RCP<BelosOperator<SCALARA>> Aptr = Teuchos::rcp(new BelosOperator<SCALARA>(A));
  auto const X = tpetra_vector(b, splitcomm);
  auto const B = tpetra_vector(b, splitcomm);
  Teuchos::RCP<BelosLinearProblem<SCALARA>> problem =
      rcp(new BelosLinearProblem<SCALARA>(Aptr, X, B));
  // problem->setRightPrec(M);
  // Tell the solver what problem you want to solve.
  if(not problem->setProblem())
    throw std::runtime_error("Could not setup up Belos problem");
  solver->setProblem(problem);

  // Print out parameters for given verbosity
  if(parameters->get<int>("Verbosity", 0) & Belos::MsgType::FinalSummary) {
    auto const out = parameters->get<Teuchos::RCP<std::ostream>>("Output Stream",
                                                                 Teuchos::rcp(&std::cout, false));
    if(splitcomm.rank() == 0)
      solver->getCurrentParameters()->print(*out);
  }
  // Attempt to solve the linear system.  result == Belos::Converged
  // means that it was solved to the desired tolerance.  This call
  // overwrites X with the computed approximate solution.

  int info = (solver->solve() == Belos::Converged ? 0 : 1);
  auto const view = as_matrix(*X, b);
  ConcreteMatrix result(b.context(), b.sizes(), b.blocks());
  result.local() = view.local();
  return std::tuple<ConcreteMatrix, int>{std::move(result), std::move(info)};
} // anonymous namespace

template <class SCALARA, class SCALARB>
std::tuple<typename Matrix<SCALARA>::ConcreteMatrix, int>
gmres_linear_system(Matrix<SCALARA> const &A, Matrix<SCALARB> const &b,
                    mpi::Communicator const &comm) {
  // Make an empty new parameter list.
  Teuchos::RCP<Teuchos::ParameterList> solverParams = Teuchos::parameterList();
  // Set some GMRES parameters.
  //
  // "Num Blocks" = Maximum number of Krylov vectors to store.  This
  // is also the restart length.  "Block" here refers to the ability
  // of this particular solver (and many other Belos solvers) to solve
  // multiple linear systems at a time, even though we are only solving
  // one linear system in this example.
  solverParams->set("Num Blocks", 1000);
  solverParams->set("Maximum Iterations", 4000);
  solverParams->set("Convergence Tolerance", 1.0e-10);
  return gmres_linear_system(A, b, solverParams, comm);
}
#endif

} // scalapack
} // optimet
#endif
