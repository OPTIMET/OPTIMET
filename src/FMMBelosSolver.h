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

#ifndef OPTIMET_FMM_BELOS_SOLVER_H
#define OPTIMET_FMM_BELOS_SOLVER_H

#include "Run.h"
#include "Solver.h"
#include "Types.h"
#include "mpi/FastMatrixMultiply.h"
#include <limits>

namespace optimet {
namespace solver {

#ifdef OPTIMET_MPI
#ifdef OPTIMET_BELOS
//! Belos optimizer using the Fast Matrix Multiply
class FMMBelos : public AbstractSolver {
public:
  FMMBelos(
      std::shared_ptr<Geometry> geometry, std::shared_ptr<Excitation const> incWave,
      mpi::Communicator const &comm = mpi::Communicator(),
      Teuchos::RCP<Teuchos::ParameterList> belos_params = Teuchos::rcp(new Teuchos::ParameterList),
      t_int subdiagonals = std::numeric_limits<t_int>::max())
      : AbstractSolver(geometry, incWave, comm), fmm_(nullptr), belos_params_(belos_params),
        subdiagonals(subdiagonals) {
    update();
  }

  FMMBelos(Run const &run)
      : FMMBelos(run.geometry, run.excitation, run.communicator, run.belos_params,
                 run.fmm_subdiagonals) {}

  ~FMMBelos(){};

  /**
   * Solve the scattered and internal coefficients using the method specified by
   * solverMethod.
   * @param X_sca_ the return vector for the scattered coefficients.
   * @param X_int_ the return vector for the internal coefficients.
   * @return 0 if successful, 1 otherwise.
   */
  void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_,Vector<t_complex> &X_sca_SH,
             Vector<t_complex> &X_int_SH, std::vector<double *> CGcoeff) const override;
  //! \brief Update after internal parameters changed externally
  //! \details Because that's how the original implementation rocked.
  virtual void update() override;

  //! \brief Parameters for Belos/Trilinos solvers
  //! \note Mere access to the parameters requires the Teuchos::ParameterList to be modifiable. So
  //! the constness is not quite respected here.
  Teuchos::RCP<Teuchos::ParameterList> belos_parameters() const { return belos_params_; }

protected:
  //! Fast-matrix multiply operator
  std::shared_ptr<mpi::FastMatrixMultiply> fmm_;
  //! Parameter list of the belos solvers
  Teuchos::RCP<Teuchos::ParameterList> belos_params_;
  //! The local field matrix Q = T*AB*a
  Vector<t_complex> Q;
  //! The number of subdiagonals when distributing calculations
  t_int subdiagonals;
};
#endif
#endif
}
}

#endif
