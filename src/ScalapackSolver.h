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

#ifndef OPTIMET_SCALAPACK_SOLVER_H
#define OPTIMET_SCALAPACK_SOLVER_H

#include "Types.h"

#ifdef OPTIMET_SCALAPACK
#include "PreconditionedMatrix.h"
#include "PreconditionedMatrixSolver.h"
#include "Solver.h"
#include "scalapack/Context.h"

namespace optimet {
namespace solver {

//! Use an actual matrix, and Eigen's Householder QR method
class Scalapack : public PreconditionedMatrix {
public:
  Scalapack(std::shared_ptr<Geometry> geometry, std::shared_ptr<Excitation const> incWave,
            mpi::Communicator const &comm = mpi::Communicator(),
            scalapack::Context const &context = scalapack::Context::Squarest(),
            scalapack::Sizes const &block_size = scalapack::Sizes{64, 64})
      :PreconditionedMatrix(geometry, incWave, comm), context_(context), block_size_(block_size) {
    update();
  }
  Scalapack(Run const &run)
      : Scalapack(run.geometry, run.excitation, run.communicator, run.context,
                  {run.parallel_params.block_size, run.parallel_params.block_size}) {}

  void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_SH, 
             Vector<t_complex> &X_int_SH, std::vector<double *> CGcoeff) const override;
  void update() override;

  //! Scalapack context used during computation
  scalapack::Context context() const { return context_; }
  //! Block size into which to separate the problem
  scalapack::Sizes const &block_size() const { return block_size_; }

protected:
  //! Scalapack context to use during communication
  scalapack::Context context_;
  //! Block-sizes for scalapack
  scalapack::Sizes block_size_;

  //! Creates scalapack matrix wrappers
  std::tuple<scalapack::Matrix<t_complex>, scalapack::Matrix<t_complex>> parallel_input() const;
  std::tuple<scalapack::Matrix<t_complex>, scalapack::Matrix<t_complex>> parallel_input_SH(Vector<t_complex> &K, int Dims) const;
};
}
}
#endif
#endif
