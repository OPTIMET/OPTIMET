#ifndef OPTIMET_SCALAPACK_SOLVER_H
#define OPTIMET_SCALAPACK_SOLVER_H

#include "Types.h"

#ifdef OPTIMET_MPI
#include "PreconditionedMatrix.h"
#include "PreconditionedMatrixSolver.h"
#include "Solver.h"
#include "scalapack/context.h"

namespace optimet {
namespace solver {

//! Use an actual matrix, and Eigen's Householder QR method
class Scalapack : public PreconditionedMatrix {
public:
  Scalapack(std::shared_ptr<Geometry> geometry, std::shared_ptr<Excitation const> incWave,
            scalapack::Context const &context = scalapack::Context::Squarest(),
            scalapack::Sizes const &block_size = scalapack::Sizes{64, 64})
      : PreconditionedMatrix(geometry, incWave), context_(context), block_size_(block_size) {
    update();
  }

  void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_,
             mpi::Communicator const &comm) const override;
  void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const override;
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
};
}
}
#endif
#endif
