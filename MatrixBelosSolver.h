#ifndef OPTIMET_MATRIX_BELOS_SOLVER_H
#define OPTIMET_MATRIX_BELOS_SOLVER_H

#include "Types.h"

#ifdef OPTIMET_MPI
#ifdef OPTIMET_BELOS
#include "PreconditionedMatrix.h"
#include "PreconditionedMatrixSolver.h"
#include "ScalapackSolver.h"
#include "Solver.h"
#include "scalapack/context.h"
#include <Teuchos_ParameterList.hpp>

namespace optimet {
namespace solver {

//! Use an actual matrix, and Eigen's Householder QR method
class MatrixBelos : public Scalapack {
public:
  MatrixBelos(
      std::shared_ptr<Geometry> geometry, std::shared_ptr<Excitation const> incWave,
      mpi::Communicator const &communicator = mpi::Communicator(),
      scalapack::Context const &context = scalapack::Context::Squarest(),
      scalapack::Sizes const &block_size = scalapack::Sizes{64, 64},
      Teuchos::RCP<Teuchos::ParameterList> belos_params = Teuchos::rcp(new Teuchos::ParameterList))
      : Scalapack(geometry, incWave, communicator, context, block_size),
        belos_params_(belos_params) {}
  MatrixBelos(Run const &run)
      : MatrixBelos(run.geometry, run.excitation, run.communicator, run.context,
                    {run.parallel_params.block_size, run.parallel_params.block_size},
                    run.belos_params) {}

  void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const override;

  //! \brief Parameters for Belos/Trilinos solvers
  //! \note Mere access to the parameters requires the Teuchos::ParameterList to be modifiable. So
  //! the constness is not quite respected here.
  Teuchos::RCP<Teuchos::ParameterList> belos_parameters() const { return belos_params_; }

protected:
  //! Parameter list of the belos solvers
  Teuchos::RCP<Teuchos::ParameterList> belos_params_;
};
}
}
#endif
#endif
#endif
