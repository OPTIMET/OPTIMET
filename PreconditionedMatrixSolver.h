#ifndef OPTIMET_PRECONDITIONNED_MATRIX_SOLVER_H
#define OPTIMET_PRECONDITIONNED_MATRIX_SOLVER_H

#include "PreconditionedMatrix.h"
#include "Solver.h"
#include "Types.h"
#include <Eigen/Dense>

namespace optimet {
namespace solver {

//! Use an actual matrix, and Eigen's Householder QR method
class PreconditionedMatrix : public AbstractSolver {
public:
  PreconditionedMatrix(std::shared_ptr<Geometry> geometry,
                       std::shared_ptr<Excitation const> incWave)
      : AbstractSolver(geometry, incWave) {
    update();
  }

  //! Ignores communicator and solves serially
  void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_,
             mpi::Communicator const &) const override {
    return solve(X_sca_, X_int_);
  }

  void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const override {
    X_sca_ = S.colPivHouseholderQr().solve(Q);
    unprecondition(X_sca_, X_int_);
  }

  void update() override {
    Q = source_vector(*geometry, incWave);
    S = preconditioned_scattering_matrix(*geometry, incWave);
  }

protected:
  //! The scattering matrix S = I - T*AB
  Matrix<t_complex> S;
  //! The local field matrix Q = T*AB*a
  Vector<t_complex> Q;

  //! Unpreconditions the result of preconditioned computation
  void unprecondition(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const {
    X_sca_ = AbstractSolver::convertIndirect(X_sca_);
    X_int_ = AbstractSolver::solveInternal(X_sca_);
  }
};
}
}
#endif
