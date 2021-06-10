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
                       std::shared_ptr<Excitation const> incWave,
                       mpi::Communicator const &communicator = mpi::Communicator())
      : AbstractSolver(geometry, incWave, communicator) {
  #ifndef OPTIMET_MPI
    update();
   #endif
  }


  PreconditionedMatrix(Run const &run)
      : PreconditionedMatrix(run.geometry, run.excitation, run.communicator) {}

  void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_SH,

              Vector<t_complex> &X_int_SH, std::vector<double *> CGcoeff) const override {
  
    // fundamental frequency
    
    X_sca_ = S.colPivHouseholderQr().solve(Q);
    
    unprecondition(X_sca_, X_int_);

    Vector<t_complex> K, K1;
    
    // SH frequency
    if(incWave->SH_cond){

    K = source_vectorSH(*geometry, incWave, X_int_, X_sca_, CGcoeff);

    K1 = source_vectorSHarb1(*geometry, incWave, X_int_, X_sca_, CGcoeff);

    X_sca_SH = V.colPivHouseholderQr().solve(K);
    
    unprecondition_SH(X_sca_SH, X_int_SH, K1);
    }
  }
  
    

  void update() override {
  
    Q = source_vector(*geometry, incWave);

    S = preconditioned_scattering_matrix(*geometry, incWave);
    if(incWave->SH_cond)
    V = preconditioned_scattering_matrixSH(*geometry, incWave);

  }
  

protected:
  //! The scattering matrix S = I - T*AB
  Matrix<t_complex> S;
  //! The local field matrix Q = AB*a sources fundamental frequency
  Vector<t_complex> Q;
  // Second harmonic scattering matrix, outer coefficients
  Matrix<t_complex> V;
  //! The local field matrix K= AB*a sources SH frequency
  Vector<t_complex> K;

  //! Unpreconditions the result of preconditioned computation
  void unprecondition(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const {
    X_sca_ = AbstractSolver::convertIndirect(X_sca_);
    X_int_ = AbstractSolver::solveInternal(X_sca_);
    }
     
    
    void unprecondition_SH(Vector<t_complex> &X_sca_SH, Vector<t_complex> &X_int_SH, Vector<t_complex> &K1) const {
    X_sca_SH = AbstractSolver::convertIndirect_SH_outer(X_sca_SH);
    X_int_SH = AbstractSolver::solveInternal_SH(X_sca_SH, K1);
    
    }

   
};
}
}

#endif
