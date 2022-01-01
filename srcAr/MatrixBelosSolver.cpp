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

#include "MatrixBelosSolver.h"
#include "scalapack/LinearSystemSolver.h"
#include "scalapack/BroadcastToOutOfContext.h"
#include <chrono>
#include <Eigen/Dense>
using namespace std::chrono;
namespace optimet {
namespace solver {

void MatrixBelos::solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_,Vector<t_complex> &X_sca_SH,
                         Vector<t_complex> &X_int_SH, std::vector<double *> CGcoeff) const {
  //FF
  if(context().is_valid()) {
   
    Matrix<t_complex> TmatrixFF, RgQmatrixFF;
    int nMax = geometry->nMax();
    int pMax = nMax * (nMax + 2);
    TmatrixFF = S.block(0 ,0 , 2*pMax, 2*pMax);
    RgQmatrixFF = S.block(0 , 2*pMax , 2*pMax, 2*pMax);

    X_sca_ = Q;

    PreconditionedMatrix::unprecondition(X_sca_, X_int_, TmatrixFF, RgQmatrixFF);
 }
  //SH
  if(incWave->SH_cond){
  Vector<t_complex> KmNOD, K1;
  Matrix<t_complex> TmatrixSH, RgQmatrixSH;

   int nMaxS = geometry->nMaxS();
   int pMax = nMaxS * (nMaxS + 2);
  TmatrixSH = V.block(0 ,0 , 2*pMax, 2*pMax);
  RgQmatrixSH = V.block(0 , 2*pMax , 2*pMax, 2*pMax);

  KmNOD = distributed_source_vector_SH_Mnode(*geometry, incWave, X_int_, X_sca_, TmatrixSH);
  MPI_Barrier(MPI_COMM_WORLD);

  K1 =  distributed_vector_SH_AR1(*geometry, incWave, X_int_, X_sca_);
  MPI_Barrier(MPI_COMM_WORLD);

  if(context().is_valid()) {
  
    X_sca_SH = KmNOD;

    PreconditionedMatrix::unprecondition_SH(X_sca_SH, X_int_SH, K1, RgQmatrixSH);
 
   
  }
  
 }
}
}
}
