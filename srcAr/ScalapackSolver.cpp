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

#include "ScalapackSolver.h"
#include "scalapack/BroadcastToOutOfContext.h"
#include "scalapack/LinearSystemSolver.h"
#include <chrono>
#include <Eigen/Dense>
using namespace std::chrono;
namespace optimet {
namespace solver {

std::tuple<scalapack::Matrix<t_complex>, scalapack::Matrix<t_complex>>
Scalapack::parallel_input() const {
  auto const nMax = geometry->nMax();
  auto const N = 2 * nMax * (nMax + 2) * geometry->objects.size();
  scalapack::Matrix<t_complex> Aparallel(context(), {N, N}, block_size());
  if(Aparallel.size() > 0)
    Aparallel.local() = S;
  scalapack::Matrix<t_complex> bparallel(context(), {N, 1}, block_size());
  if(bparallel.local().size() > 0)
    bparallel.local() = Q;
  return std::make_tuple(Aparallel, bparallel);
}

std::tuple<scalapack::Matrix<t_complex>, scalapack::Matrix<t_complex>>
Scalapack::parallel_input_SH(Vector<t_complex> &K, int Dims) const {

   t_uint N = Dims;

  scalapack::Matrix<t_complex> Aparallel(context(), {N, N}, block_size());
  if(Aparallel.size() > 0)
    Aparallel.local() = V;
  scalapack::Matrix<t_complex> bparallel(context(), {N, 1}, block_size());
  if(bparallel.local().size() > 0)
    bparallel.local() = K;
  return std::make_tuple(Aparallel, bparallel);
}


void Scalapack::solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_SH,
                      Vector<t_complex> &X_int_SH, std::vector<double *> CGcoeff) const {
    
 if(context().is_valid()) {
    //FF
     auto const nobj = geometry->objects.size();
    Matrix<t_complex> TmatrixFF, RgQmatrixFF, SCATmatFF;
    int nMax = geometry->nMax();
    int pMax = nMax * (nMax + 2);
    TmatrixFF = S.block(0 ,0 , 2*pMax, nobj*2*pMax);
    RgQmatrixFF = S.block(0 , nobj*2*pMax , 2*pMax, nobj*2*pMax);

    SCATmatFF = ScatteringMatrixFF(TmatrixFF, *geometry, incWave);
    
    X_sca_ = SCATmatFF.colPivHouseholderQr().solve(Q);   
    PreconditionedMatrix::unprecondition(X_sca_, X_int_, TmatrixFF, RgQmatrixFF);

  }

  if(context().size() != communicator().size()) {
    broadcast_to_out_of_context(X_sca_, context(), communicator());
    broadcast_to_out_of_context(X_int_, context(), communicator());
  }

  //SH
  if(incWave->SH_cond){
  auto const nobj = geometry->objects.size();
  Vector<t_complex> KmNOD, K1;
  Matrix<t_complex> TmatrixSH, RgQmatrixSH, SCATmatSH;
  int nMaxS = geometry->nMaxS();
  int pMax = nMaxS * (nMaxS + 2);
  TmatrixSH = V.block(0 ,0 , 2*pMax, nobj*2*pMax);
  RgQmatrixSH = V.block(0 , nobj*2*pMax , 2*pMax, nobj*2*pMax);
   
  KmNOD = distributed_source_vector_SH_Mnode(*geometry, incWave, X_int_, X_sca_, TmatrixSH);
 
  MPI_Barrier(MPI_COMM_WORLD);
  K1 =  distributed_vector_SH_AR1(*geometry, incWave, X_int_, X_sca_);
  MPI_Barrier(MPI_COMM_WORLD);
      
  if(context().is_valid()) {
    
    SCATmatSH = ScatteringMatrixSH(TmatrixSH, *geometry, incWave);
    X_sca_SH = SCATmatSH.colPivHouseholderQr().solve(KmNOD);

    PreconditionedMatrix::unprecondition_SH(X_sca_SH, X_int_SH, K1, RgQmatrixSH);

  }

if(context().size() != communicator().size()) {
    broadcast_to_out_of_context(X_sca_SH, context(), communicator());
    broadcast_to_out_of_context(X_int_SH, context(), communicator());
  }

}
 
}

void Scalapack::update() {

  Q = source_vector(*geometry, incWave);
  MPI_Barrier(MPI_COMM_WORLD);
  
  S = getTRgQmatrix_FF_parr(*geometry, incWave);
  MPI_Barrier(MPI_COMM_WORLD);

 if(incWave->SH_cond){
  V = getTRgQmatrix_SH_parr(*geometry, incWave);
  MPI_Barrier(MPI_COMM_WORLD);
 }


}
}
}
