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
    
    // parameters for ACA-gmres solver
    double tol = 1e-6;
    int maxit = 250;
    int no_rest = 3;
    // FF part
    Vector<t_complex> Q;
    
    if (geometry->ACA_cond_){
    Q = source_vector(*geometry, incWave);
    X_sca_ = Gmres_Zcomp(S_comp_FF, Q, tol, maxit, no_rest, *geometry);
    PreconditionedMatrix::unprecondition(X_sca_, X_int_);
    }
    else {
    if(context().is_valid()) {
    auto input = parallel_input();
    // Now the actual work
    auto const gls_result =
        scalapack::general_linear_system(std::get<0>(input), std::get<1>(input));
    
    if(std::get<1>(gls_result) != 0)
      throw std::runtime_error("Error encountered while solving the linear system");
    // Transfer back to root
    X_sca_ = gather_all_source_vector(std::get<0>(gls_result)); 
    PreconditionedMatrix::unprecondition(X_sca_, X_int_);

  }
  }
  if(context().size() != communicator().size()) {
    broadcast_to_out_of_context(X_sca_, context(), communicator());
    broadcast_to_out_of_context(X_int_, context(), communicator());
  }

  if(incWave->SH_cond){
  Vector<t_complex> KmNOD, K1, K1ana, X_int_conj;
  X_int_conj = X_int_.conjugate();
  
  KmNOD = distributed_source_vector_SH_Mnode(*geometry, incWave, X_int_conj, X_sca_, CGcoeff);
  MPI_Barrier(MPI_COMM_WORLD);
  
  K1ana = source_vectorSH_K1ana_parallel(*geometry, incWave, X_int_conj, X_sca_, CGcoeff);
  MPI_Barrier(MPI_COMM_WORLD);

  if (geometry->ACA_cond_){
  X_sca_SH = Gmres_Zcomp(S_comp_SH, KmNOD, tol, maxit, no_rest, *geometry);
  PreconditionedMatrix::unprecondition_SH(X_sca_SH, X_int_SH, K1ana);
  }
  else{   
  if(context().is_valid()) {
    //SH part
    int Dims = KmNOD.size();
    Vector<t_complex> K;
 
    K = distributed_source_vector_SH(*geometry, KmNOD, context(), block_size());
     auto input_SH = parallel_input_SH(K, Dims);
    // Now the actual work
    auto const gls_result_SH =
        scalapack::general_linear_system(std::get<0>(input_SH), std::get<1>(input_SH));

    if(std::get<1>(gls_result_SH) != 0)
      throw std::runtime_error("Error encountered while solving the linear system");
    // Transfer back to root
    X_sca_SH = gather_all_source_vector(std::get<0>(gls_result_SH));
       
    PreconditionedMatrix::unprecondition_SH(X_sca_SH, X_int_SH, K1ana);
    
  }
}
if(context().size() != communicator().size()) {
    broadcast_to_out_of_context(X_sca_SH, context(), communicator());
    broadcast_to_out_of_context(X_int_SH, context(), communicator());
  }
}
 
}

void Scalapack::update() {

  Q = distributed_source_vector(source_vector(*geometry, incWave), context(), block_size());
  MPI_Barrier(MPI_COMM_WORLD);  
  
  if (geometry->ACA_cond_)
  Scattering_matrix_ACA_FF_parallel(*geometry, incWave, S_comp_FF);
  else
  S = preconditioned_scattering_matrix(*geometry, incWave, context(), block_size());
  MPI_Barrier(MPI_COMM_WORLD);

 if(incWave->SH_cond){
  
  if (geometry->ACA_cond_)
    Scattering_matrix_ACA_SH_parallel(*geometry, incWave, S_comp_SH);
  else
  V = preconditioned_scattering_matrix_SH(*geometry, incWave, context(), block_size());
  MPI_Barrier(MPI_COMM_WORLD);
 
 }
 
}
}
}
