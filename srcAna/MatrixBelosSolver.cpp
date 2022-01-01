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

  if(belos_parameters()->get<std::string>("Solver", "GMRES") == "scalapack")
    return Scalapack::solve(X_sca_, X_int_, X_sca_SH, X_int_SH, CGcoeff);
  auto const splitcomm = communicator().split(context().is_valid());
    int nMaxS = geometry->nMaxS();
    int N = nMaxS * (nMaxS + 2);
    int nobj = geometry->objects.size();
    double tol = 1e-6;
    int maxit = 340;
    int no_rest = 1;
    Vector<t_complex> Q;

    if (geometry->ACA_cond_){
    Q = source_vector(*geometry, incWave);
    X_sca_ = Gmres_Zcomp(S_comp_FF, Q, tol, maxit, no_rest, *geometry);
    PreconditionedMatrix::unprecondition(X_sca_, X_int_);
    }
  else{
  if(context().is_valid()) {
   auto const solver = belos_parameters()->get<std::string>("Solver");
    if(solver == "scalapack") {
      Scalapack::solve(X_sca_, X_int_, X_sca_SH, X_int_SH, CGcoeff);
      return;
    }

    // FF part
    auto input = parallel_input();

    // Now the actual work
    auto const gls_result = scalapack::gmres_linear_system(std::get<0>(input), std::get<1>(input),
                                                           belos_parameters(), splitcomm);

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
  Vector<t_complex> KmNOD, K1, X_int_conj, K1ana;
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
   
     Vector<t_complex> K;
     int Dims = KmNOD.size();

     K = distributed_source_vector_SH(*geometry, KmNOD,  context(), block_size());
 
     auto input_SH = parallel_input_SH(K, Dims);

   // Now the actual work
    auto const gls_result_SH =
    scalapack::gmres_linear_system(std::get<0>(input_SH), std::get<1>(input_SH),
                                                           belos_parameters(), splitcomm);

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
}
}
