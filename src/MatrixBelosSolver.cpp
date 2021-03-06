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
    int N = 2 * nMaxS * (nMaxS + 2);
    int nobj = geometry->objects.size();
  if(context().is_valid()) {
   auto const solver = belos_parameters()->get<std::string>("Solver");
    if(solver == "scalapack") {
      Scalapack::solve(X_sca_, X_int_, X_sca_SH, X_int_SH, CGcoeff);
      return;
    }
 
    int uppLIM =nobj * 2 * N;
    belos_parameters()->set("Maximum Iterations", uppLIM);
    belos_parameters()->set("Num Blocks", 1000);

    // FF part
    auto input = parallel_input();

    // Now the actual work
    auto start3 = high_resolution_clock::now();
    auto const gls_result = scalapack::gmres_linear_system(std::get<0>(input), std::get<1>(input),
                                                           belos_parameters(), splitcomm);
    auto stop3 = high_resolution_clock::now();
     auto duration3 = duration_cast<microseconds>(stop3 - start3);
   if(communicator().rank() == 0) {
     std::cout << "FF solve iterate-" << std::endl;
     std::cout << duration3.count()/1e6 <<"e-"<<  std::endl;
              }

    if(std::get<1>(gls_result) != 0)
      throw std::runtime_error("Error encountered while solving the linear system");
    // Transfer back to root
    X_sca_ = gather_all_source_vector(std::get<0>(gls_result));
    PreconditionedMatrix::unprecondition(X_sca_, X_int_);
   }
   if(context().size() != communicator().size()) {
    broadcast_to_out_of_context(X_sca_, context(), communicator());
    broadcast_to_out_of_context(X_int_, context(), communicator());
   }

  if(incWave->SH_cond){
  Vector<t_complex> KmNOD, K1;
 
  KmNOD = distributed_source_vector_SH_Mnode(*geometry, incWave, X_int_, X_sca_, CGcoeff);
  MPI_Barrier(MPI_COMM_WORLD);
  K1 = distributed_vector_SH_AR1(*geometry, incWave, X_int_, X_sca_, CGcoeff);
  MPI_Barrier(MPI_COMM_WORLD);

  if(context().is_valid()) {
  //SH part
   
     Vector<t_complex> K;
     int Dims = KmNOD.size();

     K = distributed_source_vector_SH(*geometry, KmNOD,  context(), block_size());
 
     auto input_SH = parallel_input_SH(K, Dims);

   // Now the actual work
    auto start4 = high_resolution_clock::now();
    auto const gls_result_SH =
    scalapack::gmres_linear_system(std::get<0>(input_SH), std::get<1>(input_SH),
                                                           belos_parameters(), splitcomm);
    auto stop4 = high_resolution_clock::now();
     auto duration4 = duration_cast<microseconds>(stop4 - start4);
   if(communicator().rank() == 0) {
     std::cout << "SH solve iterate-" << std::endl;
     std::cout << duration4.count()/1e6 <<"e-"<<  std::endl;
              }

    if(std::get<1>(gls_result_SH) != 0)
      throw std::runtime_error("Error encountered while solving the linear system");
    // Transfer back to root
    X_sca_SH = gather_all_source_vector(std::get<0>(gls_result_SH));

    PreconditionedMatrix::unprecondition_SH(X_sca_SH, X_int_SH, K1);


  }
  if(context().size() != communicator().size()) {

    broadcast_to_out_of_context(X_sca_SH, context(), communicator());
    broadcast_to_out_of_context(X_int_SH, context(), communicator());
   
  }
 }
}
}
}
