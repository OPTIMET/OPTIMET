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

#include "Solver.h"
#include "ElectroMagnetic.h"
#include "MatrixBelosSolver.h"
#include "PreconditionedMatrixSolver.h"
#include "ScalapackSolver.h"
#include "Scatterer.h"
#include "Solver.h"
#include "Types.h"

namespace optimet {
namespace solver {


std::shared_ptr<AbstractSolver> factory(Run const &run) {

#ifndef OPTIMET_MPI

  return std::make_shared<PreconditionedMatrix>(run);
#elif defined(OPTIMET_SCALAPACK) && !defined(OPTIMET_BELOS)
  
  return std::make_shared<Scalapack>(run);
#elif defined(OPTIMET_BELOS) && defined(OPTIMET_SCALAPACK)
   
  if(run.belos_params()->get<std::string>("Solver") == "eigen")
    return std::make_shared<PreconditionedMatrix>(run);
  if(run.belos_params()->get<std::string>("Solver") == "scalapack")
    return std::make_shared<Scalapack>(run);
  
  return std::make_shared<MatrixBelos>(run);
#elif defined(OPTIMET_BELOS)

  if(run.belos_params()->get<std::string>("Solver") == "scalapack")
    throw std::runtime_error("Optimet was not compiled with scalapack");
  
#else
#error Need at least Belos to run MPI solvers
#endif
}
} // namespace solver

Vector<t_complex> convertInternal(Vector<t_complex> const &scattered, Matrix<t_complex> const &RgQ, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {
  Vector<t_complex> result(scattered.size());
  auto const nobj = objects.size();
  auto const N = 2 * objects[0].nMax * (objects[0].nMax + 2);
  Matrix<t_complex> Qmatsingle;

  for (int ii = 0; ii != nobj; ++ii ){
    
    Qmatsingle = RgQ.block(0, ii*N, N, N);

    result.segment(ii*N, N) = Qmatsingle.inverse() * scattered.segment(ii*N, N);
}   
  return result;
}



Vector<t_complex> convertIndirect(Vector<t_complex> const &scattered, Matrix<t_complex> const &Tmat, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {
  Vector<t_complex> result(scattered.size());
  auto const nobj = objects.size();
  auto const N = 2 * objects[0].nMax * (objects[0].nMax + 2);
  Matrix<t_complex> Tmatsingle;
  
 for (int ii = 0; ii != nobj; ++ii ){

    Tmatsingle = Tmat.block(0, ii*N, N, N);
    
    result.segment(ii*N, N) = Tmatsingle * scattered.segment(ii*N, N);
     
 }

  return result;
}


Vector<t_complex> convertIndirect_SH_outer(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {
  
   return scattered;
}


Vector<t_complex> convertInternal_SH(Vector<t_complex> const &scattered, Vector<t_complex> const &K_1, Matrix<t_complex> const &RgQ, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {

  Vector<t_complex> result(scattered.size());
  auto const N = 2 * objects[0].nMaxS * (objects[0].nMaxS + 2);
  auto const nobj = objects.size();
  auto const k_b_SH = 2 * omega * std::sqrt(bground.epsilon * bground.mu);
 
  Matrix<t_complex> Qmatsingle;

  for (int ii = 0; ii != nobj; ++ii ){

    Qmatsingle = RgQ.block(0, ii*N, N, N);
   
    result.segment(ii*N, N) = Qmatsingle.inverse()*((-consCi/k_b_SH)*scattered.segment(ii*N, N) + K_1.segment(ii*N, N));
  
 } 
   
  return result;
}

} // optimet namespace
