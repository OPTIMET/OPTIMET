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

Vector<t_complex> convertInternal(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {
  Vector<t_complex> result(scattered.size());
  size_t i = 0;
  auto const N = 2 * objects[0].nMax * (objects[0].nMax + 2);
  Matrix<t_complex> Intrmatrix(N, N);

  for(auto const &object : objects) {
  
   if (object.scatterer_type == "sphere"){
    result.segment(i, N).array() =
        scattered.segment(i, N).array() * object.getIaux(omega, bground).array();
    }

    i += N;

  }
   
  return result;
}

Vector<t_complex> convertIndirect(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {
 
  return scattered;
}


Vector<t_complex> convertIndirect_SH_outer(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {
  
   return scattered;
}


Vector<t_complex> convertInternal_SH(Vector<t_complex> const &scattered, Vector<t_complex> const &K_1ana, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {
  Vector<t_complex> result(scattered.size());
  auto const N = 2 * objects[0].nMaxS * (objects[0].nMaxS + 2);
  Matrix<t_complex> Intrmatrix(N, N);
  size_t i = 0; 
 auto const k_b_SH = 2 * omega * std::sqrt(bground.epsilon * bground.mu);

  if (objects[0].scatterer_type == "sphere"){
  for(auto const &object : objects) {
    
    result.segment(i, N).array() =
        scattered.segment(i, N).array() * object.getIauxSH1(omega, bground).array();
        
    result.segment(i, N) = result.segment(i, N)  - K_1ana.segment(i, N);  
  
    i += N;
  }
 } 
  return result;
}

} // optimet namespace
