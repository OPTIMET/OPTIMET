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

#ifndef SOLVER_H_
#define SOLVER_H_

#include "Coupling.h"
#include "Excitation.h"
#include "Geometry.h"
#include "Result.h"
#include "Run.h"
#include "Types.h"
#include "mpi/Collectives.h"
#include "mpi/Communicator.h"
#include "scalapack/Context.h"
#include "scalapack/Matrix.h"
#include "scalapack/Parameters.h"
#include <complex>
#include <exception>
#include <memory>

#ifdef OPTIMET_BELOS
#include <Teuchos_ParameterList.hpp>
#endif

namespace optimet {
//! Computes coeffs scattered from spheres FF
Vector<t_complex> convertIndirect(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground, std::vector<Scatterer> const &);
                                  
                                  
//! Computes coeffs internal to spheres FF
Vector<t_complex> convertInternal(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground, std::vector<Scatterer> const &);
  
  
  //! Computes coeffs scattered from spheres SH
Vector<t_complex> convertIndirect_SH_outer(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground, std::vector<Scatterer> const &); 
                                  
//  Computes coeffs internal to spheres SH                                
Vector<t_complex> convertInternal_SH(Vector<t_complex> const &scattered, Vector<t_complex> const &K_1, Vector<t_complex> const &K_1ana, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &);                                  
                                  
                                                                
                                  
                                  
namespace solver {

//! Abstract base class for all solvers
class AbstractSolver {
public:
  /**
   * Initialization constructor for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param method_ the solver method to be used.
   */
  AbstractSolver(std::shared_ptr<Geometry> geometry, std::shared_ptr<Excitation const> incWave,
                 mpi::Communicator const &communicator = mpi::Communicator())
      : geometry(geometry), incWave(incWave), communicator_(communicator), nMax(0) {}

  AbstractSolver(Run const &run) : AbstractSolver(run.geometry, run.excitation, run.communicator) {}


  ~AbstractSolver(){};

  /**
   * Solve the scattered and internal coefficients using the method specified by
   * solverMethod.
   * @param X_sca_ the return vector for the scattered coefficients.
   * @param X_int_ the return vector for the internal coefficients.
   * @return 0 if successful, 1 otherwise.
   */
  virtual void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_, Vector<t_complex> &X_sca_SH, 
                    Vector<t_complex> &X_int_SH, std::vector<double *> CGcoeff) const = 0;
  /**
   * Update method for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param nMax_ the maximum value for the n iterator.
   */
  virtual void
  update(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation const> incWave_) {
    geometry = geometry_;
    incWave = incWave_;
    update();
  }
  void update(Run const &run) { return update(run.geometry, run.excitation); }
  //! \brief Update after internal parameters changed externally
  //! \details Because that's how the original implementation rocked.
  virtual void update() = 0;

  //! Converts back to the scattered result from the indirect calculation
  Vector<t_complex> convertIndirect(Vector<t_complex> const &scattered) const {
    return optimet::convertIndirect(scattered, incWave->omega(), geometry->bground,
                                    geometry->objects);
  }
  

  //! Solves for the internal coefficients fundamental frequency
  Vector<t_complex> solveInternal(Vector<t_complex> const &scattered) const {
    return optimet::convertInternal(scattered, incWave->omega(), geometry->bground,
                                    geometry->objects);
  }
  
  
    //! Converts back to the scattered result from the indirect calculation SH frequency outer coefficients
  Vector<t_complex> convertIndirect_SH_outer(Vector<t_complex> const &scattered) const {
    return optimet::convertIndirect_SH_outer(scattered, incWave->omega(), geometry->bground,
                                    geometry->objects);
  }
  
      //! Solves for the internal coefficients SH frequency
  Vector<t_complex> solveInternal_SH(Vector<t_complex> const &scattered, Vector<t_complex> const &K_1, Vector<t_complex> const &K_1ana) const {
    return optimet::convertInternal_SH(scattered, K_1, K_1ana, incWave->omega(), geometry->bground,
                                    geometry->objects);
  }

  
  
  //! Number of spherical harmonics in expansion
  t_uint scattering_size() const {
    auto const n = nMax == 0 ? geometry->nMax() : nMax;
    return 2 * n * (n + 2) * geometry->objects.size();
  }

  mpi::Communicator const &communicator() const { return communicator_; }

protected:
  std::shared_ptr<Geometry> geometry;        /**< Pointer to the geometry. */
  std::shared_ptr<Excitation const> incWave; /**< Pointer to the incoming excitation. */
  mpi::Communicator communicator_;
  t_uint nMax;
};

//! A factory function for solvers
std::shared_ptr<AbstractSolver> factory(Run const &run);
}
} // namespace optimet
#endif /* SOLVER_H_ */
