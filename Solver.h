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
//! Computes coeffs scattered from spheres
Vector<t_complex> convertIndirect(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground, std::vector<Scatterer> const &);
//! Computes coeffs internal to spheres
Vector<t_complex> convertInternal(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground, std::vector<Scatterer> const &);
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
  AbstractSolver(std::shared_ptr<Geometry> geometry, std::shared_ptr<Excitation const> incWave)
      : geometry(geometry), incWave(incWave) {}

  ~AbstractSolver(){};

  /**
   * Solve the scattered and internal coefficients using the method specified by
   * solverMethod.
   * @param X_sca_ the return vector for the scattered coefficients.
   * @param X_int_ the return vector for the internal coefficients.
   * @return 0 if successful, 1 otherwise.
   */
  virtual void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const = 0;
  //! \brief Solves linear system of equations
  //! \details Makes sure all procs in comm have access to result.
  //! The communicator should contain all the procs in the scalapack context of the solver.
  virtual void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_,
                     mpi::Communicator const &comm) const = 0;
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

  //! Solves for the internal coefficients.
  Vector<t_complex> solveInternal(Vector<t_complex> const &scattered) const {
    return optimet::convertInternal(scattered, incWave->omega(), geometry->bground,
                                    geometry->objects);
  }

  //! Number of spherical harmonics in expansion
  t_uint scattering_size() const { return 2 * nMax * (nMax + 2) * geometry->objects.size(); }

protected:
  std::shared_ptr<Geometry> geometry;        /**< Pointer to the geometry. */
  std::shared_ptr<Excitation const> incWave; /**< Pointer to the incoming excitation. */
  t_uint nMax;
};

//! A factory function for solvers
std::shared_ptr<AbstractSolver> factory(Run const &run);
}
} // namespace optimet
#endif /* SOLVER_H_ */
