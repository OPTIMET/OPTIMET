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

protected:
  std::shared_ptr<Geometry> geometry;        /**< Pointer to the geometry. */
  std::shared_ptr<Excitation const> incWave; /**< Pointer to the incoming excitation. */
  t_uint nMax;
};

/**
 * The Solver class builds and solves the scattering matrix equation.
 */
class Solver : public AbstractSolver {
public:
  using AbstractSolver::update;
#if defined(OPTIMET_BELOS)
  /**
   * Initialization constructor for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param method_ the solver method to be used.
   * @param nMax_ the maximum value for the n iterator.
   * @param belos_params parameters to setup the belos solvers
   * @param context Scalapack context associated with this solver instance
   */
  Solver(
      std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation const> incWave_, int method_,
      scalapack::Context const &context = scalapack::Context::Squarest(),
      Teuchos::RCP<Teuchos::ParameterList> belos_params = Teuchos::rcp(new Teuchos::ParameterList));
#else
  /**
   * Initialization constructor for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param method_ the solver method to be used.
   * @param nMax_ the maximum value for the n iterator.
   * @param context Scalapack context associated with this solver instance. Unused if compiled
   * without MPI.
   */
  Solver(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation const> incWave_,
         int method_, scalapack::Context const &context = scalapack::Context::Squarest());
#endif

  /**
   * Default destructor for the Solver class.
   */
  virtual ~Solver(){};

  /**
   * Switches the Solver object to the SH case.
   * Must have been created and initialized in a FF case.
   * @param incWave_ pointer to the new SH wave excitation.
   * @param result_FF_ pointer to the Fundamental Frequency result.
   * @param nMax_ the maximum value for the n iterator.
   * @return 0 if successful, 1 otherwise.
   */
  Solver &SH(Result *r);
  bool SH() const { return result_FF != nullptr; }

  /**
   * Solve the scattered and internal coefficients using the method specified by
   * solverMethod.
   * @param X_sca_ the return vector for the scattered coefficients.
   * @param X_int_ the return vector for the internal coefficients.
   * @return 0 if successful, 1 otherwise.
   */
  void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const override;
  //! \brief Solves linear system of equations
  //! \details Makes sure all procs in comm have access to result.
  //! The communicator should contain all the procs in the scalapack context of the solver.
  void solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_,
             mpi::Communicator const &comm) const override;
  //! \brief Update after internal parameters changed externally
  //! \details Because that's how the original implementation rocked.
  void update() override { populate(); }
  void
  update(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation const> incWave_) override;

  scalapack::Context context() const { return context_; }
  scalapack::Sizes const &block_size() const { return block_size_; }
  Solver &block_size(scalapack::Sizes const &c) {
    if(c.rows != c.cols)
      throw std::invalid_argument("ScaLAPACK solvers require a square block size");
    block_size_ = c;
    return *this;
  }

#ifdef OPTIMET_BELOS
  //! \brief Parameters for Belos/Trilinos solvers
  //! \note Mere access to the parameters requires the Teuchos::ParameterList to be modifiable. So
  //! the constness is not quite respected here.
  Teuchos::RCP<Teuchos::ParameterList> belos_parameters() const { return belos_params_; }
#endif

  //! Order of the scattering matrix
  t_uint scattering_size() const;

protected:
  //! Populate the S and Q matrices using the solverMethod option
  void populate();

  //! Populate the S and Q matrices using the Direct (Mischenko1996) method.
  void populateDirect();

  /**
   * Populate the S and Q matrices using the Indirect (Stout2002) method.
   * @return 0 if succesful, 1 otherwise.
   */
  void populateIndirect();
  //! Solves linear system of equations
  void solveLinearSystem(Matrix<t_complex> const &A, Vector<t_complex> const &b,
                         Vector<t_complex> &x, mpi::Communicator const &comm) const;
  //! Solves linear system of equations
  void solveLinearSystem(Matrix<t_complex> const &A, Vector<t_complex> const &b,
                         Vector<t_complex> &x) const {
    return solveLinearSystem(A, b, x, mpi::Communicator());
  }

private:
#ifdef OPTIMET_MPI
  void solveLinearSystemScalapack(Matrix<t_complex> const &A, Vector<t_complex> const &b,
                                  Vector<t_complex> &x, mpi::Communicator const &comm) const;
#endif
  //! Results for the fundamental frequency
  Result *result_FF;
  //! Solver method
  int solverMethod;
#ifdef OPTIMET_BELOS
  //! Parameter list of the belos solvers
  Teuchos::RCP<Teuchos::ParameterList> belos_params_;
#endif
  //! \brief MPI commnunicator
  //! \details Fake if not compiled with MPI
  scalapack::Context context_;
  scalapack::Sizes block_size_;

  Matrix<t_complex> S; /**< The scattering matrix S = I - T*AB. */
  Vector<t_complex> Q; /**< The local field matrix Q = T*AB*a. */
  long nMax;           /**< The maximum n order. */
};
}
} // namespace optimet
#endif /* SOLVER_H_ */
