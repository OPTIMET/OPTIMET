#ifndef SOLVER_H_
#define SOLVER_H_

#include "Coupling.h"
#include "Excitation.h"
#include "Geometry.h"
#include "Result.h"
#include "Types.h"
#include "scalapack/Context.h"
#include "scalapack/Parameters.h"
#include <complex>
#include <complex>
#include <exception>

#ifdef OPTIMET_BELOS
#include <Teuchos_ParameterList.hpp>
#endif

namespace optimet {
/**
 * The Solver class builds and solves the scattering matrix equation.
 */
class Solver {
public:
  Matrix<t_complex> S; /**< The scattering matrix S = I - T*AB. */
  Vector<t_complex> Q; /**< The local field matrix Q = T*AB*a. */

#if defined(OPTIMET_MPI)
  /**
   * Initialization constructor for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param method_ the solver method to be used.
   * @param nMax_ the maximum value for the n iterator.
   * @param context Scalapack context associated with this solver instance
   */
  Solver(Geometry *geometry_, Excitation const *incWave_, int method_, long nMax_,
         scalapack::Context const &context = scalapack::Context::Squarest());
#elif defined(OPTIMET_BELOS)
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
      Geometry *geometry_, Excitation const *incWave_, int method_, long nMax_,
      Teuchos::RCP<Teuchos::ParameterList> belos_params = Teuchos::rcp(new Teuchos::ParameterList),
      scalapack::Context const &context = scalapack::Context::Squarest());
#else
  /**
   * Initialization constructor for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param method_ the solver method to be used.
   * @param nMax_ the maximum value for the n iterator.
   */
  Solver(Geometry *geometry_, Excitation const *incWave_, int method_, long nMax_);
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
  int solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_);
  /**
   * Update method for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param nMax_ the maximum value for the n iterator.
   */
  void update(Geometry *geometry_, Excitation const *incWave_, long nMax_);
  //! \brief Update after internal parameters changed externally
  //! \details Because that's how the original implementation rocked.
  void update() { populate(); }

  //! Converts back to the scattered result from the indirect calculation
  Vector<t_complex> convertIndirect(Vector<t_complex> const &scattered);

  //! Solves for the internal coefficients.
  Vector<t_complex> solveInternal(Vector<t_complex> const &X_sca_);

  scalapack::Context context() const { return context_; }
  Solver &context(scalapack::Context const &c) {
    context_ = c;
    return *this;
  }
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
                         Vector<t_complex> &x) const;

private:
#ifdef OPTIMET_MPI
  void solveLinearSystemScalapack(Matrix<t_complex> const &A, Vector<t_complex> const &b,
                                  Vector<t_complex> &x) const;
#endif
  Geometry *geometry;        /**< Pointer to the geometry. */
  Excitation const *incWave; /**< Pointer to the incoming excitation. */
  long nMax;                 /**< The maximum n order. */
  Result *result_FF;         /**< The fundamental frequency results. */
  int solverMethod;          /**< Solver method: Direct = Mischenko1996, Indirect =
                                Stout2002 */
#ifdef OPTIMET_BELOS
  Teuchos::RCP<Teuchos::ParameterList> belos_params_;
#endif
  //! \brief MPI commnunicator
  //! \details Fake if not compiled with MPI
  scalapack::Context context_;
  scalapack::Sizes block_size_;
};
} // namespace optimet
#endif /* SOLVER_H_ */
