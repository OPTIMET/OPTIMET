#ifndef SOLVER_H_
#define SOLVER_H_

#include "Geometry.h"
#include "Excitation.h"
#include "Coupling.h"
#include "Result.h"
#include <complex>

using std::complex;

/**
 * The Solver class builds and solves the scattering matrix equation.
 */
class Solver
{
private:
  Geometry *geometry;   /**< Pointer to the geometry. */
  Excitation *incWave;  /**< Pointer to the incoming excitation. */
  Coupling AB;      /**< The coupling coefficients. */
  bool initDone;      /**< Specifies initialization status. */
  bool flagSH;      /**< Specifies if we have switched to the SH case. */
  long nMax;        /**< The maximum n order. */
  Result *result_FF;    /**< The fundamental frequency results. */
  int solverMethod;   /**< Solver method: Direct = Mischenko1996, Indirect = Stout2002 */
public:
  complex<double> **S;  /**< The scattering matrix S = I - T*AB. */
  complex<double> *Q;   /**< The local field matrix Q = T*AB*a. */

  /**
   * Default constructor for the Solver class.
   * Does NOT initialize the object.
   */
  Solver();

  /**
   * Initialization constructor for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param method_ the solver method to be used.
   * @param nMax_ the maximum value for the n iterator.
   */
  Solver(Geometry *geometry_, Excitation *incWave_, int method_, long nMax_);

  /**
   * Default destructor for the Solver class.
   */
  virtual ~Solver();

  /**
   * Initialization method for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param method_ the solver method to be used.
   * @param nMax_ the maximum value for the n iterator.
   */
  void init(Geometry *geometry_, Excitation *incWave_, int method_, long nMax_);

  /**
   * Switches the Solver object to the SH case.
   * Must have been created and initialized in a FF case.
   * @param incWave_ pointer to the new SH wave excitation.
   * @param result_FF_ pointer to the Fundamental Frequency result.
   * @param nMax_ the maximum value for the n iterator.
   * @return 0 if successful, 1 otherwise.
   */
  int switchSH(Excitation *incWave_, Result *result_FF_, long nMax_);

  /**
   * Populate the S and Q matrices using the solverMethod option.
   * @return 0 if successful, 1 otherwise.
   */
  int populate();

  /**
   * Populate the S and Q matrices using the Direct (Mischenko1996) method.
   * Default to Direct (Mischenko1996).
   * @return 0 if succesful, 1 otherwise.
   */
  int populateDirect();

  /**
   * Populate the S and Q matrices using the Indirect (Stout2002) method.
   * @return 0 if succesful, 1 otherwise.
   */
  int populateIndirect();

  /**
   * Solve the scattered and internal coefficients using the method specified by solverMethod.
   * @param X_sca_ the return vector for the scattered coefficients.
   * @param X_int_ the return vector for the internal coefficients.
   * @return 0 if successful, 1 otherwise.
   */
  int solve(complex<double> *X_sca_, complex<double> *X_int_);

  /**
   * Solve the S*X=Q equation using the Direct (Mischenko1996) method.
   * @param X_sca_ the return vector for the scattered coefficients.
   * @return 0 if successful, 1 otherwise.
   */
  int solveScatteredDirect(complex<double> *X_sca_);

  /**
   * Solve the S*X=Q equation using the Indirect (Stout2002) method.
   * @param X_sca_ the return vector for the scattered coefficients.
   * @return 0 if successful, 1 otherwise.
   */
  int solveScatteredIndirect(complex<double> *X_sca_);

  /**
   * Solve the internal coefficients.
   * @warning Must be called AFTER the solveScattered.
   * @param X_sca_ the input vector for the scattered coefficients.
   * @param X_int_ the return vector for the internal coefficients.
   * @return 0 if successful, 1 otherwise.
   */
  int solveInternal(complex<double> *X_sca_, complex<double> *X_int_);

  /**
   * Update method for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param nMax_ the maximum value for the n iterator.
   */
  void update(Geometry *geometry_, Excitation *incWave_, long nMax_);

};

#endif /* SOLVER_H_ */
