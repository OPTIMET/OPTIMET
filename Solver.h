#ifndef SOLVER_H_
#define SOLVER_H_

#include <complex>
#include "Types.h"
#include "Geometry.h"
#include "Excitation.h"
#include "Coupling.h"
#include "Result.h"
#include "mpi/Communicator.h"
#include "scalapack/Parameters.h"

namespace optimet {
/**
 * The Solver class builds and solves the scattering matrix equation.
 */
class Solver {
public:
  Matrix<t_complex> S; /**< The scattering matrix S = I - T*AB. */
  Vector<t_complex> Q; /**< The local field matrix Q = T*AB*a. */

  /**
   * Initialization constructor for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param method_ the solver method to be used.
   * @param nMax_ the maximum value for the n iterator.
   */
  Solver(Geometry *geometry_, Excitation const *incWave_, int method_, long nMax_,
         mpi::Communicator const &communicator = mpi::Communicator());

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

protected:
  /**
   * Populate the S and Q matrices using the solverMethod option.
    @return 0 if successful, 1 otherwise.
   */
  int populate();

  /**
   * Populate the S and Q matrices using the Direct (Mischenko1996) method.
   * Default to Direct (Mischenko1996).
   * @return 0 if succesful, 1 otherwise.
   */
  void populateDirect();
  int populateDirectOld();

  /**
   * Populate the S and Q matrices using the Indirect (Stout2002) method.
   * @return 0 if succesful, 1 otherwise.
   */
  void populateIndirect();
  /**
   * Populate the S and Q matrices using the Indirect (Stout2002) method.
   * @return 0 if succesful, 1 otherwise.
   */
  int populateIndirectOld();

  /**
   * Solve the S*X=Q equation using the Direct (Mischenko1996) method.
   * @param X_sca_ the return vector for the scattered coefficients.
   * @return 0 if successful, 1 otherwise.
   */
  int solveScatteredDirect(Vector<t_complex> &X_sca_);

  /**
   * Solve the S*X=Q equation using the Indirect (Stout2002) method.
   * @param X_sca_ the return vector for the scattered coefficients.
   * @return 0 if successful, 1 otherwise.
   */
  int solveScatteredIndirect(Vector<t_complex> &X_sca_);

  /**
   * Solve the internal coefficients.
   * @warning Must be called AFTER the solveScattered.
   * @param X_sca_ the input vector for the scattered coefficients.
   * @param X_int_ the return vector for the internal coefficients.
   * @return 0 if successful, 1 otherwise.
   */
  int solveInternal(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_);

  mpi::Communicator const &communicator() const { return communicator_; }
  Solver &communicator(mpi::Communicator const &c) {
    communicator_ = c;
    return *this;
  }
  scalapack::Parameters const &parallel_params() const { return parallel_params_; }
  Solver &parallel_params(scalapack::Parameters const &c) {
    parallel_params_ = c;
    return *this;
  }

  //! Solves linear system of equations
  void solveLinearSystem(Matrix<t_complex> const &A, Vector<t_complex> const &b,
                         Vector<t_complex> &x) const;

private:
  Geometry *geometry;        /**< Pointer to the geometry. */
  Excitation const *incWave; /**< Pointer to the incoming excitation. */
  long nMax;                 /**< The maximum n order. */
  Result *result_FF;         /**< The fundamental frequency results. */
  int solverMethod;          /**< Solver method: Direct = Mischenko1996, Indirect =
                                Stout2002 */
  //! \brief MPI commnunicator
  //! \details Fake if not compiled with MPI
  mpi::Communicator communicator_;
  scalapack::Parameters parallel_params_;
};
} // namespace optimet
#endif /* SOLVER_H_ */
