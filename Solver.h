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
class SolverBase {
public:
  /**
   * Initialization constructor for the Solver class.
   * @param geometry_ the geometry of the simulation.
   * @param incWave_ the incoming wave excitation.
   * @param method_ the solver method to be used.
   * @param nMax_ the maximum value for the n iterator.
   */
  SolverBase(std::shared_ptr<Geometry> geometry, std::shared_ptr<Excitation const> incWave,
             long nMax)
      : geometry(geometry), incWave(incWave), nMax(nMax) {}

  ~SolverBase(){};

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
  virtual void update(std::shared_ptr<Geometry> geometry_,
                      std::shared_ptr<Excitation const> incWave_, long nMax_) {
    geometry = geometry_;
    incWave = incWave_;
    nMax = nMax_;
    update();
  }
  void update(Run const &run) { return update(run.geometry, run.excitation, run.nMax); }
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
  long nMax;                                 /**< The maximum n order. */
};

/**
 * The Solver class builds and solves the scattering matrix equation.
 */
class Solver : public SolverBase {
public:
  using SolverBase::update;
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
      long nMax_, scalapack::Context const &context = scalapack::Context::Squarest(),
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
         int method_, long nMax_,
         scalapack::Context const &context = scalapack::Context::Squarest());
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
  void update(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation const> incWave_,
              long nMax_) override;

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
};

//! \brief Computes source vector
Vector<t_complex>
source_vector(Geometry const &geometry, std::shared_ptr<Excitation const> incWave);
//! \brief Computes source vector from fundamental frequency
Vector<t_complex> local_source_vector(Geometry const &geometry,
                                      std::shared_ptr<Excitation const> incWave,
                                      Vector<t_complex> const &input_coeffs);

//! Computes preconditioned scattering matrix
Matrix<t_complex> preconditioned_scattering_matrix(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave);

//! Computes preconditioned scattering matrix in paralllel
Matrix<t_complex> preconditioned_scattering_matrix(Geometry const &geometry,
                                                   std::shared_ptr<Excitation const> incWave,
                                                   scalapack::Context const &context,
                                                   scalapack::Sizes const &blocks);
//! Distributes the source vectors
Vector<t_complex> distributed_source_vector(Vector<t_complex> const &input,
                                            scalapack::Context const &context,
                                            scalapack::Sizes const &blocks);
#ifdef OPTIMET_MPI
//! Gather the distributed vector into a single vector
Vector<t_complex> gather_all_source_vector(t_uint n, Vector<t_complex> const &input,
                                           scalapack::Context const &context,
                                           scalapack::Sizes const &blocks);
//! Gather the distributed vector into a single vector
Vector<t_complex> gather_all_source_vector(scalapack::Matrix<t_complex> const &matrix);

//! \brief Broadcast data from a proc in the context to procs outside the context
//! \details Usefull if some procs are not part of the context but still require the data.
template <class T>
void broadcast_to_out_of_context(T &inout, scalapack::Context const &context,
                                 mpi::Communicator const &comm) {
  auto const is_in_context = comm.all_gather<int>(context.is_valid());
  auto is_true = [](int input) { return input; };
  auto is_false = [](int input) { return not input; };
  // All procs in context, nothing to do
  if(std::all_of(is_in_context.begin(), is_in_context.end(), is_true))
    return;
  bool const is_root =
      is_in_context[comm.rank()] and
      std::all_of(is_in_context.begin(), is_in_context.begin() + comm.rank(), is_false);
  auto const is_in_group = is_root or not is_in_context[comm.rank()];
  auto const split = comm.split(is_in_group, is_root ? 0 : 1);
  if(is_in_group)
    inout = split.broadcast(inout, 0);
}
#else
//! Broadcast data (vectors, matrices from a proc in the context to procs outside the context)
template <class T>
void broadcast_to_out_of_context(T &, scalapack::Context const &, mpi::Communicator const &) {}
#endif
} // namespace optimet
#endif /* SOLVER_H_ */
