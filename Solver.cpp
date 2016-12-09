#include "Solver.h"

#include "Algebra.h"
#include "Aliases.h"
#include "CompoundIterator.h"
#include "HarmonicsIterator.h"
#include "PreconditionedMatrix.h"
#include "Tools.h"
#include "constants.h"
#include "mpi/Communicator.h"
#include "scalapack/BroadcastToOutOfContext.h"
#include "scalapack/LinearSystemSolver.h"
#include "scalapack/Matrix.h"
#include <Eigen/Dense>
#include <cstdlib>
#include <iostream>

namespace optimet {
namespace solver {
#if defined(OPTIMET_BELOS)
Solver::Solver(std::shared_ptr<Geometry> geometry, std::shared_ptr<Excitation const> incWave,
               int method, scalapack::Context const &context,
               Teuchos::RCP<Teuchos::ParameterList> belos_params)
    : AbstractSolver(geometry, incWave), result_FF(nullptr), solverMethod(method),
      belos_params_(belos_params), context_(context), block_size_{64, 64}, nMax(0)

{
  populate();
}
#else
Solver::Solver(std::shared_ptr<Geometry> geometry, std::shared_ptr<Excitation const> incWave,
               int method, scalapack::Context const &context)
    : AbstractSolver(geometry, incWave), result_FF(nullptr), solverMethod(method),
      context_(context), block_size_{64, 64}

{
  populate();
}
#endif

void Solver::populate() {
  nMax = geometry->nMax();
  assert(solverMethod == O3DSolverIndirect);
  populateIndirect();
}

void Solver::populateDirect() {
  auto const flatMax = HarmonicsIterator::max_flat(nMax) - 1;
  Matrix<t_complex> T_AB(2 * flatMax, 2 * flatMax);

  if(result_FF) // if SH simulation, set the local source first
    geometry->setSourcesSingle(incWave, result_FF->internal_coef.data(), nMax);

  for(size_t i = 0; i < geometry->objects.size(); i++) {
    Vector<t_complex> Q_local(2 * flatMax);

    // Get the T and IncLocal matrices first
    auto const T = geometry->getTLocal(incWave->omega(), i, nMax);
    // we are in the SH case -> get the local sources from the geometry
    if(result_FF)
      geometry->getSourceLocal(i, incWave, nMax, Q_local.data());
    // we are in the FF case -> get the incoming excitation from the geometry
    else
      incWave->getIncLocal(geometry->objects[i].vR, Q_local.data(), nMax);
    Q.segment(i * 2 * flatMax, 2 * flatMax) = T.array() * Q_local.array();

    for(size_t j = 0; j < geometry->objects.size(); j++)
      if(i == j)
        S.block(i * 2 * flatMax, i * 2 * flatMax, 2 * flatMax, 2 * flatMax) =
            Matrix<t_complex>::Identity(2 * flatMax, 2 * flatMax);
      else {
        // Build the T_AB matrix
        Coupling const AB(geometry->objects[i].vR - geometry->objects[j].vR, incWave->waveK, nMax);

        T_AB.topLeftCorner(flatMax, flatMax) = AB.diagonal.transpose();
        T_AB.bottomRightCorner(flatMax, flatMax) = AB.diagonal.transpose();
        T_AB.topRightCorner(flatMax, flatMax) = AB.offdiagonal.transpose();
        T_AB.bottomLeftCorner(flatMax, flatMax) = AB.offdiagonal.transpose();

        S.block(i * 2 * flatMax, j * 2 * flatMax, 2 * flatMax, 2 * flatMax) = -T * T_AB;
      }
  }
}

void Solver::solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const {
  // If the context is invalid, we cannot ensure that this proc will receive the solution.
  // This context is the only one we have, and what we need is a context where the solution exists
  // or can be computed. We do not have that context.
  if(not context().is_valid())
    throw std::runtime_error("Scalapack context is invalid");
  solveLinearSystem(S, Q, X_sca_);
  if(solverMethod == O3DSolverIndirect)
    X_sca_ = convertIndirect(X_sca_);

  X_int_ = solveInternal(X_sca_);
}

void Solver::solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_,
                   mpi::Communicator const &comm) const {
  assert(comm.size() >= context().size());
#ifdef OPTIMET_MPI
  auto const splitcomm = comm.split(context().is_valid());
#else
  auto const &splitcomm = comm;
#endif
  if(context().is_valid()) {
    solveLinearSystem(S, Q, X_sca_, splitcomm);
    if(solverMethod == O3DSolverIndirect)
      X_sca_ = convertIndirect(X_sca_);

    X_int_ = solveInternal(X_sca_);
  }
  broadcast_to_out_of_context(X_sca_, context(), comm);
  broadcast_to_out_of_context(X_int_, context(), comm);
}

Solver &Solver::SH(Result *r) {

  if(r != result_FF) {
    result_FF = r;
    populate();
  }

  return *this;
}

void Solver::populateIndirect() {
  Q = source_vector(*geometry, incWave);
#ifdef OPTIMET_MPI
  Q = distributed_source_vector(Q, context(), block_size());
#endif

  S = preconditioned_scattering_matrix(*geometry, incWave, context(), block_size());
}

void Solver::update(std::shared_ptr<Geometry> geometry_,
                    std::shared_ptr<Excitation const> incWave_) {
  result_FF = nullptr;
  AbstractSolver::update(geometry_, incWave_);
}

void Solver::solveLinearSystem(Matrix<t_complex> const &A, Vector<t_complex> const &b,
                               Vector<t_complex> &x, mpi::Communicator const &comm) const {
#ifdef OPTIMET_BELOS
  if(belos_parameters()->get<std::string>("Solver", "scalapack") != "eigen")
    solveLinearSystemScalapack(A, b, x, comm);
  else
#elif defined(OPTIMET_MPI)
  if(comm.size() > 1)
    solveLinearSystemScalapack(A, b, x, comm);
  else
#endif
    x = A.colPivHouseholderQr().solve(b);
}

#ifdef OPTIMET_MPI
template <class SCALARA, class SCALARB>
typename scalapack::Matrix<SCALARA>::ConcreteMatrix
solveLinearSystem(Solver const &solver, scalapack::Matrix<SCALARA> const &A,
                  scalapack::Matrix<SCALARB> const &b, mpi::Communicator const &comm) {
#ifdef OPTIMET_BELOS
  auto const result =
      solver.belos_parameters()->get<std::string>("Solver", "scalapack") != "scalapack" ?
          scalapack::gmres_linear_system(A, b, solver.belos_parameters(), comm) :
          scalapack::general_linear_system(A, b);
#else
  auto const result = scalapack::general_linear_system(A, b);
#endif
  if(std::get<1>(result) != 0)
    throw std::runtime_error("Error encountered while solving the linear system");
  return std::get<0>(result);
}

void Solver::solveLinearSystemScalapack(Matrix<t_complex> const &A, Vector<t_complex> const &b,
                                        Vector<t_complex> &x, mpi::Communicator const &comm) const {
  auto const N = 2 * (HarmonicsIterator::max_flat(nMax) - 1) * geometry->objects.size();
  scalapack::Matrix<t_complex> Aparallel(context(), {N, N}, block_size());
  if(Aparallel.size() > 0)
    Aparallel.local() = A;
  scalapack::Matrix<t_complex> bparallel(context(), {N, 1}, block_size());
  if(bparallel.local().size() > 0)
    bparallel.local() = b;

  // Now the actual work
  auto Xparallel = optimet::solver::solveLinearSystem(*this, Aparallel, bparallel, comm);
  // Transfer back to root
  x = gather_all_source_vector(Xparallel);
}
#endif

t_uint Solver::scattering_size() const {
  return 2 * (HarmonicsIterator::max_flat(nMax) - 1) * geometry->objects.size();
}
} // namespace solver

Vector<t_complex> convertInternal(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {
  Vector<t_complex> result(scattered.size());
  size_t i = 0;
  for(auto const &object : objects) {
    auto const N = 2 * object.nMax * (object.nMax + 2);
    result.segment(i, N).array() =
        scattered.segment(i, N).array() * object.getIaux(omega, bground).array();
    i += N;
  }
  return result;
}

Vector<t_complex> convertIndirect(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {
  Vector<t_complex> result(scattered.size());
  size_t i(0);
  for(auto const &object : objects) {
    auto const N = 2 * object.nMax * (object.nMax + 2);
    result.segment(i, N) =
        object.getTLocal(omega, bground).array() * scattered.segment(i, N).array();
    i += N;
  }
  return result;
}
} // optimet namespace
