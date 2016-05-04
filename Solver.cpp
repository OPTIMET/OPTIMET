#include "Solver.h"

#include "Algebra.h"
#include "Aliases.h"
#include "CompoundIterator.h"
#include "HarmonicsIterator.h"
#include "Tools.h"
#include "constants.h"
#include "mpi/Communicator.h"
#include "scalapack/LinearSystemSolver.h"
#include "scalapack/Matrix.h"
#include <Eigen/Dense>
#include <cstdlib>
#include <iostream>

namespace optimet {
#if defined(OPTIMET_BELOS)
Solver::Solver(Geometry *geometry, std::shared_ptr<Excitation const> incWave, int method, long nMax,
               Teuchos::RCP<Teuchos::ParameterList> belos_params, scalapack::Context const &context)
    : geometry(geometry), incWave(incWave), nMax(nMax), result_FF(nullptr), solverMethod(method),
      belos_params_(belos_params), context_(context), block_size_{64, 64}

{
  populate();
}
#elif defined(OPTIMET_MPI)
Solver::Solver(Geometry *geometry, std::shared_ptr<Excitation const> incWave, int method, long nMax,
               scalapack::Context const &context)
    : geometry(geometry), incWave(incWave), nMax(nMax), result_FF(nullptr), solverMethod(method),
      context_(context), block_size_{64, 64}

{
  populate();
}
#else
Solver::Solver(Geometry *geometry, std::shared_ptr<Excitation const> incWave, int method, long nMax)
    : geometry(geometry), incWave(incWave), nMax(nMax), result_FF(nullptr), solverMethod(method),
      block_size_{64, 64}

{
  populate();
}
#endif

void Solver::populate() {
  auto const N = 2 * (HarmonicsIterator::max_flat(nMax) - 1) * geometry->objects.size();
  S.resize(N, N);
  Q.resize(N);

  if(solverMethod == O3DSolverIndirect)
    populateIndirect();
  else // Default
    populateDirect();
}

void Solver::populateDirect() {
  auto const flatMax = HarmonicsIterator::max_flat(nMax) - 1;
  Matrix<t_complex> T_AB(2 * flatMax, 2 * flatMax);

  if(result_FF) // if SH simulation, set the local source first
    geometry->setSourcesSingle(incWave, result_FF->internal_coef.data(), nMax);

  for(size_t i = 0; i < geometry->objects.size(); i++) {
    Vector<t_complex> Q_local(2 * flatMax);

    // Get the T and IncLocal matrices first
    auto const T = geometry->getTLocal(incWave->omega, i, nMax);
    // we are in the SH case -> get the local sources from the geometry
    if(result_FF)
      geometry->getSourceLocal(i, incWave, result_FF->internal_coef.data(), nMax, Q_local.data());
    // we are in the FF case -> get the incoming excitation from the geometry
    else
      incWave->getIncLocal(geometry->objects[i].vR, Q_local.data(), nMax);
    Q.segment(i * 2 * flatMax, 2 * flatMax) = T * Q_local;

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

int Solver::solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) {
  if(not context().is_valid())
    throw std::runtime_error("Scalapack context is invalid");
  solveLinearSystem(S, Q, X_sca_);
  if(solverMethod == O3DSolverIndirect)
    X_sca_ = convertIndirect(X_sca_);

  X_int_ = solveInternal(X_sca_);
  return 0;
}

Vector<t_complex> Solver::convertIndirect(Vector<t_complex> const &scattered) {
  auto const N = 2 * (HarmonicsIterator::max_flat(nMax) - 1);
  Vector<t_complex> result(N * geometry->objects.size());
  for(size_t i = 0; i < geometry->objects.size(); i++)
    result.segment(i * N, N) =
        geometry->getTLocal(incWave->omega, i, nMax) * scattered.segment(i * N, N);
  return result;
}

Solver &Solver::SH(Result *r) {

  if(r != result_FF) {
    result_FF = r;
    populate();
  }

  return *this;
}

Vector<t_complex> Solver::solveInternal(Vector<t_complex> const &scattered) {
  auto const N = 2 * (HarmonicsIterator::max_flat(nMax) - 1);
  Vector<t_complex> result(scattered.size());
  // sets each result.segment(j * N, N) to something
  for(size_t j = 0; j < geometry->objects.size(); j++)
    geometry->getIaux(incWave->omega, j, nMax, result.data() + j * N);
  // then multiply by incomming array
  result.array() *= scattered.array();
  return result;
}

void Solver::populateIndirect() {
  auto const flatMax = HarmonicsIterator::max_flat(nMax) - 1;
  Matrix<t_complex> T_AB(2 * flatMax, 2 * flatMax);

  if(result_FF) // if SH simulation, set the local source first
    geometry->setSourcesSingle(incWave, result_FF->internal_coef.data(), nMax);

  for(size_t i = 0; i < geometry->objects.size(); i++) {
    Vector<t_complex> Q_local(2 * flatMax);

    // Get the IncLocal matrices
    // These correspond directly to the Beta*a in Stout2002 Eq. 10 as
    // they are already translated.
    // we are in the SH case -> get the local sources from the geometry
    if(result_FF)
      geometry->getSourceLocal(i, incWave, result_FF->internal_coef.data(), nMax, Q_local.data());
    // we are in the FF case -> get the incoming excitation from the geometry
    else
      incWave->getIncLocal(geometry->objects[i].vR, Q_local.data(), nMax);

    // Push Q_local into Q (direct equivalence)
    Q.segment(i * 2 * flatMax, 2 * flatMax) = Q_local;

    for(size_t j = 0; j < geometry->objects.size(); j++)
      if(i == j)
        S.block(i * 2 * flatMax, i * 2 * flatMax, 2 * flatMax, 2 * flatMax) =
            Matrix<>::Identity(2 * flatMax, 2 * flatMax);
      else {
        // Build the T_AB matrix (non-regular corresponding to alpha(i,j))
        Coupling const AB(geometry->objects[i].vR - geometry->objects[j].vR, incWave->waveK, nMax);

        T_AB.topLeftCorner(flatMax, flatMax) = AB.diagonal.transpose();
        T_AB.bottomRightCorner(flatMax, flatMax) = AB.diagonal.transpose();
        T_AB.topRightCorner(flatMax, flatMax) = AB.offdiagonal.transpose();
        T_AB.bottomLeftCorner(flatMax, flatMax) = AB.offdiagonal.transpose();

        S.block(i * 2 * flatMax, j * 2 * flatMax, 2 * flatMax, 2 * flatMax) =
            -T_AB * geometry->getTLocal(incWave->omega, j, nMax);
      }
  }
}

void Solver::update(Geometry *geometry_, std::shared_ptr<Excitation const> incWave_, long nMax_) {
  geometry = geometry_;
  incWave = incWave_;
  nMax = nMax_;

  result_FF = nullptr;

  populate();
}

void Solver::solveLinearSystem(Matrix<t_complex> const &A, Vector<t_complex> const &b,
                               Vector<t_complex> &x) const {
#ifdef OPTIMET_MPI
  if(scalapack::global_size() > 1)
    solveLinearSystemScalapack(A, b, x);
  else
#endif
    x = A.colPivHouseholderQr().solve(b);
}

#ifdef OPTIMET_MPI
void Solver::solveLinearSystemScalapack(Matrix<t_complex> const &A, Vector<t_complex> const &b,
                                        Vector<t_complex> &x) const {
  // scalapack parameters: matrix size, grid of processors, block size
  scalapack::Sizes const size{static_cast<t_uint>(A.rows()), static_cast<t_uint>(A.cols())};

  // "serial" version to distribute from A to root.
  scalapack::Matrix<t_complex> Aserial(context().serial(), size, block_size());
  scalapack::Matrix<t_complex> bserial(Aserial.context(), {size.rows, 1}, block_size());
  if(Aserial.context().is_valid()) {
    Aserial.local() = A;
    bserial.local() = b;
  }
  // Transfer to grid
  auto Aparallel = Aserial.transfer_to(context(), block_size());
  auto bparallel = bserial.transfer_to(context(), block_size());

  // Now the actual work
  auto Xparallel = std::get<0>(
#ifdef OPTIMET_BELOS
      belos_parameters()->get<std::string>("Solver", "scalapack") != "scalapack" ?
          scalapack::gmres_linear_system(Aparallel, bparallel, belos_parameters()) :
#endif
          scalapack::general_linear_system(Aparallel, bparallel));

  // Transfer back to root
  auto Xserial = Xparallel.transfer_to(Aserial.context(), block_size());
  // Broadcast from root
  x = context().broadcast(Xserial.local(), 0, 0);
}
#endif

} // optimet namespace
