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
  auto const N = geometry->scatterer_size();
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

void Solver::solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const {
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
  auto const splitcomm = comm.split(context().is_valid());
  if(context().is_valid()) {
    solveLinearSystem(S, Q, X_sca_, splitcomm);
    if(solverMethod == O3DSolverIndirect)
      X_sca_ = convertIndirect(X_sca_);

    X_int_ = solveInternal(X_sca_);
  }
  broadcast_to_out_of_context(X_sca_, context(), comm);
  broadcast_to_out_of_context(X_int_, context(), comm);
}

Vector<t_complex> Solver::convertIndirect(Vector<t_complex> const &scattered) const {
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

Vector<t_complex> Solver::solveInternal(Vector<t_complex> const &scattered) const {
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
  }
  assert(geometry->objects.size() == 0 or nMax == geometry->objects.front().nMax);
  S = preconditioned_scattering_matrix(*geometry, *incWave);
}

void Solver::update(Geometry *geometry_, std::shared_ptr<Excitation const> incWave_, long nMax_) {
  geometry = geometry_;
  incWave = incWave_;
  nMax = nMax_;

  result_FF = nullptr;

  populate();
}

void Solver::solveLinearSystem(Matrix<t_complex> const &A, Vector<t_complex> const &b,
                               Vector<t_complex> &x, mpi::Communicator const &comm) const {
#ifdef OPTIMET_BELOS
  if(belos_parameters()->get<std::string>("Solver", "scalapack") != "eigen")
    solveLinearSystemScalapack(A, b, x, comm);
  else
#elif OPTIMET_MPI
  if(scalapack::global_size() > 1)
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
  // scalapack parameters: matrix size, grid of processors, block size
  scalapack::Sizes const size{static_cast<t_uint>(A.rows()), static_cast<t_uint>(A.cols())};

  // "serial" version to distribute from A to root.
  auto const serial_context = context().serial();
  auto const nrows = serial_context.is_valid() ? A.rows() : 0;
  auto const ncols = serial_context.is_valid() ? A.cols() : 0;
  assert(A.rows() == b.rows());
  Eigen::Map<Matrix<t_complex> const> const amap(A.data(), nrows, ncols);
  Eigen::Map<Matrix<t_complex> const> const bmap(b.data(), nrows, ncols != 0 ? 1 : 0);
  scalapack::Matrix<t_complex const *> Aserial(amap, serial_context, size, block_size());
  scalapack::Matrix<t_complex const *> bserial(bmap, Aserial.context(), {size.rows, 1},
                                               block_size());
  if(context().size() == 1) {
    auto const X = ::optimet::solveLinearSystem(*this, Aserial, bserial, comm);
    x = context().broadcast(X.local(), 0, 0);
  } else {
    // Transfer to grid
    auto Aparallel = Aserial.transfer_to(context(), block_size());
    auto bparallel = bserial.transfer_to(context(), block_size());
    auto const X = ::optimet::solveLinearSystem(*this, Aparallel, bparallel, comm);
    auto Xserial = X.transfer_to(Aserial.context(), block_size());
    x = context().broadcast(Xserial.local(), 0, 0);
  }
}
#endif

Matrix<t_complex>
preconditioned_scattering_matrix(std::vector<Scatterer>::const_iterator const &first,
                                 std::vector<Scatterer>::const_iterator const &end_first,
                                 std::vector<Scatterer>::const_iterator const &second,
                                 std::vector<Scatterer>::const_iterator const &end_second,
                                 ElectroMagnetic const &bground, Excitation const &incWave) {
  auto const nMax = first->nMax;
  auto const n = HarmonicsIterator::max_flat(nMax) - 1;

  Matrix<t_complex> result =
      Matrix<t_complex>::Identity(2 * n * (end_first - first), 2 * n * (end_second - second));
  size_t y(0);
  for(auto iterj(second); iterj != end_second; ++iterj, y += 2 * n) {
    Matrix<t_complex> const factor = -iterj->getTLocal(incWave.omega, bground);
    size_t x(0);
    for(auto iteri(first); iteri != end_first; ++iteri, x += 2 * n) {
      if(iteri == iterj)
        continue;
      Coupling const AB(iteri->vR - iterj->vR, incWave.waveK, nMax);
      result.block(x, y, n, n) = AB.diagonal.transpose();
      result.block(x + n, y + n, n, n) = AB.diagonal.transpose();
      result.block(x, y + n, n, n) = AB.offdiagonal.transpose();
      result.block(x + n, y, n, n) = AB.offdiagonal.transpose();
      result.block(x, y, 2 * n, 2 * n) *= factor;
    }
  }
  return result;
}

Matrix<t_complex> preconditioned_scattering_matrix(std::vector<Scatterer> const &objects,
                                                   ElectroMagnetic const &bground,
                                                   Excitation const &incWave) {
  return preconditioned_scattering_matrix(objects.begin(), objects.end(), objects.begin(),
                                          objects.end(), bground, incWave);
}

Matrix<t_complex>
preconditioned_scattering_matrix(Geometry const &geometry, Excitation const &incWave) {
  if(geometry.objects.size() == 0)
    return Matrix<t_complex>(0, 0);
  // Check nMax is same accross all objects
  auto const nMax = geometry.objects.front().nMax;
  for(auto const &scatterer : geometry.objects)
    if(scatterer.nMax != nMax)
      throw std::runtime_error("All objects must have same number of harmonics");
  return preconditioned_scattering_matrix(geometry.objects, geometry.bground, incWave);
}

#ifdef OPTIMET_MPI
scalapack::Matrix<Scalar> preconditioned_scattering_matrix(Geometry const &geometry,
                                                           Excitation const &incWave,
                                                           scalapack::Context const &context) {}
#endif
} // optimet namespace
