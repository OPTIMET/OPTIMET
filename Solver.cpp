#include "Solver.h"

#include <iostream>
#include <cstdlib>
#include "CompoundIterator.h"
#include "Tools.h"
#include "Algebra.h"
#include "constants.h"
#include "Aliases.h"
#include "mpi/Communicator.h"
#include "scalapack/Matrix.h"
#include <Eigen/Dense>

Solver::Solver(Geometry *geometry_, Excitation *incWave_, int method_, long nMax_,
               optimet::mpi::Communicator const &c)
    : flagSH(false), communicator_(c) {
  geometry = geometry_;
  incWave = incWave_;
  nMax = nMax_;

  flagSH = false;

  solverMethod = method_;

  CompoundIterator p;

  S.resize(2 * p.max(nMax) * geometry->objects.size(), 2 * p.max(nMax) * geometry->objects.size());
  Q.resize(2 * p.max(nMax) * geometry->objects.size());

  result_FF = NULL;

  populate();
}

int Solver::populate() {
  if(solverMethod == O3DSolverDirect)
    populateDirect();
  else if(solverMethod == O3DSolverIndirect)
    populateIndirect();
  else // Default
    populateDirect();

  return 0;
}


int Solver::populateDirect() {
  using namespace optimet;

  CompoundIterator p;
  CompoundIterator q;

  // Local matrices for storing T, T_fin = -T*A and T_AB = AB to be multiplied
  // by BLAS.
  // Also local matrix Q_local to be added to S.

  std::complex<double> **T = new std::complex<double> *[2 * p.max(nMax)];
  std::complex<double> **Ti = new std::complex<double> *[2 * p.max(nMax)];
  std::complex<double> **T_fin = new std::complex<double> *[2 * p.max(nMax)];
  std::complex<double> **T_AB = new std::complex<double> *[2 * (p.max(nMax))];

  for(p = 0; p < (int)(2 * p.max(nMax)); p++) {
    T[p] = new std::complex<double>[2 * p.max(nMax)];
    Ti[p] = new std::complex<double>[2 * p.max(nMax)];
    T_fin[p] = new std::complex<double>[2 * p.max(nMax)];
    T_AB[p] = new std::complex<double>[2 * p.max(nMax)];
  }

  if(flagSH) // if SH simulation, set the local source first
  {
    geometry->setSourcesSingle(incWave, result_FF->internal_coef, nMax);
  }

  for(size_t i = 0; i < geometry->objects.size(); i++) {
    int pMax = p.max(nMax);
    int qMax = q.max(nMax);

    std::complex<double> *Q_local = new std::complex<double>[2 * pMax];

    // Get the T and IncLocal matrices first
    geometry->getTLocal(incWave->omega, i, nMax, T);
    if(flagSH) // we are in the SH case -> get the local sources from the
               // geometry
    {
      geometry->getSourceLocal(i, incWave, result_FF->internal_coef, nMax, Q_local);
    } else // we are in the FF case -> get the incoming excitation from the
           // geometry
    {
      incWave->getIncLocal(geometry->objects[i].vR, Q_local, nMax);
    }

    // Allocate space for a "result" matrix C, multiply T [A] with Q_local and
    // push to Q.
    std::complex<double> *C = new std::complex<double>[2 * pMax];
    Algebra::multiplyVectorMatrix(T, 2 * pMax, 2 * pMax, Q_local, C, consC1, consC0);

    for(size_t j = 0; j < static_cast<size_t>(2 * pMax); j++) {
      Q(j + i * 2 * pMax) = C[j];
    }

    delete[] Q_local;
    delete[] C;

    for(size_t j = 0; j < geometry->objects.size(); j++) {
      if(i == j) {
        // Build and push an I matrix
        Tools::makeUnitMatrix(2 * pMax, Ti);
        pushToMatrix(Ti, 2 * pMax, 2 * pMax, S, i * 2 * pMax, i * 2 * pMax);
      } else {
        // Build the T_AB matrix
        Coupling const AB(geometry->objects[i].vR - geometry->objects[j].vR, incWave->waveK, nMax);

        for(p = 0; p < pMax; p++)
          for(q = 0; q < qMax; q++) {
            T_AB[p][q] = AB.diagonal(q, p);
            T_AB[p + pMax][q + qMax] = AB.diagonal(q, p);
            T_AB[p + pMax][q] = AB.offdiagonal(q, p);
            T_AB[p][q + qMax] = AB.offdiagonal(q, p);
          }

        // Multiply T and T_AB to get T_fin = -T*T_AB
        Algebra::multiplyMatrixMatrix(T, 2 * pMax, 2 * pMax, T_AB, 2 * pMax, 2 * pMax, T_fin,
                                      consCm1, consC0);

        // Push T into the S matrix
        pushToMatrix(T_fin, 2 * pMax, 2 * pMax, S, i * 2 * pMax, j * 2 * pMax);
      }
    }
  }

  // Clean up memory

  for(p = 0; p < (int)(2 * p.max(nMax)); p++) {
    delete[] T[p];
    delete[] Ti[p];
    delete[] T_fin[p];
    delete[] T_AB[p];
  }

  delete[] T;
  delete[] Ti;
  delete[] T_fin;
  delete[] T_AB;

  return 0;
}

int Solver::solve(std::complex<double> *X_sca_, std::complex<double> *X_int_) {
  if(solverMethod == O3DSolverDirect)
    solveScatteredDirect(X_sca_);
  else if(solverMethod == O3DSolverIndirect)
    solveScatteredIndirect(X_sca_);
  else // Default
    solveScatteredDirect(X_sca_);

  solveInternal(X_sca_, X_int_);

  return 0;
}

int Solver::solveScatteredDirect(std::complex<double> *X_sca_) {
  using namespace optimet;

  solveLinearSystem(S, Q, X_sca_);

  return 0;
}

int Solver::solveScatteredIndirect(std::complex<double> *X_sca_) {
  using namespace optimet;

  CompoundIterator p;

  std::complex<double> *X_sca_local =
      new std::complex<double>[2 * p.max(nMax) * geometry->objects.size()];
  std::complex<double> *X_sca_part =
      new std::complex<double>[2 * p.max(nMax) * geometry->objects.size()];
  std::complex<double> *X_sca_part_fin =
      new std::complex<double>[2 * p.max(nMax) * geometry->objects.size()];
  std::complex<double> **T = new std::complex<double> *[2 * p.max(nMax)];

  for(p = 0; p < (int)(2 * p.max(nMax)); p++) {
    T[p] = new std::complex<double>[2 * p.max(nMax)];
  }

  // Solve the equation, here Q and S correspond to Eq. 10 in (Stout2002). Store
  // result in X_sca_local
  solveLinearSystem(S, Q, X_sca_local);

  for(size_t i = 0; i < geometry->objects.size(); i++) {
    // Get the local scattering matrix for object i
    geometry->getTLocal(incWave->omega, i, nMax, T);

    // Build a partial X_sca_part vector
    for(p = 0; p < (int)(2 * p.max(nMax)); p++) {
      X_sca_part[p] = X_sca_local[i * 2 * p.max(nMax) + p.compound];
    }

    // Multiply T with X_sca_part to get X_sca_fin
    Algebra::multiplyVectorMatrix(T, 2 * p.max(nMax), 2 * p.max(nMax), X_sca_part, X_sca_part_fin,
                                  consC1, consC0);

    // Push new X_sca_fin into the X_sca_ final solution
    for(p = 0; p < (int)(2 * p.max(nMax)); p++) {
      X_sca_[i * 2 * p.max(nMax) + p.compound] = X_sca_part_fin[p];
    }
  }

  // Clean up the memory
  delete[] X_sca_local;
  delete[] X_sca_part;
  delete[] X_sca_part_fin;
  for(p = 0; p < (int)(2 * p.max(nMax)); p++) {
    delete[] T[p];
  }
  delete[] T;

  return 0;
}

int Solver::switchSH(Excitation *incWave_, Result *result_FF_, long nMax_) {

  incWave = incWave_;
  nMax = nMax_;

  result_FF = result_FF_;

  flagSH = true;

  return 0;
}

int Solver::solveInternal(std::complex<double> *X_sca_, std::complex<double> *X_int_) {

  CompoundIterator p, q;

  std::complex<double> *Iaux = new std::complex<double>[2 * Tools::iteratorMax(nMax)];

  for(size_t j = 0; j < geometry->objects.size(); j++) {
    int pMax = p.max(nMax);
    geometry->getIaux(incWave->omega, j, nMax, Iaux);

    for(p = 0; p < pMax; p++) {
      // scattered
      X_int_[j * 2 * pMax + p.compound] = X_sca_[j * 2 * pMax + p.compound] * Iaux[p];
      X_int_[pMax + j * 2 * pMax + p.compound] =
          X_sca_[pMax + j * 2 * pMax + p.compound] * Iaux[p + pMax];
    }
  }

  delete[] Iaux;

  return 0;
}

int Solver::populateIndirect() {
  using namespace optimet;

  CompoundIterator p;
  CompoundIterator q;

  // Local matrices for storing T, T_fin = -T*A and T_AB = AB to be multiplied
  // by BLAS.
  // Also local matrix Q_local to be added to S.

  std::complex<double> **T = new std::complex<double> *[2 * p.max(nMax)];
  std::complex<double> **Ti = new std::complex<double> *[2 * p.max(nMax)];
  std::complex<double> **T_fin = new std::complex<double> *[2 * p.max(nMax)];
  std::complex<double> **T_AB = new std::complex<double> *[2 * (p.max(nMax))];

  for(p = 0; p < (int)(2 * p.max(nMax)); p++) {
    T[p] = new std::complex<double>[2 * p.max(nMax)];
    Ti[p] = new std::complex<double>[2 * p.max(nMax)];
    T_fin[p] = new std::complex<double>[2 * p.max(nMax)];
    T_AB[p] = new std::complex<double>[2 * p.max(nMax)];
  }

  if(flagSH) // if SH simulation, set the local source first
  {
    geometry->setSourcesSingle(incWave, result_FF->internal_coef, nMax);
  }

  for(size_t i = 0; i < geometry->objects.size(); i++) {
    int const pMax = p.max(nMax);
    int const qMax = q.max(nMax);

    Vector<t_complex> Q_local(2 * pMax);

    // Get the IncLocal matrices
    // These correspond directly to the Beta*a in Stout2002 Eq. 10 as
    // they are already translated.
    if(flagSH) // we are in the SH case -> get the local sources from the
               // geometry
    {
      geometry->getSourceLocal(i, incWave, result_FF->internal_coef, nMax, Q_local.data());
    } else // we are in the FF case -> get the incoming excitation from the
           // geometry
    {
      incWave->getIncLocal(geometry->objects[i].vR, Q_local.data(), nMax);
    }

    // Push Q_local into Q (direct equivalence) (NOTE: j is unused at this
    // point)
    Q.segment(i * 2 * pMax, 2 * pMax) = Q_local;

    for(size_t j = 0; j < geometry->objects.size(); j++) {
      if(i == j)
        S.block(i * 2 * pMax, i * 2 * pMax, 2 * pMax, 2 * pMax) =
            Matrix<>::Identity(2 * pMax, 2 * pMax);
      else {
        // Build the T_AB matrix (non-regular corresponding to alpha(i,j))
        Coupling const AB(geometry->objects[i].vR - geometry->objects[j].vR, incWave->waveK, nMax);

        for(p = 0; p < pMax; p++)
          for(q = 0; q < qMax; q++) {
            T_AB[p][q] = AB.diagonal(q, p);
            T_AB[p + pMax][q + qMax] = AB.diagonal(q, p);
            T_AB[p + pMax][q] = AB.offdiagonal(q, p);
            T_AB[p][q + qMax] = AB.offdiagonal(q, p);
          }

        // In this case, T needs to correspond to the SECOND object; so get it
        // here
        geometry->getTLocal(incWave->omega, j, nMax, T);

        // Multiply T_AB and T to get T_fin = -T_AB*T
        Algebra::multiplyMatrixMatrix(T_AB, 2 * pMax, 2 * pMax, T, 2 * pMax, 2 * pMax, T_fin,
                                      consCm1, consC0);

        // Push T_fin into the S matrix
        pushToMatrix(T_fin, 2 * pMax, 2 * pMax, S, i * 2 * pMax, j * 2 * pMax);
      }
    }
  }

  // Clean up memory

  for(p = 0; p < (int)(2 * p.max(nMax)); p++) {
    delete[] T[p];
    delete[] Ti[p];
    delete[] T_fin[p];
    delete[] T_AB[p];
  }

  delete[] T;
  delete[] Ti;
  delete[] T_fin;
  delete[] T_AB;

  return 0;
}

void Solver::update(Geometry *geometry_, Excitation *incWave_, long nMax_) {
  geometry = geometry_;
  incWave = incWave_;
  nMax = nMax_;

  result_FF = NULL;

  populate();
}

void Solver::solveLinearSystem(optimet::Matrix<optimet::t_complex> const &A,
                               optimet::Vector<optimet::t_complex> const &b,
                               optimet::t_complex *x) const {
  using namespace optimet;
#ifdef OPTIMET_MPI
// scalapack::Sizes const size{static_cast<t_uint>(A.rows()), static_cast<t_uint>(A.cols())};
// scalapack::Sizes const block{parallel_params().block_size, parallel_params().block_size};
// scalapack::Sizes const grid = parallel_params().grid.rows * parallel_params().grid.cols != 0 ?
//   parallel_params().grid: scalapack::squarest_largest_grid(communicator().size());
// scalapack::Matrix<t_complex> Aserial({1, 1}, size, block);
// scalapack::Matrix<t_complex> bserial({1, 1}, size, block);
// Aserial.local() = A;
// bserial.local() = b;
#endif
  Vector<t_complex>::Map(x, A.cols()) = A.colPivHouseholderQr().solve(b);
}
