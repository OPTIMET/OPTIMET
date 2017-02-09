// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#include "Algebra.h"

#include <stddef.h>
#include <gsl/gsl_cblas.h>

Algebra::Algebra() {
  //
}

Algebra::~Algebra() {
  //
}

void Algebra::multiplyMatrixMatrix(std::complex<double> **A, int rows_A_,
                                   int cols_A_, std::complex<double> **B,
                                   int rows_B_, int cols_B_,
                                   std::complex<double> **C,
                                   std::complex<double> alpha_,
                                   std::complex<double> beta_) {

  // GNU Scientific Library CBLas implementation
  std::complex<double> *A_cblas = new std::complex<double>[rows_A_ * cols_A_];
  std::complex<double> *B_cblas = new std::complex<double>[rows_B_ * cols_B_];
  std::complex<double> *C_cblas = new std::complex<double>[rows_A_ * cols_B_];

  matrixToVector(rows_A_, cols_A_, A, A_cblas);
  matrixToVector(rows_B_, cols_B_, B, B_cblas);

  cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows_A_, cols_B_,
              cols_A_, &alpha_, A_cblas, rows_A_, B_cblas, rows_B_, &beta_,
              C_cblas, rows_A_);

  vectorToMatrix(rows_A_, cols_B_, C_cblas, C);

  delete[] A_cblas;
  delete[] B_cblas;
  delete[] C_cblas;
}

void Algebra::multiplyVectorMatrix(std::complex<double> **A, int rows_A_,
                                   int cols_A_, std::complex<double> const *X,
                                   std::complex<double> *Y,
                                   std::complex<double> alpha_,
                                   std::complex<double> beta_) {
  std::complex<double> *A_cblas = new std::complex<double>[rows_A_ * cols_A_];

  matrixToVector(rows_A_, cols_A_, A, A_cblas);

  cblas_zgemv(CblasRowMajor, CblasNoTrans, rows_A_, cols_A_, &alpha_, A_cblas,
              rows_A_, X, 1, &beta_, Y, 1);

  delete[] A_cblas;
}

void Algebra::matrixToVector(long rows_, long columns_,
                             std::complex<double> **T_,
                             std::complex<double> *V_) {
  int i, j;

  for (i = 0; i < rows_; i++)
    for (j = 0; j < columns_; j++) {
      V_[i * rows_ + j] = T_[i][j];
    }
}

void Algebra::vectorToMatrix(long rows_, long columns_,
                             std::complex<double> *V_,
                             std::complex<double> **T_) {
  int i, j;
  for (i = 0; i < rows_; i++)
    for (j = 0; j < columns_; j++) {
      T_[i][j] = V_[i * rows_ + j];
    }
}
