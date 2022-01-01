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

#ifndef ALGEBRA_H_
#define ALGEBRA_H_

#include <complex>

/**
 * The Algebra class provides wrappers to common Linear Algebra operations.
 * Algebra is used to ensure portability between different LAPACK, SCALAPACK and
 * BLAS implementations. Algebra is released under the GSL. To view a copy
 * of the licence, look in the documentation.
 */
class Algebra {
public:
  /**
   * Default constructor for the Algebra class.
   */
  Algebra();

  /**
   * Default destructor for the Algebra class.
   */
  virtual ~Algebra();

  /**
   * Implements the multiplication C = alpha*A*B+beta*C.
   * @param A the A matrix.
   * @param rows_A_ the number of rows of A.
   * @param cols_A_ the number of columns of A.
   * @param B the B matrix.
   * @param rows_B_ the number of rows of B.
   * @param cols_B_ the number of columns of B.
   * @param C the C matrix (return matrix).
   * @param alpha the scalar alpha.
   * @param beta the scalar beta.
   */
  static void multiplyMatrixMatrix(std::complex<double> **A, int rows_A_,
                                   int cols_A_, std::complex<double> **B,
                                   int rows_B_, int cols_B_,
                                   std::complex<double> **C,
                                   std::complex<double> alpha_,
                                   std::complex<double> beta_);

  /**
   * Implements the multiplication Y = alpha*A*X + beta*Y
   * @param A the A matrix.
   * @param rows_A_ the number of rows of A.
   * @param cols_A_ the number of columns of A.
   * @param X the X vector.
   * @param Y the Y vector (return vector).
   * @param alpha_ the scalar alpha.
   * @param beta_ the scalar beta.
   */
  static void multiplyVectorMatrix(std::complex<double> **A, int rows_A_,
                                   int cols_A_, std::complex<double> const *X,
                                   std::complex<double> *Y,
                                   std::complex<double> alpha_,
                                   std::complex<double> beta_);

  /**
   * Convert a Matrix into a Vector.
   * @param rows_ the number of rows of the matrix.
   * @param columns_ the number of columns of the matrix.
   * @param T_ the matrix.
   * @param V_ the vector.
   */
  static void matrixToVector(long rows_, long columns_,
                             std::complex<double> **T_,
                             std::complex<double> *V_);

  /**
   * Convert a Vector into a Matrix.
   * @param rows_ the number of rows of the matrix.
   * @param columns_ the number of columns of the matrix.
   * @param V_ the vector.
   * @param T_ the matrix.
   */
  static void vectorToMatrix(long rows_, long columns_,
                             std::complex<double> *V_,
                             std::complex<double> **T_);
};

#endif /* ALGEBRA_H_ */
