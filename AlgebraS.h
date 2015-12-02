#ifndef ALGEBRAS_H_
#define ALGEBRAS_H_

#include <complex>

/**
 * The AlgebraS class provides wrappers to common Linear Algebra operations.
 * AlgebraS is used to ensure portability between different LAPACK, SCALAPACK and
 * BLAS implementations. AlgebraS separates from Algebra to ensure GSL compatibility.
 */
class AlgebraS
{
public:
  /**
   * Default constructor for the Algebra_ class.
   */
  AlgebraS();

  /**
   * Default destructor for the Algebra class.
   */
  virtual ~AlgebraS();

  /**
   * Solves the A*x = b equation.
   * @param A the A matrix.
   * @param rows_A_ the number of rows of A.
   * @param cols_A_ the number of columns of A.
   * @param b the b vector (free terms).
   * @param x the x vector (solution return).
   */
  static void solveMatrixVector(std::complex<double> **A, int rows_A_, int cols_A_, std::complex<double> *b, std::complex<double> *x);
};

#endif /* ALGEBRA__H_ */
