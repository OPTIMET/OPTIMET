#ifndef ALGEBRAS_H_
#define ALGEBRAS_H_
#include "Tools.h"
//#include <armadillo>
#include <complex>
#define MKL_Complex16 std::complex<double>
#include "externals.h"

using std::complex;

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
	static void solveMatrixVector		(complex<double> **A, int rows_A_, int cols_A_, complex<double> *b, complex<double> *x);

	static void solveMatrixVectorBP_AJ	(complex<double> **A, int rows_A_, int cols_A_, complex<double> *b, complex<double> *x, int pMax_, int noRowBlocks_, int noColBlocks_);

	static void solveMatrixVectorBP_CB(complex<double> **A, int rows_A_, int cols_A_, complex<double> *b, complex<double> *x, int nMax_, int noColBlocks_, int context_);
};

#endif /* ALGEBRA__H_ */
