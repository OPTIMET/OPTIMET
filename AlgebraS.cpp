#include "AlgebraS.h"

#include <armadillo>

using arma::solve;
using arma::Mat;
using arma::Col;
using arma::cx_double;

AlgebraS::AlgebraS()
{
	//

}

AlgebraS::~AlgebraS()
{
	//
}

void AlgebraS::solveMatrixVector(complex<double> **A, int rows_A_, int cols_A_, complex<double> *b, complex<double> *x)
{
	//Armadillo implementation (fast LAPACK++ prototyping tool)
	Mat<cx_double> A_arma;
	Col<cx_double> x_arma;
	Col<cx_double> b_arma;

	//Convert A to A_arma
	A_arma.set_size(rows_A_, cols_A_);
	for(int i=0; i<rows_A_; i++)
		for(int j=0; j<cols_A_; j++)
			A_arma(i,j) = A[i][j];

	//Convert b to b_arma
	b_arma.set_size(rows_A_);
	for(int i=0; i<rows_A_; i++)
		b_arma(i) = b[i];

	//Solve the A*x = b
	//x_arma = solve(A_arma, b_arma, true);
	Mat<cx_double> inv_A = arma::inv(A_arma);

	x_arma = inv_A * b_arma;

	//Convert x_arma to x
	for(int i=0; i<rows_A_; i++)
		x[i] = x_arma(i);
}
