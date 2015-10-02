#include "Algebra.h"

Algebra::Algebra()
{
	//

}

Algebra::~Algebra()
{
	//
}

void Algebra::multiplyMatrixMatrix(complex<double>** A, int rows_A_, int cols_A_,
		complex<double>** B, int rows_B_, int cols_B_,
		complex<double>** C, complex<double> alpha_, complex<double> beta_)
{

	//GNU Scientific Library CBLas implementation
	complex<double> *A_cblas = new complex<double>[rows_A_ * cols_A_];
	complex<double> *B_cblas = new complex<double>[rows_B_ * cols_B_];
	complex<double> *C_cblas = new complex<double>[rows_A_ * cols_B_];

	matrixToVector(rows_A_, cols_A_, A, A_cblas);
	matrixToVector(rows_B_, cols_B_, B, B_cblas);

	cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, rows_A_, cols_B_, cols_A_, &alpha_, A_cblas, rows_A_, B_cblas, rows_B_, &beta_, C_cblas, rows_A_);

	vectorToMatrix(rows_A_, cols_B_, C_cblas, C);

	delete [] A_cblas;
	delete [] B_cblas;
	delete [] C_cblas;
}

void Algebra::multiplyVectorMatrix(complex<double>** A, int rows_A_,
		int cols_A_, complex<double>* X, complex<double>* Y,
		complex<double> alpha_, complex<double> beta_)
{
	complex<double> *A_cblas = new complex<double>[rows_A_ * cols_A_];

	matrixToVector(rows_A_, cols_A_, A, A_cblas);

	cblas_zgemv(CblasRowMajor, CblasNoTrans, rows_A_, cols_A_, &alpha_, A_cblas, rows_A_, X, 1, &beta_, Y, 1);

	delete [] A_cblas;
}

void Algebra::matrixToVector(long rows_, long columns_, complex<double>** T_,
		complex<double>* V_)
{
	int i,j;

	for(i=0; i<rows_; i++)
		for(j=0; j<columns_; j++)
		{
			V_[i*rows_ + j] = T_[i][j];
		}
}

void Algebra::vectorToMatrix(long rows_, long columns_, complex<double>* V_,
		complex<double>** T_)
{
	int i,j;
	for(i=0; i<rows_; i++)
		for(j=0; j<columns_; j++)
		{
			T_[i][j] = V_[i*rows_ + j];
		}
}

