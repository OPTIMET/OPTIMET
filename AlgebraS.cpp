#include "AlgebraS.h"

//#include "mkl.h"
//#include "mkl_lapacke.h"
using namespace std;

//using arma::solve;
//using arma::Mat;
//using arma::Col;
//using arma::cx_double;
//#include "mpi.h"

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
/*
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
*/

	int matrix_order = LAPACK_ROW_MAJOR;
	char fact = 'E';
	char trans = 'N';
	int n = rows_A_;
	int nrhs = 1;
	int lda = cols_A_;
	complex<double> *af = new complex<double>[rows_A_ * cols_A_];
	int ldaf = cols_A_;
	int *ipiv = new int[rows_A_ * cols_A_];
	char equed;
	double *r = new double[rows_A_];
	double *c = new double[rows_A_];
	int ldb = 1;
	int ldx = 1;
	double rcond;
	double ferr;
	double berr;
	double rpivot;

//	LAPACKE_zgesvx(matrix_order, fact, trans, n, nrhs, &A[0][0], lda, af, ldaf, ipiv, &equed, r, c, b, ldb, x, ldx, &rcond, &ferr, &berr, &rpivot);

	// AJ --------------------------------------------------------
	//Convert A from 2D to 1D matrix
	int N 		= rows_A_;				// Total number of unknowns
	int NN 		= rows_A_*cols_A_;		// Total number of unknowns squared
	int NRHS	= 1;					// number of RHS to be solved
	int LDA 	= cols_A_;				// leading dimension in square matrix - array A
	int LDB 	= 1;					// leading dimension in vector array  - array B

	int *IPIV = new int[N];
	complex<double> *oneD_A;
	oneD_A = new complex<double>[NN];

	int ii =0;
//	cout<<" ALgebraS op SS matrix values ------------------------------------------------------------------ "<<endl;
	for(int i=0; i<rows_A_; i++){
		for(int j=0; j<cols_A_; j++){
			oneD_A[ii] = A[i][j];
//			cout<<A[i][j]<<" ";
			ii++;
		}
	//	cout<<endl;
	}
//	cout<<endl;

//	cout<<" ALgebraS op B matrix values ------------------------------------------------------------------ "<<endl;
//	for(ii=0; ii<rows_A_; ii++){
//		cout<<b[ii]<<" ";
//	}
//	cout<<endl;

	LAPACKE_zgesv(matrix_order, N, NRHS, oneD_A, LDA, IPIV, b, LDB);
	for(ii=0; ii<rows_A_; ii++){
		x[ii]=b[ii];
	}

//	cout<<" ALgebraS op E matrix values ------------------------------------------------------------------ "<<endl;
//	for(ii=0; ii<rows_A_; ii++){
//		cout<<b[ii]<<" ";
//	}
//	cout<<endl;

	delete [] IPIV;
	// AJ --------------------------------------------------------

	delete [] af;
	delete [] ipiv;
	delete [] r;
	delete [] c;
}

void AlgebraS::solveMatrixVectorBP_AJ(complex<double> **A, int rows_A_, int cols_A_, complex<double> *b, complex<double> *x, int pMax_, int noRowBlocks_, int noColBlocks_)
{

}

// Parallel implemntation of solveMatrixVectorBP
void AlgebraS::solveMatrixVectorBP_CB(complex<double> **A, int rows_A_, int cols_A_, complex<double> *b, complex<double> *x, int nMax_, int noColBlocks_, int context_)
{
    //Block first to ensure all processes are lined up
    Cblacs_barrier(context_, "All");

    int myrow, mycol, myid, numproc;
    Cblacs_pinfo(&myid, &numproc);
    Cblacs_pcoord(context_, myid, &myrow, &mycol);

    int desca[9], descb[9];

    int block_rows_A = rows_A_;
    int block_cols_A = cols_A_;

    //if(mycol == numproc-1)
    //{
        //block_cols_A = nMax_ * (noColBlocks_ / numproc + noColBlocks_ % numproc);
    //}
    //else
    //{
        //block_cols_A = nMax_ * (noColBlocks_ / numproc);
    //}

	int ZERO = 0;
    int ONE = 1;

    int IERR;

    int *ipiv = new int[block_rows_A * block_cols_A];
    
    descinit_(desca, &rows_A_, &cols_A_, &block_rows_A, &block_cols_A, &ZERO, &ZERO, &context_, &block_cols_A, &IERR);
    descinit_(descb, &rows_A_, &ONE, &block_rows_A, &ONE, &ZERO, &ZERO, &context_, &block_cols_A, &IERR);

    int aux_rows = myrow * block_rows_A;
    int aux_cols = mycol * block_cols_A;

    pzgesv(&block_rows_A, &ONE, &A[0][0], &ONE, &ONE, desca, ipiv, b, &ONE, &ONE, descb, &IERR);
    
    //exit(0);
 
    //Block to make sure all processes are here
    Cblacs_barrier(context_, "All");

    //Print the results
    //for(int id=0; id<numproc; id++)
    //{
        //if(myid == id)
        //{
            //std::cout << "b on node " << myid <<"\n";
            //for(int i=0; i<block_rows_A; i++)
                //std::cout << b[i] << " ";
            //std::cout << "\n";
        //}
        //Cblacs_barrier(context_, "All");
    //}

    //Gather the data in the first process
    
    int sendr = 0;

    for(sendr = 0; sendr < (int)sqrt(numproc); sendr++)
    {
            if(myrow == sendr)
            {
                Czgesd2d(context_, 1, block_cols_A, b, 1, 0, 0);
            }

            if(myid == 0)
            {
                Czgerv2d(context_, 1, block_cols_A, &x[sendr * block_cols_A], 1, sendr, 0);
            }
    }
}