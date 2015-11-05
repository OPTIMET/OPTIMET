//#include <complex>
//#include "mkl.h"
#include "mkl_lapacke.h"
#include "mkl_scalapack.h"
#include "mkl_blacs.h"

//using std::complex;

extern "C" {
    /* Cblacs declarations */
    void Cblacs_pinfo(int*, int*);
    void Cblacs_get(int, int, int*);
    void Cblacs_gridinit(int*, const char*, int, int);
    void Cblacs_pcoord(int, int, int*, int*);
    void Cblacs_gridexit(int);
    void Cblacs_barrier(int, const char*);
    void Czgerv2d(int, int, int, complex<double>*, int, int, int);
    void Czgesd2d(int, int, int, complex<double>*, int, int, int);
 
	// The numroc function (NUMber of Rows Or Columns) computes how many rows and columns a given process owns
	// Then nrows*ncols will be the number of total local entries (i.e. the entries of A_loc).
    int numroc_(int*, int*, int*, int*, int*);

    void descinit_(int *, int *, int *, int *, int *, int *, int *, int *, int *, int *);
}

