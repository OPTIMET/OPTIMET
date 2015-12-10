#include "AlgebraS.h"
#include "types.h"
#include <Eigen/Dense>

namespace optimet { namespace algebra {

void solveMatrixVector(t_complex **A, t_uint rows_A_, t_uint cols_A_, t_complex *b, t_complex *x) {
  //Convert A explicitly, since we don't know whether it is contiguous
  Matrix<> A_(rows_A_, cols_A_);
  for(t_uint i=0; i<rows_A_; i++)
    for(t_uint j=0; j<cols_A_; j++)
      A_(i,j) = A[i][j];

  //Solve A*x = b
  Vector<t_complex>::Map(x, cols_A_)
    = A_.colPivHouseholderQr().solve(Vector<t_complex>::Map(b, cols_A_));
}
}} // namespace optimet::algebra
