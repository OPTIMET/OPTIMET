#include "AlgebraS.h"
#include "types.h"
#include <Eigen/Dense>

namespace optimet {
namespace algebra {

void solveMatrixVector(Matrix<t_complex> const &A, Vector<t_complex> const &b,
                       t_complex *x) {
  // Solve A*x = b
  Vector<t_complex>::Map(x, A.cols()) = A.colPivHouseholderQr().solve(b);
}

} // namespace algebra
} // namespace optimet
