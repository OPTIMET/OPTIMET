#ifndef ALGEBRAS_H_
#define ALGEBRAS_H_

#include "Types.h"

namespace optimet { namespace algebra {
  //! Adapter for eigen solve
  void solveMatrixVector(Matrix<t_complex> const &A, Vector<t_complex> const &b, t_complex* x);
}}
#endif /* ALGEBRA__H_ */
