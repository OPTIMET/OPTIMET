#ifndef ALGEBRAS_H_
#define ALGEBRAS_H_

#include "types.h"

namespace optimet { namespace algebra {
  //! Adapter for eigen solve
  void solveMatrixVector(
      t_complex **A, t_uint rows_A_, t_uint cols_A_, t_complex *b, t_complex *x);
}}
#endif /* ALGEBRA__H_ */
