#include <iostream>
#include "scalapack/Matrix.h"
#include "scalapack/Blacs.h"

namespace optimet {
namespace scalapack {


template <>
void pdgemm<double>(double alpha, Matrix<double> const &a, Matrix<double> const &b, double beta,
                    Matrix<double> &c, char opa, char opb) {
  if((opa != 'N' ? a.rows() : a.cols()) != (opb != 'N' ? b.cols() : b.rows()))
    throw std::runtime_error("Incompatible a and b matrices");
  if((opa != 'N' ? a.cols() : a.rows()) != c.rows())
    throw std::runtime_error("Incompatible a and c matrices");
  if((opb != 'N' ? b.rows() : b.cols()) != c.cols())
    throw std::runtime_error("Incompatible b and c matrices");
  if(a.context() != b.context()) {
    auto b_in_a = b.transfer_to(a.context());
    pdgemm<double>(alpha, a, b_in_a, beta, c, opa, opb);
    return;
  }
  if(a.context() != c.context()) {
    auto c_in_a = c.transfer_to(a.context());
    pdgemm<double>(alpha, a, b, beta, c_in_a, opa, opb);
    c = c_in_a;
    return;
  }
  int m(c.rows()), n(c.cols()), k(opa != 'N' ? a.rows() : a.cols());
  int zero = 1;
  return OPTIMET_FC_GLOBAL(pdgemm, pdgemm)(
      &opa, &opb, &m, &n, &k, &alpha, const_cast<double *>(a.local().data()), &zero, &zero,
      const_cast<int *>(a.blacs().data()), const_cast<double *>(b.local().data()), &zero, &zero,
      const_cast<int *>(b.blacs().data()), &beta, const_cast<double *>(c.local().data()), &zero,
      &zero, const_cast<int *>(c.blacs().data())
  );
}

} // scalapack
} // optimet
