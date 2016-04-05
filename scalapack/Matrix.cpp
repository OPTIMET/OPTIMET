#include "scalapack/Blacs.h"
#include "scalapack/Matrix.h"
#include <iostream>

namespace optimet {
namespace scalapack {
namespace {
// Helpers so we can call different functions at the same place in the code
// note: c++14's template variable would make this simpler.
template <class SCALAR> struct fortran_pdgemm;
template <> struct fortran_pdgemm<double> {
  const static decltype(&OPTIMET_FC_GLOBAL(pdgemm, PDGEMM)) pointer;
};
const auto fortran_pdgemm<double>::pointer = &OPTIMET_FC_GLOBAL(pdgemm, PDGEMM);
template <> struct fortran_pdgemm<std::complex<double>> {
  const static decltype(&OPTIMET_FC_GLOBAL(pzgemm, PZGEMM)) pointer;
};
const auto fortran_pdgemm<std::complex<double>>::pointer = &OPTIMET_FC_GLOBAL(pzgemm, PZGEMM);

template <class SCALAR_A, class SCALAR_B, class SCALAR_C>
void pdgemm_(typename MatrixTraits<SCALAR_A>::Scalar alpha, Matrix<SCALAR_A> const &a,
             Matrix<SCALAR_B> const &b, typename MatrixTraits<SCALAR_A>::Scalar beta,
             Matrix<SCALAR_C> &c, char opa = 'N', char opb = 'N') {
  if((opa != 'N' ? a.rows() : a.cols()) != (opb != 'N' ? b.cols() : b.rows()))
    throw std::runtime_error("Incompatible a and b matrices");
  if((opa != 'N' ? a.cols() : a.rows()) != c.rows())
    throw std::runtime_error("Incompatible a and c matrices");
  if((opb != 'N' ? b.rows() : b.cols()) != c.cols())
    throw std::runtime_error("Incompatible b and c matrices");
  if(a.context() != b.context()) {
    auto b_in_a = b.transfer_to(a.context());
    pdgemm_(alpha, a, b_in_a, beta, c, opa, opb);
    return;
  }
  if(a.context() != c.context()) {
    auto c_in_a = c.transfer_to(a.context());
    pdgemm_(alpha, a, b, beta, c_in_a, opa, opb);
    c_in_a.transfer_to(c);
    return;
  }
  int m(c.rows()), n(c.cols()), k(opa != 'N' ? a.rows() : a.cols());
  int zero = 1;
  typedef typename MatrixTraits<SCALAR_A>::Scalar Scalar;
  (*fortran_pdgemm<Scalar>::pointer)(
      &opa, &opb, &m, &n, &k, &alpha, const_cast<Scalar *>(a.local().data()), &zero, &zero,
      const_cast<int *>(a.blacs().data()), const_cast<Scalar *>(b.local().data()), &zero, &zero,
      const_cast<int *>(b.blacs().data()), &beta, const_cast<Scalar *>(c.local().data()), &zero,
      &zero, const_cast<int *>(c.blacs().data()));
}
}

#define OPTIMET_MACRO(TYPE, CONST, POINTER)                                                        \
  template <>                                                                                      \
  void pdgemm<TYPE>(TYPE alpha, Matrix<TYPE> const &a, Matrix<TYPE CONST POINTER> const &b,        \
                    TYPE beta, Matrix<TYPE POINTER> &c, char opa, char opb) {                      \
    pdgemm_(alpha, a, b, beta, c, opa, opb);                                                       \
  }

OPTIMET_MACRO(double, , )
OPTIMET_MACRO(double, const, *)
OPTIMET_MACRO(std::complex<double>, , )
OPTIMET_MACRO(std::complex<double>, const, *)
#undef OPTIMET_MACRO
} // scalapack
} // optimet
