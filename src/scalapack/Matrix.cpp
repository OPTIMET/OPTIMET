// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

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
const decltype(&OPTIMET_FC_GLOBAL(pdgemm, PDGEMM)) fortran_pdgemm<double>::pointer =
    &OPTIMET_FC_GLOBAL(pdgemm, PDGEMM);
template <> struct fortran_pdgemm<std::complex<double>> {
  const static decltype(&OPTIMET_FC_GLOBAL(pzgemm, PZGEMM)) pointer;
};
const decltype(&OPTIMET_FC_GLOBAL(pzgemm, PZGEMM)) fortran_pdgemm<std::complex<double>>::pointer =
    &OPTIMET_FC_GLOBAL(pzgemm, PZGEMM);

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
    c_in_a.transfer_to(a.context(), c);
    return;
  }
  int m(c.rows()), n(c.cols()), k(opa != 'N' ? a.rows() : a.cols());
  int zero = 1;
  typedef typename MatrixTraits<SCALAR_A>::Scalar Scalar;
  if(a.context().is_valid())
    (*fortran_pdgemm<Scalar>::pointer)(
        &opa, &opb, &m, &n, &k, &alpha, const_cast<Scalar *>(a.local().data()), &zero, &zero,
        const_cast<int *>(a.blacs().data()), const_cast<Scalar *>(b.local().data()), &zero, &zero,
        const_cast<int *>(b.blacs().data()), &beta, const_cast<Scalar *>(c.local().data()), &zero,
        &zero, const_cast<int *>(c.blacs().data()));
}
}

#define OPTIMET_MACRO(SCALAR, A, B, C)                                                             \
  template <>                                                                                      \
  void pdgemm<SCALAR>(SCALAR alpha, Matrix<SCALAR A> const &a, Matrix<SCALAR B> const &b,          \
                      SCALAR beta, Matrix<SCALAR C> &c, char opa, char opb) {                      \
    pdgemm_(alpha, a, b, beta, c, opa, opb);                                                       \
  }

OPTIMET_MACRO(double, , , )
OPTIMET_MACRO(double, , const *, *)
OPTIMET_MACRO(double, const *, const *, *)
OPTIMET_MACRO(std::complex<double>, , , )
OPTIMET_MACRO(std::complex<double>, , const *, *)
OPTIMET_MACRO(std::complex<double>, const *, const *, *)
#undef OPTIMET_MACRO
} // scalapack
} // optimet
