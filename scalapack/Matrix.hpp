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
#include "scalapack/InitExit.h"
#include "scalapack/Matrix.h"
#include <iostream>

namespace optimet {
namespace scalapack {

#define OPTIMET_MACRO(NAME)                                                                        \
  template <class SCALAR>                                                                          \
  typename Matrix<SCALAR>::EigenMatrix::Index Matrix<SCALAR>::local_##NAME##s(                     \
      Context const &context, Sizes size, Sizes blocks, Index index) {                             \
    if(not context.is_valid())                                                                     \
      return 0;                                                                                    \
    int n = static_cast<int>(size.NAME##s);                                                        \
    int nb = static_cast<int>(blocks.NAME##s);                                                     \
    int iproc = static_cast<int>(context.NAME());                                                  \
    int isrc = static_cast<int>(index.NAME);                                                       \
    int nprocs = static_cast<int>(context.NAME##s());                                              \
    auto const result = OPTIMET_FC_GLOBAL(numroc, NUMROC)(&n, &nb, &iproc, &isrc, &nprocs);        \
    return static_cast<typename Matrix<SCALAR>::EigenMatrix::Index>(result);                       \
  }
OPTIMET_MACRO(row);
OPTIMET_MACRO(col);
#undef OPTIMET_MACRO

namespace {
#define OPTIMET_MACRO(LETTER, TYPE)                                                                \
  inline void gemr2d(int m, int n, TYPE *ptrmyblock, int ia, int ja, int *ma, TYPE *ptrmynewblock, \
                     int ib, int jb, int *mb, int globcontext) {                                   \
    Cp##LETTER##gemr2d(m, n, ptrmyblock, ia, ja, ma, ptrmynewblock, ib, jb, mb, globcontext);      \
  }
OPTIMET_MACRO(s, float);
OPTIMET_MACRO(d, double);
OPTIMET_MACRO(c, std::complex<float>);
OPTIMET_MACRO(z, std::complex<double>);
#undef OPTIMET_MACRO
}

template <class SCALAR>
void Matrix<SCALAR>::transfer_to(Context const &un, MapMatrix &other) const {
  if(context().is_valid() or un.is_valid() or other.context().is_valid()) {
    Scalar dummy;
    Scalar const *input = local().size() > 0 ? local().data() : &dummy;
    Scalar *output = other.local().size() > 0 ? other.local().data() : &dummy;
    gemr2d(rows(), cols(), const_cast<Scalar *>(input), 1, 1, const_cast<int *>(blacs().data()),
           output, 1, 1, const_cast<int *>(other.blacs().data()), *un);
  }
}
template <class SCALAR>
void Matrix<SCALAR>::transfer_to(Context const &un, ConcreteMatrix &other) const {
  if(context().is_valid() or un.is_valid() or other.context().is_valid()) {
    Scalar dummy;
    Scalar const *input = local().size() > 0 ? local().data() : &dummy;
    Scalar *output = other.local().size() > 0 ? other.local().data() : &dummy;
    gemr2d(rows(), cols(), const_cast<Scalar *>(input), 1, 1, const_cast<int *>(blacs().data()),
           output, 1, 1, const_cast<int *>(other.blacs().data()), *un);
  }
}
template <class SCALAR>
typename Matrix<SCALAR>::ConcreteMatrix
Matrix<SCALAR>::transfer_to(Context const &un, Context const &other, Sizes const &_blocks,
                            Index const &_index) const {
  ConcreteMatrix result(
      other, sizes(),
      {_blocks.rows == std::numeric_limits<t_uint>::max() ? blocks().rows : _blocks.rows,
       _blocks.cols == std::numeric_limits<t_uint>::max() ? blocks().cols : _blocks.cols},
      {_index.row == std::numeric_limits<t_uint>::max() ? index().row : _index.row,
       _index.col == std::numeric_limits<t_uint>::max() ? index().col : _index.col});
  transfer_to(un, result);
  return result;
}

template <class SCALAR>
std::tuple<t_uint, t_uint, t_uint, t_uint>
Matrix<SCALAR>::local_indices(std::tuple<t_uint, t_uint> const &i) const {
  int np_row(context().rows()), np_col(context().cols());
  int nb_row(blocks().rows), nb_col(blocks().cols);
  int i_row(std::get<0>(i)), i_col(std::get<1>(i));
  int f_row(first_row()), f_col(first_col());
  int dummy = 0;
  return std::tuple<t_uint, t_uint, t_uint, t_uint>(
      OPTIMET_FC_GLOBAL(indxg2l, INDXG2L)(&i_row, &nb_row, &dummy, &f_row, &np_row),
      OPTIMET_FC_GLOBAL(indxg2l, INDXG2L)(&i_col, &nb_col, &dummy, &f_col, &np_col),
      OPTIMET_FC_GLOBAL(indxg2p, INDXG2P)(&i_row, &nb_row, &dummy, &f_row, &np_row),
      OPTIMET_FC_GLOBAL(indxg2p, INDXG2P)(&i_col, &nb_col, &dummy, &f_col, &np_col));
}

template <class SCALAR> void Matrix<SCALAR>::operator=(Matrix<SCALAR> const &other) {
  if(rows() != other.rows() or cols() != other.cols())
    throw std::runtime_error("Matrices have different sizes.");
  if(context() != other.context())
    *this = other.transfer_to(context());
  else
    local() = other.local();
}

template <class SCALAR> void Matrix<SCALAR>::operator=(Matrix<SCALAR> &&other) {
  if(rows() != other.rows() or cols() != other.cols())
    throw std::runtime_error("Matrices have different sizes.");
  if(context() != other.context())
    *this = other.transfer_to(context());
  else
    local().swap(other.local());
}

// Explicit declaration of full specialization
template <>
void pdgemm<double>(double alpha, Matrix<double> const &a, Matrix<double> const &b, double beta,
                    Matrix<double> &c, char opa, char opb);
template <>
void pdgemm<double>(double alpha, Matrix<double const *> const &a, Matrix<double const *> const &b,
                    double beta, Matrix<double *> &c, char opa, char opb);
template <>
void pdgemm<double>(double alpha, Matrix<double> const &a, Matrix<double const *> const &b,
                    double beta, Matrix<double *> &c, char opa, char opb);
template <>
void pdgemm<std::complex<double>>(std::complex<double> alpha, Matrix<std::complex<double>> const &a,
                                  Matrix<std::complex<double>> const &b, std::complex<double> beta,
                                  Matrix<std::complex<double>> &c, char opa, char opb);
template <>
void pdgemm<std::complex<double>>(std::complex<double> alpha, Matrix<std::complex<double>> const &a,
                                  Matrix<std::complex<double> const *> const &b,
                                  std::complex<double> beta, Matrix<std::complex<double> *> &c,
                                  char opa, char opb);
template <>
void pdgemm<std::complex<double>>(std::complex<double> alpha,
                                  Matrix<std::complex<double> const *> const &a,
                                  Matrix<std::complex<double> const *> const &b,
                                  std::complex<double> beta, Matrix<std::complex<double> *> &c,
                                  char opa, char opb);

template <class SCALAR>
Matrix<SCALAR *> map_matrix(optimet::Matrix<SCALAR> &matrix, Context const &context,
                            Sizes const &sizes, Sizes const &blocks) {
  Eigen::Map<optimet::Matrix<SCALAR>> const map(matrix.data(), matrix.rows(), matrix.cols());
  return Matrix<SCALAR *>(map, context, sizes, blocks);
}
template <class SCALAR>
Matrix<SCALAR const *> map_matrix(optimet::Matrix<SCALAR> const &matrix, Context const &context,
                                  Sizes const &sizes, Sizes const &blocks) {
  Eigen::Map<optimet::Matrix<SCALAR> const> const map(matrix.data(), matrix.rows(), matrix.cols());
  return Matrix<SCALAR const*>(map, context, sizes, blocks);
}
template <class SCALAR>
Matrix<SCALAR *> map_matrix(optimet::Vector<SCALAR> &matrix, Context const &context,
                            Sizes const &sizes, Sizes const &blocks) {
  auto const nrows = Matrix<SCALAR>::local_rows(context, sizes, blocks, {0, 0});
  auto const ncols = Matrix<SCALAR>::local_cols(context, sizes, blocks, {0, 0});
  if(nrows * ncols != matrix.size())
    throw std::runtime_error("Incorrect local size");
  Eigen::Map<optimet::Matrix<SCALAR>> const map(matrix.data(), nrows, ncols);
  return Matrix<SCALAR *>(map, context, sizes, blocks);
}
template <class SCALAR>
Matrix<SCALAR const *> map_matrix(optimet::Vector<SCALAR> const &matrix, Context const &context,
                                  Sizes const &sizes, Sizes const &blocks) {
  auto const nrows = Matrix<SCALAR>::local_rows(context, sizes, blocks, {0, 0});
  auto const ncols = Matrix<SCALAR>::local_cols(context, sizes, blocks, {0, 0});
  if(nrows * ncols != matrix.size())
    throw std::runtime_error("Incorrect local size");
  Eigen::Map<optimet::Matrix<SCALAR> const> const map(matrix.data(), nrows, ncols);
  return Matrix<SCALAR const *>(map, context, sizes, blocks);
}
} // scalapack
} // optimet
