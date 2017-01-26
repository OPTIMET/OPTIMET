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

#include <cassert>
#include "scalapack/Collectives.h"
#include "scalapack/Context.h"
#include "scalapack/Blacs.h"

namespace optimet {
namespace scalapack {

namespace {

#define OPTIMET_MACRO(letter, LETTER, TYPE)                                                                                                     \
  void gebs2d(int context, TYPE const &value) {                             \
    char const *scope = "All";                                              \
    char const *topology = " ";                                             \
    int one = 1;                                                            \
    OPTIMET_FC_GLOBAL(letter ## gebs2d, LETTER ## gebs2d)(                  \
        &context, scope, topology, &one, &one, &value, &one);               \
  }                                                                         \
  void gebr2d(int context, TYPE &value, int row, int col) {                 \
    char const *scope = "All";                                              \
    char const *topology = " ";                                             \
    int one = 1;                                                            \
    OPTIMET_FC_GLOBAL(letter ## gebr2d, LETTER ## gebr2d)(                  \
        &context, scope, topology, &one, &one, &value, &one, &row, &col);   \
  }                                                                         \
  void gebs2d(int context, int m, int n, TYPE const *value) {               \
    char const *scope = "All";                                              \
    char const *topology = " ";                                             \
    OPTIMET_FC_GLOBAL(letter ## gebs2d, LETTER ## gebs2d)(                  \
        &context, scope, topology, &m, &n, value, &m);                      \
  }                                                                         \
  void gebr2d(int context, int m, int n, TYPE *value, int row, int col) {   \
    char const *scope = "All";                                              \
    char const *topology = " ";                                             \
    OPTIMET_FC_GLOBAL(letter ## gebr2d, LETTER ## gebr2d)(                  \
        &context, scope, topology, &m, &n, value, &m, &row, &col);          \
  }
OPTIMET_MACRO(i, I, int);
OPTIMET_MACRO(s, S, float);
OPTIMET_MACRO(d, D, double);
OPTIMET_MACRO(c, C, std::complex<float>);
OPTIMET_MACRO(z, Z, std::complex<double>);
#undef OPTIMET_MACRO

template <class T>
typename std::enable_if<details::is_fundamental<T>::value, T>::type
tm_broadcast(T const &value, Context const &context, t_uint row, t_uint col) {
  if(not context.is_valid())
    return value;
  assert(context.rows() > row);
  assert(context.cols() > col);
  auto const is_root =
      static_cast<t_uint>(context.row()) == row and static_cast<t_uint>(context.col()) == col;
  T result(value);
  if(is_root)
    gebs2d(*context, value);
  else
    gebr2d(*context, result, row, col);
  return result;
}

template <class T>
typename std::enable_if<details::is_fundamental<T>::value, Matrix<T>>::type
tm_broadcast(Matrix<T> const &value, Context const &context, t_uint row, t_uint col) {
  if(not context.is_valid())
    return value;
  assert(context.rows() > row);
  assert(context.cols() > col);
  auto const is_root =
      static_cast<t_uint>(context.row()) == row and static_cast<t_uint>(context.col()) == col;
  if(is_root) {
    gebs2d(*context, static_cast<int>(value.rows()));
    gebs2d(*context, static_cast<int>(value.cols()));
    gebs2d(*context, value.rows(), value.cols(), value.data());
    return value;
  } else {
    int matrix_rows, matrix_cols;
    gebr2d(*context, matrix_rows, row, col);
    gebr2d(*context, matrix_cols, row, col);
    Matrix<T> result(matrix_rows, matrix_cols);
    gebr2d(*context, matrix_rows, matrix_cols, result.data(), row, col);
    return result;
  }
}

template <class T>
typename std::enable_if<details::is_fundamental<T>::value, Vector<T>>::type
tm_broadcast(Vector<T> const &value, Context const &context, t_uint row, t_uint col) {
  if(not context.is_valid())
    return value;
  assert(context.rows() > row);
  assert(context.cols() > col);
  auto const is_root =
      static_cast<t_uint>(context.row()) == row and static_cast<t_uint>(context.col()) == col;
  if(is_root) {
    gebs2d(*context, static_cast<int>(value.size()));
    gebs2d(*context, value.size(), 1, value.data());
    return value;
  } else {
    int vector_size;
    gebr2d(*context, vector_size, row, col);
    Vector<T> result(vector_size);
    gebr2d(*context, vector_size, 1, result.data(), row, col);
    return result;
  }
}
} // anonymous namespace

#define OPTIMET_MACRO(TYPE)                                                                        \
  TYPE broadcast(TYPE const &value, Context const &context, t_uint row, t_uint col) {              \
    return tm_broadcast(value, context, row, col);                                                 \
  }                                                                                                \
  Matrix<TYPE> broadcast(Matrix<TYPE> const &value, Context const &context,                        \
      t_uint row, t_uint col) {                                                                    \
    return tm_broadcast(value, context, row, col);                                                 \
  }                                                                                                \
  Vector<TYPE> broadcast(Vector<TYPE> const &value, Context const &context,                        \
      t_uint row, t_uint col) {                                                                    \
    return tm_broadcast(value, context, row, col);                                                 \
  }
OPTIMET_MACRO(int);
OPTIMET_MACRO(float);
OPTIMET_MACRO(double);
OPTIMET_MACRO(std::complex<float>);
OPTIMET_MACRO(std::complex<double>);
#undef OPTIMET_MACRO

} /* scalapack  */
} /* optimet  */
