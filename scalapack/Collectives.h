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

#ifndef OPTIMET_SCALAPACK_COLLECTIVES_H
#define OPTIMET_SCALAPACK_COLLECTIVES_H
#include "Types.h"
#ifdef OPTIMET_SCALAPACK

#include <mpi.h>
#include <type_traits>
#include <vector>

namespace optimet {
namespace scalapack {
namespace details {
//! True if a type scalapack knows about
template <class T>
struct is_fundamental
    : std::integral_constant<bool, std::is_same<T, int>::value or std::is_same<T, float>::value or
                                       std::is_same<T, double>::value or
                                       std::is_same<T, std::complex<float>>::value or
                                       std::is_same<T, std::complex<double>>::value> {};
}

class Context;
#define OPTIMET_MACRO(TYPE)                                                                        \
  /** Broadcasts from given process to all others **/                                              \
  TYPE broadcast(TYPE const &, Context const &context, t_uint row, t_uint col);                    \
  Matrix<TYPE> broadcast(Matrix<TYPE> const &, Context const &context, t_uint row, t_uint col);    \
  Vector<TYPE> broadcast(Vector<TYPE> const &, Context const &context, t_uint row, t_uint col);
OPTIMET_MACRO(int);
OPTIMET_MACRO(float);
OPTIMET_MACRO(double);
OPTIMET_MACRO(std::complex<double>);
OPTIMET_MACRO(std::complex<float>);
#undef OPTIMET_MACRO
} /* optime::mpi */
} /* optimet */
#endif
#endif
