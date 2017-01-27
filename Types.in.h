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

#ifndef OPTIMET_TYPES_H
#define OPTIMET_TYPES_H

#include <complex>
#include <functional>
#include <Eigen/Core>

#cmakedefine OPTIMET_BELOS
#cmakedefine OPTIMET_MPI
#ifdef OPTIMET_MPI
#cmakedefine OPTIMET_SCALAPACK
#endif
#cmakedefine OPTIMET_CHAR_ARCH
#cmakedefine OPTIMET_LONG_ARCH
#cmakedefine OPTIMET_ULONG_ARCH

namespace optimet {
//! Root of the type hierarchy for signed integers
typedef int t_int;
//! Root of the type hierarchy for unsigned integers
typedef std::size_t t_uint;
//! Root of the type hierarchy for real numbers
typedef double t_real;
//! Root of the type hierarchy for (real) complex numbers
typedef std::complex<t_real> t_complex;

//! \brief A vector of a given type
//! \details Operates as mathematical vector.
//! \note Scalapack probably expects column-major. Best not to offend it.
template <class T = t_complex>
using Vector = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;
//! \brief A matrix of a given type
//! \details Operates as mathematical matrix.
//! \note Scalapack probably expects column-major. Best not to offend it.
template <class T = t_complex>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
}
#endif
