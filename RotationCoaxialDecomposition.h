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

#ifndef OPTIMET_ROTATION_COAXIAL_DECOMPOSITION_H
#define OPTIMET_ROTATION_COAXIAL_DECOMPOSITION_H

#include "Coefficients.h"
#include "Types.h"
#include "constants.h"
#include <iostream>

#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace optimet {
//! \brief Rotation-coaxial decomposition of the vector potentials
//! \param[in] wavenumber
//! \param[in] tz: translation alongst the z axis
//! \param[in] input: input vector potential as a 2-column matrix (Φ, Ψ) with spherical harmonic
//!                   coefficients from order 1 to n. The maximum order is determined from the
//!                   number of rows in the matrix.
//| \param[out] out: output vector potential
template <class T0, class T1>
void rotation_coaxial_decomposition(t_real wavenumber, t_real tz,
                                    Eigen::MatrixBase<T0> const &input,
                                    Eigen::MatrixBase<T1> const &out) {
  using coefficient::a;
  auto const nr = input.rows();
  auto const with_n0 = std::abs(std::sqrt(nr) - std::lround(std::sqrt(nr))) <
                       std::abs(std::sqrt(nr + 1) - std::lround(std::sqrt(nr + 1)));
  t_int const N = std::lround(std::sqrt(with_n0 ? nr : nr + 1)) - 1;
  assert(input.cols() == 2);
  assert((with_n0 and (N + 1) * (N + 1) == input.rows()) or N * (N + 2) == input.rows());
  int const min_n = with_n0 ? 0 : 1;
  auto const index = with_n0 ? [](t_int n, t_int m) { return n * (n + 1) + m; } :
                               [](t_int n, t_int m) { return n * (n + 1) + m - 1; };
  assert(index(min_n, -min_n) == 0);
  assert(index(N, N) + 1 == input.rows());
  const_cast<Eigen::MatrixBase<T1> &>(out).resize(input.rows(), input.cols());
  auto const in_phi = [&input, N, index, min_n](t_int n, t_int m) -> t_complex {
    return n > N or std::abs(m) > n or n < min_n ? 0 : input(index(n, m), 0);
  };
  auto const in_psi = [&input, N, index, min_n](t_int n, t_int m) -> t_complex {
    return n > N or std::abs(m) > n or n < min_n ? 0 : input(index(n, m), 1);
  };
  auto const out_phi = [&out, index](t_int n, t_int m) -> t_complex & {
    return const_cast<Eigen::MatrixBase<T1> &>(out)(index(n, m), 0);
  };
  auto const out_psi = [&out, index](t_int n, t_int m) -> t_complex & {
    return const_cast<Eigen::MatrixBase<T1> &>(out)(index(n, m), 1);
  };
  if(with_n0)
    const_cast<Eigen::MatrixBase<T1> &>(out).row(0).fill(0);
  // n ≥ 1
  for(t_int n(1); n <= N; ++n) {
    auto const factor = tz * wavenumber / static_cast<t_real>(n * n + n);
    for(t_int m(-n); m <= n; ++m) {
      t_complex const cm(0, m * factor);
      auto const c0 = n * a<t_real>(n, m) * factor;
      auto const c1 = (n + 1) * a<t_real>(n - 1, m) * factor;
      out_phi(n, m) =
          in_phi(n, m) + cm * in_psi(n, m) + c0 * in_phi(n + 1, m) + c1 * in_phi(n - 1, m);
      out_psi(n, m) =
          in_psi(n, m) + cm * in_phi(n, m) + c0 * in_psi(n + 1, m) + c1 * in_psi(n - 1, m);
    }
  }
}

//! \brief Tranpose of the rotation-coaxial decomposition of the vector potentials
//! \param[in] wavenumber
//! \param[in] tz: translation alongst the z axis
//! \param[in] input: input vector potential as a 2-column matrix (Φ, Ψ) with spherical harmonic
//!                   coefficients from order 1 to n. The maximum order is determined from the
//!                   number of rows in the matrix.
//| \param[out] out: output vector potential
template <class T0, class T1>
void rotation_coaxial_decomposition_transpose(t_real wavenumber, t_real tz,
                                              Eigen::MatrixBase<T0> const &input,
                                              Eigen::MatrixBase<T1> const &out) {
  using coefficient::a;
  auto const nr = input.rows();
  auto const with_n0 = std::abs(std::sqrt(nr) - std::lround(std::sqrt(nr))) <
                       std::abs(std::sqrt(nr + 1) - std::lround(std::sqrt(nr + 1)));
  t_int const N = std::lround(std::sqrt(with_n0 ? nr : nr + 1)) - 1;
  assert((with_n0 and (N + 1) * (N + 1) == input.rows()) or N * (N + 2) == input.rows());
  auto const index = with_n0 ? [](t_int n, t_int m) { return n * (n + 1) + m; } :
                               [](t_int n, t_int m) { return n * (n + 1) + m - 1; };
  assert(index(with_n0 ? 0: 1, -with_n0 ? 0: 1) == 0);
  assert(index(N, N) + 1 == input.rows());
  const_cast<Eigen::MatrixBase<T1> &>(out).resize(input.rows(), input.cols());
  auto const in_phi = [&input, N, index](t_int n, t_int m) -> t_complex {
    return n > N or std::abs(m) > n or n <= 0 ? 0 : input(index(n, m), 0);
  };
  auto const in_psi = [&input, N, index](t_int n, t_int m) -> t_complex {
    return n > N or std::abs(m) > n or n <= 0 ? 0 : input(index(n, m), 1);
  };
  auto const out_phi = [&out, index](t_int n, t_int m) -> t_complex & {
    return const_cast<Eigen::MatrixBase<T1> &>(out)(index(n, m), 0);
  };
  auto const out_psi = [&out, index](t_int n, t_int m) -> t_complex & {
    return const_cast<Eigen::MatrixBase<T1> &>(out)(index(n, m), 1);
  };
  for(t_int n(1); n <= N; ++n) {
    auto const factor = tz * wavenumber / static_cast<t_real>(n * n + n);
    for(t_int m(-n); m <= n; ++m) {
      t_complex const cm(0, m * factor);
      auto const c0 = (n + 1) * a<t_real>(n - 1, m) * factor;
      auto const c1 = n * a<t_real>(n, m) * factor;
      out_phi(n, m) =
          in_phi(n, m) + cm * in_psi(n, m) + c0 * in_phi(n - 1, m) + c1 * in_phi(n + 1, m);
      out_psi(n, m) =
          in_psi(n, m) + cm * in_phi(n, m) + c0 * in_psi(n - 1, m) + c1 * in_psi(n + 1, m);
    }
  }
  if(with_n0) {
    out_phi(0, 0) = tz * wavenumber * a<t_real>(0, 0) * in_phi(1, 0);
    out_psi(0, 0) = tz * wavenumber * a<t_real>(0, 0) * in_psi(1, 0);
  }
}

//! Computes rotation-coaxial decomposition and returns result
template <class T0>
Matrix<typename T0::Scalar>
rotation_coaxial_decomposition(t_real wavenumber, t_real tz, Eigen::MatrixBase<T0> const &input) {
  Matrix<typename T0::Scalar> out = Matrix<typename T0::Scalar>::Zero(input.rows(), input.cols());
  rotation_coaxial_decomposition(wavenumber, tz, input, out);
  return out;
}

//! Computes rotation-coaxial decomposition transpose and returns result
template <class T0>
Matrix<typename T0::Scalar>
rotation_coaxial_decomposition_transpose(t_real wavenumber, t_real tz,
                                         Eigen::MatrixBase<T0> const &input) {
  Matrix<typename T0::Scalar> out = Matrix<typename T0::Scalar>::Zero(input.rows(), input.cols());
  rotation_coaxial_decomposition_transpose(wavenumber, tz, input, out);
  return out;
}
}
#endif
