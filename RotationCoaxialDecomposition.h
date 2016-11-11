#ifndef OPTIMET_ROTATION_COAXIAL_DECOMPOSITION_H
#define OPTIMET_ROTATION_COAXIAL_DECOMPOSITION_H

#include "Coefficients.h"
#include "Types.h"
#include "constants.h"

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
                                    Eigen::MatrixBase<T1> &out) {
  using coefficient::a;
  auto const nr = input.rows();
  auto const with_n0 = std::abs(std::sqrt(nr) - std::lround(std::sqrt(nr))) <
                       std::abs(std::sqrt(nr + 1) - std::lround(std::sqrt(nr + 1)));
  t_int const N = std::lround(std::sqrt(with_n0 ? nr : nr + 1)) - 1;
  assert((with_n0 and (N + 1) * (N + 1) == input.rows()) or N * (N + 2) == input.rows());
  int const min_n = with_n0 ? 0 : 1;
  auto const index = with_n0 ? [](t_int n, t_int m) { return n * (n + 1) + m; } :
                               [](t_int n, t_int m) { return n * (n + 1) + m - 1; };
  assert(index(min_n, -min_n) == 0);
  assert(index(N, N) + 1 == input.rows());
  out.resize(input.rows(), input.cols());
  auto const in_phi = [&input, N, index](t_int n, t_int m) -> t_complex {
    return n > N or std::abs(m) > n ? 0 : input(index(n, m), 0);
  };
  auto const in_psi = [&input, N, index](t_int n, t_int m) -> t_complex {
    return n > N or std::abs(m) > n ? 0 : input(index(n, m), 1);
  };
  auto const out_phi = [&out, index](t_int n, t_int m) -> t_complex & {
    return out(index(n, m), 0);
  };
  auto const out_psi = [&out, index](t_int n, t_int m) -> t_complex & {
    return out(index(n, m), 1);
  };
  if(with_n0)
    out.row(0).fill(0);
  // n ≥ 1
  for(t_int n(1); n <= N; ++n) {
    for(t_int m(-n); m <= n; ++m) {
      auto const c0 = n * a<t_real>(n, m);
      auto const c1 = (n + 1) * a<t_real>(n - 1, m);
      // clang-format off
      out_phi(n, m) = in_phi(n, m) + tz / static_cast<t_real>(n * n + n) * (
          t_complex(0, wavenumber * wavenumber * m) * in_psi(n, m)
          + wavenumber * (c0 * in_phi(n + 1, m) + (n > min_n ? c1 * in_phi(n - 1, m): 0))
      );
      out_psi(n, m) = in_psi(n, m) + tz / static_cast<t_real>(n * n + n) * (
          t_complex(0, m) * in_phi(n, m)
          + wavenumber * (c0 * in_psi(n + 1, m) + (n > min_n ? c1 * in_psi(n - 1, m): 0))
      );
      // clang-format on
    }
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
}
#endif
