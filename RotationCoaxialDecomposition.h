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
//! \param[in] in: input vector potential as a 2-column matrix (Φ, Ψ) with spherical harmonic
//!                coefficients from order 1 to n. The maximum order is determined from the number
//!                of rows in the matrix.
//| \param[out] out: output vector potential
template <class T0, class T1>
void rotation_coaxial_decomposition(t_real wavenumber, t_real tz, Eigen::MatrixBase<T0> const &in,
                                    Eigen::MatrixBase<T1> &out) {
  using coefficient::a;
  // solves nmax * (nmax + 2) = in.rows()
  t_int const nmax = std::lround(std::sqrt(in.rows()) - 1.0);
  assert(nmax * (nmax + 2) == in.rows());
  assert(nmax >= 1);
  out.resize(in.rows(), in.cols());
  auto const index = [](t_int n, t_int m) { return std::abs(m) > n ? 0 : n * (n + 1) + m - 1; };
  auto const in_phi = [&in, nmax, index](t_int n, t_int m) -> t_complex {
    return n > nmax ? 0 : in(index(n, m), 0);
  };
  auto const in_psi = [&in, nmax, index](t_int n, t_int m) -> t_complex {
    return n > nmax ? 0 : in(index(n, m), 1);
  };
  auto const out_phi = [&out, index](t_int n, t_int m) -> t_complex & {
    return out(index(n, m), 0);
  };
  auto const out_psi = [&out, index](t_int n, t_int m) -> t_complex & {
    return out(index(n, m), 1);
  };
  // n ≥ 1
  for(t_int n(1); n <= nmax; ++n) {
    for(t_int m(-n); m <= n; ++m) {
      auto const c0 = n * a<t_real>(n, m);
      auto const c1 = (n + 1) * a<t_real>(n - 1, m);
      // clang-format off
      out_phi(n, m) = in_phi(n, m) + tz / static_cast<t_real>(n * n + n) * (
          t_complex(0, wavenumber * wavenumber * m) * in_psi(n, m)
          + wavenumber * (c0 * in_phi(n + 1, m) + (n > 1 ? c1 * in_phi(n - 1, m): 0))
      );
      out_psi(n, m) = in_psi(n, m) + tz / static_cast<t_real>(n * n + n) * (
          t_complex(0, m) * in_phi(n, m)
          + wavenumber * (c0 * in_psi(n + 1, m) + (n > 1 ? c1 * in_psi(n - 1, m): 0))
      );
      // clang-format on
    }
  }
}

//! Computes rotation-coaxial decomposition and returns result
template <class T0>
Matrix<typename T0::Scalar>
rotation_coaxial_decomposition(t_real wavenumber, t_real tz, Eigen::MatrixBase<T0> const &in) {
  Matrix<typename T0::Scalar> out = Matrix<typename T0::Scalar>::Zero(in.rows(), in.cols());
  rotation_coaxial_decomposition(wavenumber, tz, in, out);
  return out;
}
}
#endif
