#include "RotationCoefficients.h"
#include <algorithm>
#include <cmath>
#include <iterator>

namespace optimet {

RotationCoefficients::Real RotationCoefficients::a(t_uint n, t_int m) {
  t_uint const absm(std::abs(m));
  if(n < absm)
    return static_cast<Real>(0);
  return std::sqrt(static_cast<Real>((n + 1 + absm) * (n + 1 - absm)) /
                   static_cast<Real>((2 * n + 1) * (2 * n + 3)));
}

RotationCoefficients::Real RotationCoefficients::b(t_uint n, t_int m) {
  if(static_cast<t_uint>(std::abs(m)) > n)
    return static_cast<Real>(0);
  return (m >= 0 ? 1 : -1) * std::sqrt(static_cast<Real>((n - m - 1) * (n - m)) /
                                       static_cast<Real>((2 * n - 1) * (2 * n + 1)));
}

RotationCoefficients::Coefficients
RotationCoefficients::factors(t_uint n, t_int m, t_int mu) const {
  auto const factor = std::exp(Complex(0, chi_)) / b(n + 1, m - 1);
  auto const half_factor = static_cast<Real>(0.5) * factor;
  auto const c0 =
      half_factor * b(n + 1, -mu - 1) * std::exp(Complex(0, phi_)) * (1 - std::cos(theta_));
  auto const c1 =
      -half_factor * b(n + 1, mu - 1) * std::exp(Complex(0, -phi_)) * (1 + std::cos(theta_));
  auto const c2 = -factor * a(n, mu) * std::sin(theta_);
  return Coefficients(c0, c1, c2);
}

RotationCoefficients::Complex RotationCoefficients::with_caching(t_uint n, t_int m, t_int mu) {
  if(static_cast<t_uint>(std::abs(m)) > n or static_cast<t_uint>(std::abs(mu)) > n)
    return static_cast<Real>(0);
  if(m < 0)
    return std::conj(with_caching(n, -m, -mu));

  auto const prior = cache.find(std::make_tuple(n, m, mu));
  if(prior != cache.end())
    return prior->second;

  auto const value = recursion(n, m, mu);
  cache[std::make_tuple(n, m, mu)] = value;
  return value;
}

RotationCoefficients::Complex RotationCoefficients::recursion(t_uint n, t_int m, t_int mu) {
  if(static_cast<t_uint>(std::abs(m)) > n or static_cast<t_uint>(std::abs(mu)) > n)
    return static_cast<Real>(0);
  if(m < 0)
    return std::conj(with_caching(n, -m, -mu));

  if(m == 0)
    return initial(n, mu);

  auto const factors = this->factors(n, m, mu);
  return factors(0) * with_caching(n + 1, m - 1, mu + 1) +
         factors(1) * with_caching(n + 1, m - 1, mu - 1) +
         factors(2) * with_caching(n + 1, m - 1, mu);
}

Matrix<t_complex> RotationCoefficients::matrix(t_uint n) {
  Matrix<t_complex> result = Matrix<t_complex>::Zero(2 * n + 1, 2 * n + 1);
  for(t_uint i(0); i < 2 * n + 1; ++i) {
    t_int const m = static_cast<t_int>(i) - static_cast<t_int>(n);
    for(t_uint j(0); j < 2 * n + 1; ++j) {
      t_int const mu = static_cast<t_int>(j) - static_cast<t_int>(n);
      result(j, i) = operator()(n, m, mu);
    }
  }
  return result;
}

Rotation::Rotation(t_real const &theta, t_real const &phi, t_real const &chi, t_uint nmax)
    : theta_(theta), phi_(phi), chi_(chi), nmax_(nmax) {
  order.reserve(nmax+1);
  RotationCoefficients coeffs(theta, phi, chi);
  for(t_uint i(0); i <= nmax; ++i)
    order.push_back(coeffs.matrix(i));
}

Matrix<t_real> RotationCoefficients::basis_rotation(t_real theta, t_real phi, t_real chi) {
  Matrix<t_real> result(3, 3);
  result << -sin(phi) * sin(chi) - cos(theta) * cos(phi) * cos(chi),
         sin(phi) * cos(chi) - cos(theta) * cos(phi) * sin(chi), sin(theta) * cos(phi),
         cos(phi) * sin(chi) - cos(theta) * sin(phi) * cos(chi),
         -cos(phi) * cos(chi) - cos(theta) * sin(phi) * sin(chi), sin(theta) * sin(phi),
         sin(theta) * cos(chi), sin(theta) * sin(chi), cos(theta);
  return result;
};
}
