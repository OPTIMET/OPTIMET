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
}
