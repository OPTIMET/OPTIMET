#include "RotationCoefficients.h"
#include <algorithm>
#include <cmath>
#include <iterator>

namespace optimet {

RotationCoefficients::Real RotationCoefficients::a(t_uint n, t_int m) {
  t_uint const absm(std::abs(m));
  if(n < absm)
    return 0;
  return std::sqrt(static_cast<Real>((n + 1 + absm) * (n + 1 - absm)) /
                   static_cast<Real>((2 * n + 1) * (2 * n + 3)));
}

RotationCoefficients::Real RotationCoefficients::b(t_uint n, t_int m) {
  if(static_cast<t_uint>(std::abs(m)) > n)
    return 0;
  return (m >= 0 ? 1 : -1) * std::sqrt(static_cast<Real>((n - m - 1) * (n - m)) /
                                       static_cast<Real>((2 * n - 1) * (2 * n + 1)));
}

Eigen::Matrix<RotationCoefficients::Complex, 3, 1>
RotationCoefficients::factors(t_uint n, t_int m, t_int mu) const {
  auto const factor = std::exp(Complex(0, chi)) / b(n + 1, m - 1);
  auto const c0 =
      factor * 0.5 * b(n + 1, -mu - 1) * std::exp(Complex(0, phi)) * (1 - std::cos(theta));
  auto const c1 =
      -factor * 0.5 * b(n + 1, mu - 1) * std::exp(Complex(0, -phi)) * (1 + std::cos(theta));
  auto const c2 = -factor * a(n, mu) * std::sin(theta);
  return Eigen::Matrix<Complex, 3, 1>(c0, c1, c2);
}

RotationCoefficients::Complex RotationCoefficients::operator()(t_uint n, t_int m, t_int mu) {
  if(static_cast<t_uint>(std::abs(m)) > n or static_cast<t_uint>(std::abs(mu)) > n)
    return 0;
  if(m < 0)
    return std::conj(operator()(n, -m, -mu));

  auto const prior = cache.find(std::make_tuple(n, m, mu));
  if(prior != cache.end())
    return prior->second;

  auto const value = m == 0 ? initial(n, mu) : recursion(n, m, mu);
  cache[std::make_tuple(n, m, mu)] = value;
  return value;
}

void RotationCoefficients::coefficients(t_uint n, t_int m, t_int mu, std::vector<Complex> &coeffs) {

  if(m == 0) {
    coeffs.push_back(initial(n, mu));
    return;
  }
  if(m < 0) {
    auto const conjgate = coefficients(n, -m, -mu);
    std::transform(conjgate.begin(), conjgate.end(), std::back_inserter(coeffs),
                   [](Complex const &c) { return std::conj(c); });
    return;
  }
  auto const factors = this->factors(n, m, mu);
  auto const coeffs0 = coefficients(n + 1, m - 1, mu + 1);
  auto const coeffs1 = coefficients(n + 1, m - 1, mu - 1);
  auto const coeffs2 = coefficients(n + 1, m - 1, mu);
  std::transform(coeffs0.begin(), coeffs0.end(), std::back_inserter(coeffs),
                 [factors](Complex const &coeff) -> Complex { return factors(0) * coeff; });
  std::transform(coeffs1.begin(), coeffs1.end(), std::back_inserter(coeffs),
                 [factors](Complex const &coeff) { return factors(1) * coeff; });
  std::transform(coeffs2.begin(), coeffs2.end(), std::back_inserter(coeffs),
                 [factors](Complex const &coeff) { return factors(2) * coeff; });
}

RotationCoefficients::Complex RotationCoefficients::recursion(t_uint n, t_int m, t_int mu) {
  assert(n >= 1);
  auto const factors = this->factors(n, m, mu);
  return factors(0) * operator()(n + 1, m - 1, mu + 1) +
         factors(1) * operator()(n + 1, m - 1, mu - 1) + factors(2) * operator()(n + 1, m - 1, mu);
}
}
