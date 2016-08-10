#include "RotationCoefficients.h"
#include <cmath>

namespace optimet {
namespace {

t_real a(t_uint n, t_uint m) {
  if(n < m)
    return 0;
  return std::sqrt(static_cast<t_real>((n + 1 + m) * (n + 1 - m)) /
                   static_cast<t_real>((2 * n + 1) * (2 * n + 3)));
}
t_real a(t_uint n, t_int m) { return a(n, static_cast<t_uint>(std::abs(m))); }

t_real b(t_uint n, t_int m) {
  if(static_cast<t_uint>(std::abs(m)) > n)
    return 0;
  return (m >= 0 ? 1 : -1) * std::sqrt(static_cast<t_real>((n - m - 1) * (n - m)) /
                                       static_cast<t_real>((2 * n - 1) * (2 * n + 1)));
}
}
t_complex RotationCoefficients::operator()(t_uint n, t_int m, t_int mu) {
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

t_complex RotationCoefficients::recursion(t_uint n, t_int m, t_int mu) {
  assert(n >= 2);
  auto const factor = std::exp(t_complex(0, chi)) / b(n + 1, m - 1);
  auto const c0 =
      factor * 0.5 * b(n + 1, -mu - 1) * std::exp(t_complex(0, phi)) * (1 - std::cos(theta));
  auto const c1 =
      -factor * 0.5 * b(n + 1, mu - 1) * std::exp(t_complex(0, -phi)) * (1 + std::cos(theta));
  auto const c2 = -factor * a(n, mu) * std::sin(theta);
  return c0 * operator()(n + 1, m - 1, mu + 1) + c1 * operator()(n + 1, m - 1, mu - 1) +
         c2 * operator()(n + 1, m - 1, mu);
}
}
