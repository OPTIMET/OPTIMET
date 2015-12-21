#include "TranslationAdditionCoefficients.h"
#include "Bessel.h"
#include "AJ_AuxFuns.h"
#include "constants.h"
#include "CompoundIterator.h"
#include <complex>
#include <iostream>
#include <cmath>

#include <boost/math/special_functions/legendre.hpp>

namespace optimet {
namespace {
constexpr bool is_valid(t_int n, t_int m) { return n >= 0 and std::abs(m) <= n; }
constexpr bool is_valid(t_int n, t_int m, t_int l, t_int k) {
  return is_valid(n, m) and is_valid(l, k);
}
constexpr t_int factorial(t_int n) { return n < 2 ? 1 : n * factorial(n - 1); }
inline t_real a_plus(t_int n, t_int m) {
  if(not is_valid(n, m))
    return 0e0;
  return std::sqrt(static_cast<t_real>((n + m + 1) * (n - m + 1)) /
                   static_cast<t_real>((2 * n + 1) * (2 * n + 3)));
}
inline t_real a_minus(t_int n, t_int m) {
  if(not is_valid(n, m))
    return 0e0;
  return std::sqrt(static_cast<t_real>((n + m) * (n - m)) /
                   static_cast<t_real>((2 * n + 1) * (2 * n - 1)));
}
inline t_real b_plus(t_int n, t_int m) {
  if(not is_valid(n, m))
    return 0e0;
  return std::sqrt(static_cast<t_real>((n + m + 2) * (n + m + 1)) /
                   static_cast<t_real>((2 * n + 1) * (2 * n + 3)));
}
inline t_real b_minus(t_int n, t_int m) {
  if(not is_valid(n, m))
    return 0e0;
  return std::sqrt(static_cast<t_real>((n - m) * (n - m - 1)) /
                   static_cast<t_real>((2 * n + 1) * (2 * n - 1)));
}
} // anonymous namespace

t_complex Ynm(Spherical<t_real> const &R, t_int n, t_int m) {
  if(not is_valid(n, m))
    return 0;
  auto const gamma =
      static_cast<t_real>(n * (n + 1) * (2 * n + 1)) /
      (constant::pi * static_cast<t_real>(4 * n * (n + 1))) *
      (static_cast<t_real>(factorial(n - m)) / static_cast<t_real>(factorial(n + m)));
  return std::sqrt(gamma) * std::exp(constant::i * (m * R.phi)) *
         boost::math::legendre_p(n, m, std::cos(R.the));
}

namespace details {
t_complex CachedRecurrence::operator()(t_int n, t_int m, t_int l, t_int k) {
  // It simplifies the recurrence if we assume zero outside the domain of validity
  if(not is_valid(n, m, l, k))
    return 0e0;

  // For simplicity, coefficients for negative m should be implemented separately using the symmetry
  // relationship.
  assert(m >= 0);

  // Check cache first
  t_indices const indices{{n, m, k, l}};
  auto const i_found = cache.find(indices);
  if(i_found != cache.end())
    return i_found->second;

  // Then check recursion that will be needed
  auto const result = impl(n, m, l, k);
  cache[indices] = result;
  return result;
}

t_complex CachedRecurrence::impl(t_int n, t_int m, t_int l, t_int k) {
  if(n == 0 and m == 0)
    return initial(l, k);
  else if(n == m)
    return diagonal_recurrence(n, l, k);
  else
    return offdiagonal_recurrence(n, m, l, k);
}

t_complex CachedRecurrence::initial(t_int l, t_int k) {
  if(l == 0 and k == 0)
    return 1e0 / std::sqrt(4e0 * constant::pi);
  auto const wave = direction.rrr * waveK;
  auto const hb = std::get<0>(regular ? bessel<Bessel>(wave, l + 1) : bessel<Hankel1>(wave, l + 1));
  auto const factor = std::sqrt(4e0 * constant::pi) * ((l + k) % 2 == 0 ? 1 : -1);
  assert(hb.size() > l + 1 and l >= 0);
  return factor * Ynm(direction, l, -k) * hb[l];
}

t_complex CachedRecurrence::diagonal_recurrence(t_int n, t_int l, t_int k) {
  return (operator()(n - 1, n - 1, l - 1, k - 1) * b_plus(l - 1, k - 1) +
          operator()(n - 1, n - 1, l + 1, k - 1) * b_minus(l + 1, k - 1)) /
         b_plus(n - 1, n - 1);
}

t_complex CachedRecurrence::offdiagonal_recurrence(t_int n, t_int m, t_int l, t_int k) {
  return (-operator()(n - 2, m, l, k) * a_minus(n - 1, m) +
          operator()(n - 1, m, l - 1, k) * a_plus(l - 1, k) +
          operator()(n - 1, m, l + 1, k) * a_minus(l + 1, k)) /
         a_plus(n - 1, m);
}
}
}
