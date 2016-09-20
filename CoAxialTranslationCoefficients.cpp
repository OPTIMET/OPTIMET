#include "CoAxialTranslationCoefficients.h"
#include "Bessel.h"
#include "constants.h"
#include <Coefficients.h>
#include <cmath>
#include <complex>

#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>

namespace optimet {
namespace {
//! True if n and m within shperical harmonics validity regime
constexpr bool is_valid(t_int n, t_int m) { return n >= 0 and std::abs(m) <= n; }
//! True if both pairs are valid
constexpr bool is_valid(t_int n, t_int m, t_int l, t_int k) {
  return is_valid(n, m) and is_valid(l, k);
}
}
namespace details {
t_complex CachedCoAxialRecurrence::operator()(t_int n, t_int m, t_int l) {
  // It simplifies the recurrence if we assume zero outside the domain of
  // validity
  if(not is_valid(n, m, l, m))
    return 0e0;

  // For simplicity, coefficients for negative m should be implemented
  // separately using the symmetry
  // relationship.
  assert(m >= 0);

  // Check cache first
  t_indices const indices{{n, m, l}};
  auto const i_found = cache.find(indices);
  if(i_found != cache.end())
    return i_found->second;

  auto const result = recurrence(n, m, l);
  cache[indices] = result;
  return result;
}

t_complex CachedCoAxialRecurrence::recurrence(t_int n, t_int m, t_int l) {
  if(n == 0 and m == 0)
    return initial(l);
  else if(l < n) {
    t_complex factor = ((l + n) % 2 == 0 ? 1 : -1);
    return operator()(l, m, n) * factor;
  } else if(m == n)
    return sectorial_recurrence(n, m, l);
  else if(m == 0)
    return zonal_recurrence(n, l);
  else if(m == n - 1)
    return sectorial_recurrence(n, m, l);
  else
    return offdiagonal_recurrence(n, m, l);
}

t_complex CachedCoAxialRecurrence::initial(t_int l) {
  assert(l >= 0);
  auto const wave = direction.rrr * waveK;
  auto const bessel = regular ? optimet::bessel<Bessel> : optimet::bessel<Hankel1>;
  auto const hb = std::get<0>(bessel(wave, l)).back();

  auto const factor = std::sqrt(2 * l + 1) * (l % 2 == 0 ? 1 : -1);

  return factor * hb;
}

t_complex CachedCoAxialRecurrence::sectorial_recurrence(t_int n, t_int m, t_int l) {
  using coefficient::b;
  assert(l > 0 and n > 0 and l >= n and (m == n or m == n - 1) and (n + m != 1));
  // This formula requires bnm = 0 which is only true from m = n and m = n-1
  // It also requires bn-m to be non zero. This is zero if n+m = 1
  // Gumerov's b coeffs are equal to b_minus from Stout for m >=0 and
  // - b_minus for m < 0. Here m = n or n-1 by definition.
  return (operator()(n - 1, m - 1, l - 1) * b(l, -m) -
          operator()(n - 1, m - 1, l + 1) * b(l + 1, m - 1)) /
         b(n, -m);
}

t_complex CachedCoAxialRecurrence::offdiagonal_recurrence(t_int n, t_int m, t_int l) {
  // gumerov 4.80
  using coefficient::b;
  assert(m != 0 and n != 0 and m != n);
  return (+operator()(n - 1, m - 1, l - 1) * b(l, -m) +
          operator()(n - 2, m, l) * (b(n - 1, m - 1)) -
          operator()(n - 1, m - 1, l + 1) * b(l + 1, m - 1)) /
         b(n, -m);
}

t_complex CachedCoAxialRecurrence::zonal_recurrence(t_int n, t_int l) {
  // Gumerov 4.79 i.e. m = 0
  using coefficient::a;
  assert(l > 0 and n > 0 and l >= n);
  return (operator()(n - 1, 0, l - 1) * a(l - 1, 0) //
          + operator()(n - 2, 0, l) * a(n - 2, 0)   //
          - operator()(n - 1, 0, l + 1) * a(l, 0)) /
         a(n - 1, 0);
}
}
t_complex CoAxialTranslationAdditionCoefficients::operator()(t_int n, t_int m, t_int l) {
  // the coaxial recurrence is independent of sign of m Gumerov (4.81)
  return cached_recurrence(n, std::abs(m), l);
}
}
