#include "Bessel.h"
#include "CoAxialTranslationCoefficients.h"
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

CachedCoAxialRecurrence::Complex CachedCoAxialRecurrence::coeff(t_int n, t_int m, t_int l) {
  // It simplifies the recurrence if we assume zero outside the domain of
  // validity
  if(not is_valid(n, m, l, m))
    return static_cast<Real>(0);

  // the coaxial recurrence is independent of sign of m Gumerov (4.81)
  m = std::abs(m);

  // Check cache first
  t_indices const indices{{n, m, l}};
  auto const i_found = cache.find(indices);
  if(i_found != cache.end())
    return i_found->second;

  CachedCoAxialRecurrence::Complex const result = recurrence(n, m, l);
  cache[indices] = result;
  return result;
}

CachedCoAxialRecurrence::Complex CachedCoAxialRecurrence::recurrence(t_int n, t_int m, t_int l) {
  if(n == 0 and m == 0)
    return initial(l);
  else if(l < n) {
    CachedCoAxialRecurrence::Complex factor = static_cast<Complex>((l + n) % 2 == 0 ? 1 : -1);
    return coeff(l, m, n) * factor;
  } else if(m == n)
    return sectorial_recurrence(n, m, l);
  else if(m == 0)
    return zonal_recurrence(n, l);
  else if(m == n - 1)
    return sectorial_recurrence(n, m, l);
  else
    return offdiagonal_recurrence(n, m, l);
}

CachedCoAxialRecurrence::Complex CachedCoAxialRecurrence::initial(t_int l) {
  assert(l >= 0);
  CachedCoAxialRecurrence::Complex const wave = distance * waveK;
  auto const bessel = regular ? optimet::bessel<Bessel> : optimet::bessel<Hankel1>;
  CachedCoAxialRecurrence::Complex const hb = static_cast<CachedCoAxialRecurrence::Complex>(
      std::get<0>(bessel(static_cast<t_complex>(wave), l)).back());

  Real const factor = static_cast<Real>(std::sqrt(2 * l + 1) * (l % 2 == 0 ? 1 : -1));
  return factor * hb;
}

CachedCoAxialRecurrence::Complex
CachedCoAxialRecurrence::sectorial_recurrence(t_int n, t_int m, t_int l) {
  using coefficient::b;
  assert(l > 0 and n > 0 and l >= n and (m == n or m == n - 1) and (n + m != 1));
  // This formula requires bnm = 0 which is only true from m = n and m = n-1
  // It also requires bn-m to be non zero. This is zero if n+m = 1
  // Gumerov's b coeffs are equal to b_minus from Stout for m >=0 and
  // - b_minus for m < 0. Here m = n or n-1 by definition.
  return (coeff(n - 1, m - 1, l - 1) * b<Real>(l, -m) -
          coeff(n - 1, m - 1, l + 1) * b<Real>(l + 1, m - 1)) /
         b<Real>(n, -m);
}

CachedCoAxialRecurrence::Complex
CachedCoAxialRecurrence::offdiagonal_recurrence(t_int n, t_int m, t_int l) {
  // gumerov 4.80
  using coefficient::b;
  assert(m != 0 and n != 0 and m != n);
  return (coeff(n - 1, m - 1, l - 1) * b<Real>(l, -m) +
          coeff(n - 2, m, l) * (b<Real>(n - 1, m - 1)) -
          coeff(n - 1, m - 1, l + 1) * b<Real>(l + 1, m - 1)) /
         b<Real>(n, -m);
}

CachedCoAxialRecurrence::Complex CachedCoAxialRecurrence::zonal_recurrence(t_int n, t_int l) {
  // Gumerov 4.79 i.e. m = 0
  using coefficient::a;
  assert(l > 0 and n > 0 and l >= n);
  return (coeff(n - 1, 0, l - 1) * a<Real>(l - 1, 0) + coeff(n - 2, 0, l) * a<Real>(n - 2, 0) -
          coeff(n - 1, 0, l + 1) * a<Real>(l, 0)) /
         a<Real>(n - 1, 0);
}

CachedCoAxialRecurrence::Functor CachedCoAxialRecurrence::functor(t_int N) {
  // now assign them
  std::vector<t_complex> coefficients;
  for(auto n = 0; n <= N; ++n)
    for(auto m = -n; m <= n; ++m)
      for(auto l = std::abs(m); l <= N; ++l)
        coefficients.push_back(operator()(n, m, l));
  return Functor(N, std::move(coefficients));
}
}
