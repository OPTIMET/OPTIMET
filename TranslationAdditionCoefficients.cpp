#include "TranslationAdditionCoefficients.h"
#include "Bessel.h"
#include "constants.h"

#include <complex>
#include <cmath>

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
//! Eases computing ratios of two factorials
inline t_real factorial_ratio(t_int n, t_int m) {
  return n == m ? 1 : std::tgamma(n + 1) / std::tgamma(m + 1);
}
//! Coefficient of Stout (2004) Appendix C recurrence relationship
inline t_real a_plus(t_int n, t_int m) {
  if(not is_valid(n, m))
    return 0e0;
  return -std::sqrt(static_cast<t_real>((n + m + 1) * (n - m + 1)) /
                    static_cast<t_real>((2 * n + 1) * (2 * n + 3)));
}
//! Coefficient of Stout (2004) Appendix C recurrence relationship
inline t_real a_minus(t_int n, t_int m) {
  if(not is_valid(n, m))
    return 0e0;
  return std::sqrt(static_cast<t_real>((n + m) * (n - m)) /
                   static_cast<t_real>((2 * n + 1) * (2 * n - 1)));
}
//! Coefficient of Stout (2004) Appendix C recurrence relationship
inline t_real b_plus(t_int n, t_int m) {
  if(not is_valid(n, m))
    return 0e0;
  return std::sqrt(static_cast<t_real>((n + m + 2) * (n + m + 1)) /
                   static_cast<t_real>((2 * n + 1) * (2 * n + 3)));
}
//! Coefficient of Stout (2004) Appendix C recurrence relationship
inline t_real b_minus(t_int n, t_int m) {
  if(not is_valid(n, m))
    return 0e0;
  return std::sqrt(static_cast<t_real>((n - m) * (n - m - 1)) /
                   static_cast<t_real>((2 * n + 1) * (2 * n - 1)));
}
} // anonymous namespace

t_complex Ynmlegacy(Spherical<t_real> const &R, t_int n, t_int m) {
  if(not is_valid(n, m))
    return 0;
  auto const gamma = static_cast<t_real>(2 * n + 1) / (constant::pi * static_cast<t_real>(4)) *
                     factorial_ratio(n - m, n + m);
  return std::sqrt(gamma) * std::exp(constant::i * (m * R.phi)) *
         boost::math::legendre_p(n, m, std::cos(R.the));
}

t_complex Ynm(Spherical<t_real> const &R, t_int n, t_int m) {
  if(not is_valid(n, m))
    return 0;
  return boost::math::spherical_harmonic(n, m, R.the, R.phi);
}

namespace details {
t_complex CachedRecurrence::operator()(t_int n, t_int m, t_int l, t_int k) {
  // It simplifies the recurrence if we assume zero outside the domain of
  // validity
  if(not is_valid(n, m, l, k))
    return 0e0;

  // For simplicity, coefficients for negative m should be implemented
  // separately using the symmetry
  // relationship.
  assert(m >= 0);

  // Check cache first
  t_indices const indices{{n, m, l, k}};
  auto const i_found = cache.find(indices);
  if(i_found != cache.end())
    return i_found->second;

  auto const result = recurrence(n, m, l, k);
  cache[indices] = result;
  return result;
}

t_complex CachedRecurrence::recurrence(t_int n, t_int m, t_int l, t_int k) {
  if(n == 0 and m == 0)
    return initial(l, k);
  else if(n == m)
    return diagonal_recurrence(n, l, k);
  else
    return offdiagonal_recurrence(n, m, l, k);
}

t_complex CachedRecurrence::initial(t_int l, t_int k) {
  assert(l >= 0);
  auto const wave = direction.rrr * waveK;
  auto const bessel = regular ? optimet::bessel<Bessel> : optimet::bessel<Hankel1>;
  auto const hb = std::get<0>(bessel(wave, l)).back();
  if(l == 0 and k == 0)
    return hb;
  auto const factor = std::sqrt(4e0 * constant::pi) * ((l + k) % 2 == 0 ? 1 : -1);
  return factor * Ynm(direction, l, -k) * hb;
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
  }
  else if(m == n)
    return sectorial_recurrence(n, m, l);
  else if(m==0)
    return zonal_recurrence(n, l);
  else if(m == n-1)
    return sectorial_recurrence(n, m, l);
  else
    return offdiagonal_recurrence(n, m, l);
}

t_complex CachedCoAxialRecurrence::initial(t_int l) {
  assert(l >= 0);
  auto const wave = direction.rrr * waveK;
  auto const bessel = regular ? optimet::bessel<Bessel> : optimet::bessel<Hankel1>;
  auto const hb = std::get<0>(bessel(wave, l)).back();

  auto const factor = std::sqrt(2*l+1) * (l % 2 == 0 ? 1 : -1);

  return factor * hb;
}

t_complex CachedCoAxialRecurrence::sectorial_recurrence(t_int n, t_int m, t_int l) {
  assert(l > 0 and n > 0 and l >= n and (m == n or m == n-1) and (n+m != 1));
  // This formula requires bnm = 0 which is only true from m = n and m = n-1
  // It also requires bn-m to be non zero. This is zero if n+m = 1
  // Gumerov's b coeffs are equal to b_minus from Stout for m >=0 and
  // - b_minus for m < 0. Here m = n or n-1 by definition.
  return (operator()(n-1, m-1, l-1) * b_minus(l, -m) +
          operator()(n-1, m-1, l+1) * b_minus(l+1, m-1)) /
          b_minus(n, -m);
}

t_complex CachedCoAxialRecurrence::offdiagonal_recurrence(t_int n, t_int m, t_int l) {
  // gumerov 4.80 
  assert(m != 0 and n != 0 and m != n);
  return (operator()(n-1, m-1, l-1) * b_minus(l, -m) - // minus because -m < 0 to match def
         operator()(n-2, m, l) * (b_minus(n-1, m-1)) +  // sign ok because m >= 1
         operator()(n-1, m-1, l+1) * b_minus(l+1, m-1)) /
         b_minus(n,-m);
}

t_complex CachedCoAxialRecurrence::zonal_recurrence(t_int n, t_int l) {
  // Gumerov 4.79 i.e. m = 0
  assert(l > 0 and n > 0 and l >= n);
  return (- operator()(n-1, 0, l-1) * a_plus(l-1, 0) +
          - operator()(n-2, 0, l) * a_plus(n-2, 0) +
          + operator()(n-1, 0, l+1) * a_plus(l, 0)) /
          - a_plus(n-1, 0);
}

} // end of detailed namespace

t_complex TranslationAdditionCoefficients::operator()(t_int n, t_int m, t_int l, t_int k) {
  if(m >= 0)
    return positive(n, m, l, k);

  auto const sign = negative.is_regular() ? k + m : n + m + l + k;
  auto const result = std::conj(negative(n, -m, l, -k));
  return sign % 2 == 0 ? result : -result;
}

t_complex CoAxialTranslationAdditionCoefficients::operator()(t_int n, t_int m, t_int l) {
  // the coaxial recurrence is independent of sign of m Gumerov (4.81)
  return cached_recurrence(n, std::abs(m), l);
}
} // end of optimet namespace
