#include "TranslationAdditionCoefficients.h"
#include "Bessel.h"
#include "AJ_AuxFuns.h"
#include "constants.h"
#include "CompoundIterator.h"
#include <complex>
#include <iostream>

#include <boost/math/special_functions/legendre.hpp>

namespace optimet {
namespace {
constexpr t_int factorial(t_int n) { return n < 2 ? 1 : n * factorial(n - 1); }
} // anonymous namespace

t_complex Ynm(Spherical<t_real> const &R, t_int n, t_int m) {
  if(n < std::abs(m))
    return 0;
  auto const gamma =
      static_cast<t_real>(n * (n + 1) * (2 * n + 1)) /
      (constant::pi * static_cast<t_real>(4 * n * (n + 1))) *
      (static_cast<t_real>(factorial(n - m)) / static_cast<t_real>(factorial(n + m)));
  return std::sqrt(gamma) * std::exp(constant::i * (m * R.phi)) *
         boost::math::legendre_p(n, m, std::cos(R.the));
}

t_complex TranslationAdditionCoefficients::operator()(t_int n, t_int m, t_int l, t_int k) {
  if(std::abs(k) > l)
    return 0e0;

  // Check cache first
  t_indices const indices{{n, m, k, l}};
  auto const i_found = cache.find(indices);
  if(i_found != cache.end())
    return i_found->second;
  if(m < 0)
    throw std::runtime_error("Not implemented yet");
  if(n != 0)
    throw std::runtime_error("Not implemented yet");

  // Then check recursion that will be needed
  auto const result = n == 0 and m == 0 ? initial(l, k) : 0e0;

  cache[indices] = result;
  return result;
}

t_complex TranslationAdditionCoefficients::initial(t_int l, t_int k) {
  if(l == 0 and k == 0)
    return 1e0 / std::sqrt(4e0 * constant::pi);
  auto const wave = direction.rrr * waveK;
  auto const hb = std::get<0>(regular ? bessel<Bessel>(wave, l + 1) : bessel<Hankel1>(wave, l));
  auto const factor = std::sqrt(4e0 * constant::pi) * ((l + k) % 2 == 0 ? 1 : -1);
  return factor * Ynm(direction, l, -k) * hb[l];
}
}
