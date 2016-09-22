#include "catch.hpp"

#include "Bessel.h"
#include "CoAxialTranslationCoefficients.h"
#include "constants.h"
#include <Coefficients.h>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cmath>
#include <iostream>
#include <random>

extern std::unique_ptr<std::mt19937_64> mersenne;

using namespace optimet;

void check_coaxial_n_recurrence(CoAxialTranslationAdditionCoefficients tca, t_int n, t_int m,
                                t_int l) {
  // This aims to check that we correctly implement formula 4.79
  using coefficient::a;
  auto left0 = a(n - 1, m) * tca(n - 1, m, l) - a(n, m) * tca(n + 1, m, l);
  auto right0 = a(l, m) * tca(n, m, l + 1) - a(l - 1, m) * tca(n, m, l - 1);
  INFO("Testing n recurrence for n " << n << " m " << m << " l " << l);
  CHECK(left0.real() == Approx(right0.real()));
  CHECK(left0.imag() == Approx(right0.imag()));
}

void check_coaxial_m_recurrence(CoAxialTranslationAdditionCoefficients tca, t_int n, t_int m,
                                t_int l) {
  // This aims to check that we correctly implement formula 4.80
  using coefficient::b;
  auto left0 = b(n, m) * tca(n - 1, m + 1, l) - b(n + 1, -m - 1) * tca(n + 1, m + 1, l);
  auto right0 = b(l + 1, m) * tca(n, m, l + 1) - b(l, -m - 1) * tca(n, m, l - 1);
  INFO("Testing m recurrence for n " << n << " m " << m << " l " << l);
  CHECK(left0.real() == Approx(right0.real()));
  CHECK(left0.imag() == Approx(right0.imag()));
}

void check_coaxial_mn_recurrence(CoAxialTranslationAdditionCoefficients tca, t_int n, t_int m,
                                 t_int l) {
  // This aims to check that we correctly implement formula 4.84
  using coefficient::b;
  auto left0 = b(n + 1, -m - 1) * tca(n + 1, m + 1, l);
  auto right0 = b(l, -m - 1) * tca(n, m, l - 1) - b(l + 1, m) * tca(n, m, l + 1);
  INFO("Testing m=n recurrence for n " << n << " m " << m << " l " << l);
  CHECK(left0.real() == Approx(right0.real()));
  CHECK(left0.imag() == Approx(right0.imag()));
}

void check_coaxial_m_symmetry(CoAxialTranslationAdditionCoefficients tca, t_int n, t_int m,
                              t_int l) {
  // check that value is independent of m
  auto left0 = tca(n, m, l);
  auto right0 = tca(n, -m, l);
  INFO("Testing m symmetry for " << n << " m " << m << " l " << l);
  CHECK(left0.real() == Approx(right0.real()));
  CHECK(left0.imag() == Approx(right0.imag()));
}

void check_coaxial_ln_symmetry(CoAxialTranslationAdditionCoefficients tca, t_int n, t_int m,
                               t_int l) {
  // check that value has the right sign change for m<->l
  auto left0 = tca(n, m, l);
  t_complex sign = ((n + l) % 2 == 0 ? 1 : -1);
  auto right0 = tca(l, m, n) * sign;
  INFO("Testing l to n symmetry n " << n << " m " << m << " l " << l);
  CHECK(left0.real() == Approx(right0.real()));
  CHECK(left0.imag() == Approx(right0.imag()));
}

TEST_CASE("CoAxial") {
  SECTION("Initial values") {
    using coefficient::a;
    using coefficient::b;
    t_real const R(1e0);
    t_complex const waveK(1e0, 1.5e0);
    CoAxialTranslationAdditionCoefficients tca(R, waveK, true);
    // Numbers are generated from the same formula in Scipy
    CHECK(tca(0, 0, 0).real() == Approx(1.1400511799225792));
    CHECK(tca(0, 0, 0).imag() == Approx(-0.55962217045848206));
    CHECK(tca(0, 0, 4).real() == Approx(-0.028191522402192234));
    CHECK(tca(0, 0, 4).imag() == Approx(-0.02162885905593049));
    CHECK(tca(1, 0, 1).real() == Approx(1.2274819687880665));
    CHECK(tca(1, 0, 1).imag() == Approx(-1.0271756758800463));
    CHECK(tca(-1, 1, 3).real() == Approx(0));
    CHECK(tca(-1, 1, 3).imag() == Approx(0));
    CHECK(tca(1, 0, -1).real() == Approx(0));
    CHECK(tca(1, 0, -1).imag() == Approx(0));
    CHECK(tca(1, 1, 3).real() == Approx(-0.085169586217943016));
    CHECK(tca(1, 1, 3).imag() == Approx(0.36331568009355053));

    t_complex expected;
    expected =
        (tca(1, 0, 1) * a(1, 0) + tca(0, 0, 2) * a(0, 0) - tca(1, 0, 3) * a(2, 0)) / (a(1, 0));
    CHECK(tca(2, 0, 2).real() == Approx(expected.real()));
    CHECK(tca(2, 0, 2).imag() == Approx(expected.imag()));

    expected = (tca(0, 0, 2) * a(2, 0) - tca(0, 0, 4) * a(3, 0)) / (a(0, 0));
    CHECK(tca(1, 0, 3).real() == Approx(expected.real()));
    CHECK(tca(1, 0, 3).imag() == Approx(expected.imag()));

    expected =
        (tca(1, 0, 3) * a(3, 0) + tca(0, 0, 4) * a(0, 0) - tca(1, 0, 5) * a(4, 0)) / (a(1, 0));
    CHECK(tca(2, 0, 4).real() == Approx(expected.real()));
    CHECK(tca(2, 0, 4).imag() == Approx(expected.imag()));

    expected = (-tca(0, 0, 2) * b(3, -1) + tca(0, 0, 4) * b(4, 0)) / (-b(1, -1));
    CHECK(tca(1, 1, 3).real() == Approx(expected.real()));
    CHECK(tca(1, 1, 3).imag() == Approx(expected.imag()));

    expected =
        (-tca(2, 0, 2) * b(3, -1) - tca(1, 1, 3) * b(2, 0) + tca(2, 0, 4) * b(4, 0)) / (-b(3, -1));
    CHECK(tca(3, 1, 3).real() == Approx(expected.real()));
    CHECK(tca(3, 1, 3).imag() == Approx(expected.imag()));
  }
  SECTION("Check n and m recurrence") {
    t_complex const waveK(1e0, 1.5e0);
    CoAxialTranslationAdditionCoefficients tca(1.0, waveK, true);
    t_int max_recur = 10;
    for(t_int l = 0; l < max_recur; ++l) {
      for(t_int n = 0; n < max_recur; ++n) {
        check_coaxial_mn_recurrence(tca, n, n, l);
        for(t_int m = -n; m <= n; ++m) {
          check_coaxial_ln_symmetry(tca, n, m, l);
          check_coaxial_m_symmetry(tca, n, m, l);
          check_coaxial_n_recurrence(tca, n, m, l);
          check_coaxial_m_recurrence(tca, n, m, l);
        }
      }
    }
  }
}

void check_coaxial_translation_onaxis(t_real expansion_pos, t_real reexpansion_pos,
                                      bool expansion_regular, bool reexpansion_regular, t_int n,
                                      t_int m, t_complex waveK) {
  assert(!(expansion_regular and !reexpansion_regular));
  assert(m <= n);
  bool coeff_regular = expansion_regular == reexpansion_regular;
  t_real const translation(expansion_pos - reexpansion_pos);
  CoAxialTranslationAdditionCoefficients tca(translation, waveK, coeff_regular);
  auto const basis_func = expansion_regular ? optimet::bessel<Bessel> : optimet::bessel<Hankel1>;
  auto const re_basis_func =
      reexpansion_regular ? optimet::bessel<Bessel> : optimet::bessel<Hankel1>;
  t_complex translated = 0;
  for(t_int l = std::abs(m); l < std::abs(m) + 25; l++) {
    translated += tca(n, m, l) * std::get<0>(re_basis_func(reexpansion_pos * waveK, l)).back() *
                  boost::math::spherical_harmonic(l, m, 0, 0);
  }
  auto expected = std::get<0>(basis_func(expansion_pos * waveK, n)).back() *
                  boost::math::spherical_harmonic(n, m, 0, 0);
  INFO("Testing translation for n: " << n << " m: " << m << " regular: " << expansion_regular
                                     << " to: " << reexpansion_regular);
  CHECK(expected.real() == Approx(translated.real()));
  CHECK(expected.imag() == Approx(translated.imag()));
}

void check_coaxial_translation_off_axis_reexpand_iregular(t_real r_p, t_real r_q, t_int n, t_int m,
                                                          t_complex waveK) {
  assert(m <= n);
  std::uniform_real_distribution<> theta_dist(0, consPi);
  std::uniform_real_distribution<> phi_dist(0, 2 * consPi);
  t_real theta_p = consPi / 10;
  t_real phi = phi_dist(*mersenne);
  t_real theta_q = asin(sin(theta_p) * r_p / r_q);
  t_real z1 = cos(theta_p) * r_p;
  t_real z2 = cos(theta_q) * r_q;
  t_real const translation(z1 - z2);
  std::cout << "translation " << translation << " z1 " << z1 << " theta_p " << theta_p << " z2 "
            << z2 << " theta_q " << theta_q << " r_p " << r_p << " r_q " << r_q << std::endl;
  CoAxialTranslationAdditionCoefficients tca(translation, waveK, false);
  auto const basis_func = optimet::bessel<Hankel1>;
  auto const re_basis_func = optimet::bessel<Bessel>;
  t_complex translated = 0;
  for(t_int l = std::abs(m); l < std::abs(m) + 25; l++) {
    translated += tca(n, m, l) * std::get<0>(re_basis_func(r_q * waveK, l)).back() *
                  boost::math::spherical_harmonic(l, m, theta_q, phi);
  }
  auto expected = std::get<0>(basis_func(r_p * waveK, n)).back() *
                  boost::math::spherical_harmonic(n, m, theta_p, phi);
  INFO("Testing translation for n: " << n << " m: " << m);
  CHECK(expected.real() == Approx(translated.real()));
  CHECK(expected.imag() == Approx(translated.imag()));
}

TEST_CASE("Coaxial translation") {
  std::uniform_real_distribution<> small_dist(0, 1);
  std::uniform_real_distribution<> large_dist(10, 50);
  std::uniform_real_distribution<> wave_dist(0, 5);
  t_real waver = wave_dist(*mersenne);
  t_real wavei = wave_dist(*mersenne);
  t_complex waveK(waver, wavei);
  for(t_int n = 0; n < 10; n++) {
    for(t_int m = -n; m <= n; m++) {
      // "Simple singular expanded in regular"
      t_real small = small_dist(*mersenne);
      t_real large = large_dist(*mersenne);
      t_real large_small_diff = large - small;
      check_coaxial_translation_onaxis(large, small, false, true, n, m, waveK);
      // "Simple singular expanded in singular"
      check_coaxial_translation_onaxis(large, large_small_diff, false, false, n, m, waveK);
      // "Simple regular expanded in regular"
      check_coaxial_translation_onaxis(large, large_small_diff, true, true, n, m, waveK);
      // "Simple regular expanded in regular"
      check_coaxial_translation_onaxis(large, small, true, true, n, m, waveK);
    }
  }
  check_coaxial_translation_off_axis_reexpand_iregular(10, 4, 0, 0, waveK);
}
