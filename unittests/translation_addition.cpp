#include "catch.hpp"

#include "Bessel.h"
#include "TranslationAdditionCoefficients.h"
#include "constants.h"
#include <Coefficients.h>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <random>

extern std::unique_ptr<std::mt19937_64> mersenne;

using namespace optimet;

TEST_CASE("Check Ynm") {
  using namespace boost::math;
  Spherical<t_real> const R(1e0, 0.42, 0.36);
  t_complex const waveK(1e0, 1.5e0);

  auto const Y00 = std::sqrt(1e0 / 4e0 / constant::pi) * legendre_p(0, 0, std::cos(R.the));
  auto const Y10 = std::sqrt(6e0 / 8e0 / constant::pi) * legendre_p(1, 0, std::cos(R.the));
  auto const Y11 = std::sqrt(6e0 / 16e0 / constant::pi) * legendre_p(1, 1, std::cos(R.the)) *
                   std::exp(constant::i * R.phi);
  auto const Y1m1 = std::sqrt(12e0 / 8e0 / constant::pi) * legendre_p(1, -1, std::cos(R.the)) *
                    std::exp(-constant::i * R.phi);
  auto const Y2m1 = std::sqrt(180e0 / 24e0 / constant::pi) * legendre_p(2, -1, std::cos(R.the)) *
                    std::exp(-constant::i * R.phi);
  auto const Y32 = std::sqrt(84e0 / 5760e0 / constant::pi) * legendre_p(3, 2, std::cos(R.the)) *
                   std::exp(constant::i * 2e0 * R.phi);
  auto const Y3m2 = std::sqrt(10080e0 / 48e0 / constant::pi) * legendre_p(3, -2, std::cos(R.the)) *
                    std::exp(-constant::i * 2e0 * R.phi);
  auto const Y40 = std::sqrt(4320e0 / 1920e0 / constant::pi) * legendre_p(4, 0, std::cos(R.the));

  CHECK(std::abs(Y00 - Ynm(R, 0, 0)) == Approx(0));
  CHECK(std::abs(Y10 - Ynm(R, 1, 0)) == Approx(0));
  CHECK(std::abs(Y11 - Ynm(R, 1, 1)) == Approx(0));
  CHECK(std::abs(Y1m1 - Ynm(R, 1, -1)) == Approx(0));
  CHECK(std::abs(Y2m1 - Ynm(R, 2, -1)) == Approx(0));
  CHECK(std::abs(Y32 - Ynm(R, 3, 2)) == Approx(0));
  CHECK(std::abs(Y3m2 - Ynm(R, 3, -2)) == Approx(0));
  CHECK(std::abs(Y40 - Ynm(R, 4, 0)) == Approx(0));
  CHECK(std::abs(Y00 - Ynmlegacy(R, 0, 0)) == Approx(0));
  CHECK(std::abs(Y10 - Ynmlegacy(R, 1, 0)) == Approx(0));
  CHECK(std::abs(Y11 - Ynmlegacy(R, 1, 1)) == Approx(0));
  CHECK(std::abs(Y1m1 - Ynmlegacy(R, 1, -1)) == Approx(0));
  CHECK(std::abs(Y2m1 - Ynmlegacy(R, 2, -1)) == Approx(0));
  CHECK(std::abs(Y32 - Ynmlegacy(R, 3, 2)) == Approx(0));
  CHECK(std::abs(Y3m2 - Ynmlegacy(R, 3, -2)) == Approx(0));
  CHECK(std::abs(Y40 - Ynmlegacy(R, 4, 0)) == Approx(0));

  std::uniform_real_distribution<> rdist(0, constant::pi);
  std::uniform_int_distribution<> ndist(1, 50);
  std::uniform_real_distribution<> mdist(-1, 1);
  Spherical<t_real> R2;
  for(int i = 0; i < 10; ++i) {
    auto const theta = rdist(*mersenne);
    auto const phi = 2 * rdist(*mersenne);
    R2.the = theta;
    R2.phi = phi;
    auto const n = ndist(*mersenne);
    auto const m = int(round(n * mdist(*mersenne)));
    CHECK(std::abs(Ynm(R2, n, m) - Ynmlegacy(R2, n, m)) == Approx(0));
  }
}

template <class RECURRENCE>
void check_recurrence(Spherical<t_real> const &R, t_complex const &waveK, bool is_regular) {
  using namespace boost::math;
  RECURRENCE ta(R, waveK, is_regular);

  SECTION("Check k > l always zero") {
    CHECK(std::abs(ta(0, 0, 0, 1)) == Approx(0));
    CHECK(std::abs(ta(0, 0, 0, -1)) == Approx(0));
    CHECK(std::abs(ta(0, 0, 1, 2)) == Approx(0));
    CHECK(std::abs(ta(0, 0, 1, -2)) == Approx(0));
    CHECK(std::abs(ta(2, 1, 1, 2)) == Approx(0));
    CHECK(std::abs(ta(2, 1, 1, -2)) == Approx(0));
  }

  SECTION("Check initial conditions") {
    auto const bess = is_regular ? bessel<Bessel> : bessel<Hankel1>;
    auto const hb = std::get<0>(bess(R.rrr * waveK, 4));
    auto const Y10 = std::sqrt(6e0 / 8e0 / constant::pi) * legendre_p(1, 0, std::cos(R.the));
    auto const Y1m1 = std::sqrt(12e0 / 8e0 / constant::pi) * legendre_p(1, -1, std::cos(R.the)) *
                      std::exp(-constant::i * R.phi);
    auto const Y3m2 = std::sqrt(10080e0 / 48e0 / constant::pi) *
                      legendre_p(3, -2, std::cos(R.the)) * std::exp(-constant::i * 2e0 * R.phi);
    auto const Y32 = std::sqrt(84e0 / 5760e0 / constant::pi) * legendre_p(3, 2, std::cos(R.the)) *
                     std::exp(constant::i * 2e0 * R.phi);
    auto const Y20 = std::sqrt(60e0 / 48e0 / constant::pi) * legendre_p(2, 0, std::cos(R.the));
    auto const Y40 = std::sqrt(4320e0 / 1920e0 / constant::pi) * legendre_p(4, 0, std::cos(R.the));
    auto const factor = std::sqrt(4e0 * constant::pi);
    CHECK(std::abs(ta(0, 0, 0, 0) - hb[0]) == Approx(0));
    CHECK(std::abs(ta(0, 0, 1, 0) + factor * Y10 * hb[1]) == Approx(0));
    CHECK(std::abs(ta(0, 0, 1, 1) - factor * Y1m1 * hb[1]) == Approx(0));
    CHECK(std::abs(ta(0, 0, 3, 2) + factor * Y3m2 * hb[3]) == Approx(0));
    CHECK(std::abs(ta(0, 0, 3, -2) + factor * Y32 * hb[3]) == Approx(0));
    CHECK(std::abs(ta(0, 0, 2, 0) - factor * Y20 * hb[2]) == Approx(0));
    CHECK(std::abs(ta(0, 0, 4, 0) - factor * Y40 * hb[4]) == Approx(0));
  }

  SECTION("Check diagonal recurrence") {
    auto const a11_31 = ta(1, 1, 3, 1);
    auto const a00_20 = ta(0, 0, 2, 0);
    auto const a00_40 = ta(0, 0, 4, 0);
    auto const left1 = std::sqrt(2e0 / 3e0) * a11_31;
    auto const right1 = std::sqrt(12e0 / 35e0) * a00_20 + std::sqrt(12e0 / 63e0) * a00_40;
    CHECK(left1.real() == Approx(right1.real()));
    CHECK(left1.imag() == Approx(right1.imag()));

    auto const a11_51 = ta(1, 1, 5, 1);
    auto const a00_60 = ta(0, 0, 6, 0);
    auto const left2 = std::sqrt(2e0 / 3e0) * a11_51;
    auto const right2 = std::sqrt(30e0 / 99e0) * a00_40 + std::sqrt(30e0 / 143e0) * a00_60;
    CHECK(left2.real() == Approx(right2.real()));
    CHECK(left2.imag() == Approx(right2.imag()));

    auto const a22_42 = ta(2, 2, 4, 2);
    auto const left3 = std::sqrt(12e0 / 15e0) * a22_42;
    auto const right3 = std::sqrt(30e0 / 63e0) * a11_31 + std::sqrt(12e0 / 99e0) * a11_51;
    CHECK(left3.real() == Approx(right3.real()));
    CHECK(left3.imag() == Approx(right3.imag()));

    auto const a11_42 = ta(1, 1, 4, 2);
    auto const a00_31 = ta(0, 0, 3, 1);
    auto const a00_51 = ta(0, 0, 5, 1);
    auto const left0 = std::sqrt(2e0 / 3e0) * a11_42;
    auto const right0 = std::sqrt(30e0 / 63e0) * a00_31 + std::sqrt(12e0 / 99e0) * a00_51;
    CHECK(left0.real() == Approx(right0.real()));
    CHECK(left0.imag() == Approx(right0.imag()));
  }

  SECTION("Check off-diagonal recurrence") {
    auto const a10_31 = ta(1, 0, 3, 1);
    auto const a00_21 = ta(0, 0, 2, 1);
    auto const a00_41 = ta(0, 0, 4, 1);
    auto const left0 = -std::sqrt(1e0 / 3e0) * a10_31;
    auto const right0 = -std::sqrt(8e0 / 35e0) * a00_21 + std::sqrt(15e0 / 63e0) * a00_41;
    CHECK(left0.real() == Approx(right0.real()));
    CHECK(left0.imag() == Approx(right0.imag()));

    auto const a52_31 = ta(5, 2, 3, 1);
    auto const a32_31 = ta(3, 2, 3, 1);
    auto const a42_21 = ta(4, 2, 2, 1);
    auto const a42_41 = ta(4, 2, 4, 1);
    auto const left1 = -std::sqrt(21e0 / 99e0) * a52_31;
    auto const right1 = -std::sqrt(12e0 / 63e0) * a32_31 + -std::sqrt(8e0 / 35e0) * a42_21 +
                        std::sqrt(15e0 / 63e0) * a42_41;
    CHECK(left1.real() == Approx(right1.real()));
    CHECK(left1.imag() == Approx(right1.imag()));
  }

  SECTION("Try and check pathologies") {
    auto const a52_00 = ta(5, 2, 0, 0);
    auto const a32_00 = ta(3, 2, 0, 0);
    auto const a42_10 = ta(4, 2, 1, 0);
    auto const left0 = -std::sqrt(21e0 / 99e0) * a52_00;
    auto const right0 = -std::sqrt(12e0 / 63e0) * a32_00 + std::sqrt(1e0 / 3e0) * a42_10;
    CHECK(left0.real() == Approx(right0.real()));
    CHECK(left0.imag() == Approx(right0.imag()));
  }
}

TEST_CASE("Translation-Addition positive m") {
  Spherical<t_real> const R(1e0, 0.42, 0.36);
  t_complex const waveK(1e0, 1.5e0);
  SECTION("Regular") { check_recurrence<details::CachedRecurrence>(R, waveK, true); }
  SECTION("Irregular") { check_recurrence<details::CachedRecurrence>(R, waveK, false); }
}

TEST_CASE("Translation-Addition all m") {
  Spherical<t_real> const R(1e0, 0.42, 0.36);
  t_complex const waveK(1e0, 1.5e0);
  SECTION("Negative regular") {
    TranslationAdditionCoefficients ta(R, waveK, true);
    TranslationAdditionCoefficients ta_conj(R, std::conj(waveK), true);
    CHECK(ta(3, -2, 5, -2).real() == Approx(ta_conj(3, 2, 5, 2).real()));
    CHECK(ta(3, -2, 5, -2).imag() == Approx(-ta_conj(3, 2, 5, 2).imag()));
    CHECK(ta(3, -2, 5, 2).real() == Approx(ta_conj(3, 2, 5, -2).real()));
    CHECK(ta(3, -2, 5, 2).imag() == Approx(-ta_conj(3, 2, 5, -2).imag()));
    CHECK(ta(3, -2, 5, -1).real() == Approx(-ta_conj(3, 2, 5, 1).real()));
    CHECK(ta(3, -2, 5, -1).imag() == Approx(ta_conj(3, 2, 5, 1).imag()));
    CHECK(ta(3, -2, 5, 1).real() == Approx(-ta_conj(3, 2, 5, -1).real()));
    CHECK(ta(3, -2, 5, 1).imag() == Approx(ta_conj(3, 2, 5, -1).imag()));
    CHECK(ta(5, -3, 3, 1).real() == Approx(ta_conj(5, 3, 3, -1).real()));
    CHECK(ta(5, -3, 3, 1).imag() == Approx(-ta_conj(5, 3, 3, -1).imag()));
    CHECK(ta(5, -3, 3, -1).real() == Approx(ta_conj(5, 3, 3, 1).real()));
    CHECK(ta(5, -3, 3, -1).imag() == Approx(-ta_conj(5, 3, 3, 1).imag()));
  }
  SECTION("Negative irregular") {
    TranslationAdditionCoefficients ta(R, waveK, false);
    TranslationAdditionCoefficients ta_conj(R, -std::conj(waveK), false);
    CHECK(ta(3, -2, 5, -2).real() == Approx(ta_conj(3, 2, 5, 2).real()));
    CHECK(ta(3, -2, 5, -2).imag() == Approx(-ta_conj(3, 2, 5, 2).imag()));
    CHECK(ta(3, -2, 5, 2).real() == Approx(ta_conj(3, 2, 5, -2).real()));
    CHECK(ta(3, -2, 5, 2).imag() == Approx(-ta_conj(3, 2, 5, -2).imag()));
    CHECK(ta(3, -2, 5, -1).real() == Approx(-ta_conj(3, 2, 5, 1).real()));
    CHECK(ta(3, -2, 5, -1).imag() == Approx(ta_conj(3, 2, 5, 1).imag()));
    CHECK(ta(3, -2, 5, 1).real() == Approx(-ta_conj(3, 2, 5, -1).real()));
    CHECK(ta(3, -2, 5, 1).imag() == Approx(ta_conj(3, 2, 5, -1).imag()));
    CHECK(ta(5, -3, 3, 1).real() == Approx(ta_conj(5, 3, 3, -1).real()));
    CHECK(ta(5, -3, 3, 1).imag() == Approx(-ta_conj(5, 3, 3, -1).imag()));
    CHECK(ta(5, -3, 3, -1).real() == Approx(ta_conj(5, 3, 3, 1).real()));
    CHECK(ta(5, -3, 3, -1).imag() == Approx(-ta_conj(5, 3, 3, 1).imag()));
  }
}

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
    Spherical<t_real> const R(1e0, 0.42, 0.36);
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
    Spherical<t_real> const R(1e0, 0.42, 0.36);
    t_complex const waveK(1e0, 1.5e0);
    CoAxialTranslationAdditionCoefficients tca(R, waveK, true);
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

void check_coaxial_translation(t_real expansion_pos, t_real reexpansion_pos, bool expansion_regular,
                               bool reexpansion_regular, t_int n, t_int m) {
  assert(!(expansion_regular and !reexpansion_regular));
  bool coeff_regular = expansion_regular == reexpansion_regular;
  Spherical<t_real> const translation(expansion_pos - reexpansion_pos, 0, 0);
  t_complex const waveK(1e0, 0);
  CoAxialTranslationAdditionCoefficients tca(translation, waveK, coeff_regular);
  auto const basis_func = expansion_regular ? optimet::bessel<Bessel> : optimet::bessel<Hankel1>;
  auto const re_basis_func =
      reexpansion_regular ? optimet::bessel<Bessel> : optimet::bessel<Hankel1>;
  t_complex translated = 0;
  for(t_int l = m; l < m + 20; l++) {
    translated += tca(n, m, l) * std::get<0>(re_basis_func(reexpansion_pos * waveK, l)).back() *
                  boost::math::spherical_harmonic(l, m, 0, 0);
  }
  auto expected = std::get<0>(basis_func(expansion_pos * waveK, n)).back() *
                  boost::math::spherical_harmonic(n, m, 0, 0);
  INFO("Testing translation for " << n << " m " << m << " regular " << expansion_regular << " to "
                                  << reexpansion_regular);
  CHECK(expected.real() == Approx(translated.real()));
  CHECK(expected.imag() == Approx(translated.imag()));
}

TEST_CASE("Coaxial translation") {
  std::uniform_real_distribution<> small_dist(0, 1);
  std::uniform_real_distribution<> large_dist(10, 100);
  for(t_int n = 0; n < 10; n++) {
    for(t_int m = -m; m <= n; m++) {
      // "Simple singular expanded in regular"
      t_real small = small_dist(*mersenne);
      t_real large = large_dist(*mersenne);
      t_real large_small_diff = large - small;
      check_coaxial_translation(large, small, false, true, n, 0);
      // "Simple singular expanded in singular"
      check_coaxial_translation(large, large_small_diff, false, false, n, 0);
      // "Simple regular expanded in regular"
      check_coaxial_translation(large, large_small_diff, true, true, n, 0);
      // Simple regular expanded in regular
      check_coaxial_translation(large, small, true, true, n, 0);
    }
  }
}
