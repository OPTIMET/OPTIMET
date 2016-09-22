#include "catch.hpp"

#include "Bessel.h"
#include "CoAxialTranslationCoefficients.h"
#include "constants.h"
#include <Coefficients.h>
#include <boost/math/special_functions/bessel.hpp>
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <cmath>
#include <iostream>
#include <random>

extern std::unique_ptr<std::mt19937_64> mersenne;

using namespace optimet;

void check_coaxial_n_recurrence(CachedCoAxialRecurrence tca, t_int n, t_int m, t_int l) {
  // This aims to check that we correctly implement formula 4.79
  using coefficient::a;
  auto left0 = a(n - 1, m) * static_cast<t_complex>(tca(n - 1, m, l)) -
               a(n, m) * static_cast<t_complex>(tca(n + 1, m, l));
  auto right0 = a(l, m) * static_cast<t_complex>(tca(n, m, l + 1)) -
                a(l - 1, m) * static_cast<t_complex>(tca(n, m, l - 1));
  INFO("Testing n recurrence for n " << n << " m " << m << " l " << l);
  CHECK(left0.real() == Approx(right0.real()));
  CHECK(left0.imag() == Approx(right0.imag()));
}

void check_coaxial_m_recurrence(CachedCoAxialRecurrence tca, t_int n, t_int m, t_int l) {
  // This aims to check that we correctly implement formula 4.80
  using coefficient::b;
  auto left0 = b(n, m) * static_cast<t_complex>(tca(n - 1, m + 1, l)) -
               b(n + 1, -m - 1) * static_cast<t_complex>(tca(n + 1, m + 1, l));
  auto right0 = b(l + 1, m) * static_cast<t_complex>(tca(n, m, l + 1)) -
                b(l, -m - 1) * static_cast<t_complex>(tca(n, m, l - 1));
  INFO("Testing m recurrence for n " << n << " m " << m << " l " << l);
  CHECK(left0.real() == Approx(right0.real()));
  CHECK(left0.imag() == Approx(right0.imag()));
}

void check_coaxial_mn_recurrence(CachedCoAxialRecurrence tca, t_int n, t_int m, t_int l) {
  // This aims to check that we correctly implement formula 4.84
  using coefficient::b;
  auto left0 = b(n + 1, -m - 1) * static_cast<t_complex>(tca(n + 1, m + 1, l));
  auto right0 = b(l, -m - 1) * static_cast<t_complex>(tca(n, m, l - 1)) -
                b(l + 1, m) * static_cast<t_complex>(tca(n, m, l + 1));
  INFO("Testing m=n recurrence for n " << n << " m " << m << " l " << l);
  CHECK(left0.real() == Approx(right0.real()));
  CHECK(left0.imag() == Approx(right0.imag()));
}

void check_coaxial_m_symmetry(CachedCoAxialRecurrence tca, t_int n, t_int m, t_int l) {
  // check that value is independent of m
  auto left0 = tca(n, m, l);
  auto right0 = tca(n, -m, l);
  INFO("Testing m symmetry for " << n << " m " << m << " l " << l);
  CHECK(left0.real() == Approx(right0.real()));
  CHECK(left0.imag() == Approx(right0.imag()));
}

void check_coaxial_ln_symmetry(CachedCoAxialRecurrence tca, t_int n, t_int m, t_int l) {
  // check that value has the right sign change for m<->l
  auto left0 = static_cast<t_complex>(tca(n, m, l));
  t_complex sign = ((n + l) % 2 == 0 ? 1 : -1);
  auto right0 = static_cast<t_complex>(tca(l, m, n)) * sign;
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
    CachedCoAxialRecurrence tca(R, waveK, true);
    // Numbers are generated from the same formula in Scipy
    CHECK(static_cast<t_complex>(tca(0, 0, 0)).real() == Approx(1.1400511799225792));
    CHECK(static_cast<t_complex>(tca(0, 0, 0)).imag() == Approx(-0.55962217045848206));
    CHECK(static_cast<t_complex>(tca(0, 0, 4)).real() == Approx(-0.028191522402192234));
    CHECK(static_cast<t_complex>(tca(0, 0, 4)).imag() == Approx(-0.02162885905593049));
    CHECK(static_cast<t_complex>(tca(1, 0, 1)).real() == Approx(1.2274819687880665));
    CHECK(static_cast<t_complex>(tca(1, 0, 1)).imag() == Approx(-1.0271756758800463));
    CHECK(static_cast<t_complex>(tca(-1, 1, 3)).real() == Approx(0));
    CHECK(static_cast<t_complex>(tca(-1, 1, 3)).imag() == Approx(0));
    CHECK(static_cast<t_complex>(tca(1, 0, -1)).real() == Approx(0));
    CHECK(static_cast<t_complex>(tca(1, 0, -1)).imag() == Approx(0));
    CHECK(static_cast<t_complex>(tca(1, 1, 3)).real() == Approx(-0.085169586217943016));
    CHECK(static_cast<t_complex>(tca(1, 1, 3)).imag() == Approx(0.36331568009355053));

    t_complex expected;
    expected = (static_cast<t_complex>(tca(1, 0, 1)) * a(1, 0) +
                static_cast<t_complex>(tca(0, 0, 2)) * a(0, 0) -
                static_cast<t_complex>(tca(1, 0, 3)) * a(2, 0)) /
               (a(1, 0));
    CHECK(static_cast<t_complex>(tca(2, 0, 2)).real() == Approx(expected.real()));
    CHECK(static_cast<t_complex>(tca(2, 0, 2)).imag() == Approx(expected.imag()));

    expected = (static_cast<t_complex>(tca(0, 0, 2)) * a(2, 0) -
                static_cast<t_complex>(tca(0, 0, 4)) * a(3, 0)) /
               (a(0, 0));
    CHECK(static_cast<t_complex>(tca(1, 0, 3)).real() == Approx(expected.real()));
    CHECK(static_cast<t_complex>(tca(1, 0, 3)).imag() == Approx(expected.imag()));

    expected = (static_cast<t_complex>(tca(1, 0, 3)) * a(3, 0) +
                static_cast<t_complex>(tca(0, 0, 4)) * a(0, 0) -
                static_cast<t_complex>(tca(1, 0, 5)) * a(4, 0)) /
               (a(1, 0));
    CHECK(static_cast<t_complex>(tca(2, 0, 4)).real() == Approx(expected.real()));
    CHECK(static_cast<t_complex>(tca(2, 0, 4)).imag() == Approx(expected.imag()));

    expected = (-static_cast<t_complex>(tca(0, 0, 2)) * b(3, -1) +
                static_cast<t_complex>(tca(0, 0, 4)) * b(4, 0)) /
               (-b(1, -1));
    CHECK(static_cast<t_complex>(tca(1, 1, 3)).real() == Approx(expected.real()));
    CHECK(static_cast<t_complex>(tca(1, 1, 3)).imag() == Approx(expected.imag()));

    expected = (-static_cast<t_complex>(tca(2, 0, 2)) * b(3, -1) -
                static_cast<t_complex>(tca(1, 1, 3)) * b(2, 0) +
                static_cast<t_complex>(tca(2, 0, 4)) * b(4, 0)) /
               (-b(3, -1));
    CHECK(static_cast<t_complex>(tca(3, 1, 3)).real() == Approx(expected.real()));
    CHECK(static_cast<t_complex>(tca(3, 1, 3)).imag() == Approx(expected.imag()));
  }
  SECTION("R zero") {
    using coefficient::a;
    using coefficient::b;
    t_real const R(0.0);
    t_complex const waveK(1e0, 0.0);
    CachedCoAxialRecurrence tca(R, waveK, true);
    // Check that the simple 0,0 coeff is as expected. This revealed
    // an issue in the bessel function implementation for (0,0)
    CHECK(static_cast<t_complex>(tca(0, 0, 0)).real() == Approx(1.0));
    CHECK(static_cast<t_complex>(tca(0, 0, 0)).imag() == Approx(0.0));
    CHECK(std::get<0>(optimet::bessel<Bessel>(0, 0))[0].real() == Approx(1.0));
  }
  SECTION("Check n and m recurrence") {
    t_complex const waveK(1e0, 1.5e0);
    CachedCoAxialRecurrence tca(1.0, waveK, true);
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

void check_coaxial_translation_zero(t_int n, t_int m, t_complex waveK, bool regular) {
  assert(m <= n);
  bool coeffs_regular;
  std::uniform_real_distribution<> theta_dist(0, consPi / 2);
  std::uniform_real_distribution<> phi_dist(0, 2 * consPi);
  std::uniform_real_distribution<> r_p_dist(1, 10);
  t_real r_p = r_p_dist(*mersenne);
  t_real theta_p = theta_dist(*mersenne);
  t_real phi = phi_dist(*mersenne);

  t_real r_pq = 0.0;

  t_real theta_q = theta_p;
  t_real r_q = r_p;

  auto basis_func = optimet::bessel<Bessel>;
  auto re_basis_func = optimet::bessel<Bessel>;
  coeffs_regular = true;
  if(regular == false) {
    basis_func = optimet::bessel<Hankel1>;
    re_basis_func = optimet::bessel<Hankel1>;
  }

  CachedCoAxialRecurrence tca(r_pq, waveK, coeffs_regular);
  t_complex translated = 0;
  for(t_int l = std::abs(m); l < std::abs(m) + 105; l++) {
    translated += static_cast<t_complex>(tca(n, m, l)) *
                  std::get<0>(re_basis_func(r_q * waveK, l)).back() *
                  boost::math::spherical_harmonic(l, m, theta_q, phi);
  }
  auto expected = std::get<0>(basis_func(r_p * waveK, n)).back() *
                  boost::math::spherical_harmonic(n, m, theta_p, phi);
  INFO("Testing translation for n: " << n << " m: " << m << " sph "
                                     << boost::math::spherical_harmonic(n, m, theta_p, phi));
  CHECK(expected.real() == Approx(translated.real()));
  CHECK(expected.imag() == Approx(translated.imag()));
}

void check_coaxial_translation_onaxis(t_real expansion_pos, t_real reexpansion_pos,
                                      bool expansion_regular, bool reexpansion_regular, t_int n,
                                      t_int m, t_complex waveK) {
  assert(!(expansion_regular and !reexpansion_regular));
  assert(m <= n);
  bool coeff_regular = expansion_regular == reexpansion_regular;
  t_real const translation(expansion_pos - reexpansion_pos);
  CachedCoAxialRecurrence tca(translation, waveK, coeff_regular);
  auto const basis_func = expansion_regular ? optimet::bessel<Bessel> : optimet::bessel<Hankel1>;
  auto const re_basis_func =
      reexpansion_regular ? optimet::bessel<Bessel> : optimet::bessel<Hankel1>;
  t_complex translated = 0;
  for(t_int l = std::abs(m); l < std::abs(m) + 25; l++) {
    translated += static_cast<t_complex>(tca(n, m, l)) *
                  std::get<0>(re_basis_func(reexpansion_pos * waveK, l)).back() *
                  boost::math::spherical_harmonic(l, m, 0, 0);
  }
  auto expected = std::get<0>(basis_func(expansion_pos * waveK, n)).back() *
                  boost::math::spherical_harmonic(n, m, 0, 0);
  INFO("Testing translation for n: " << n << " m: " << m << " regular: " << expansion_regular
                                     << " to: " << reexpansion_regular);
  CHECK(expected.real() == Approx(translated.real()));
  CHECK(expected.imag() == Approx(translated.imag()));
}

void check_coaxial_translation_off_axis_reexpand_fixed_r(t_int n, t_int m, t_complex waveK,
                                                         bool regular, bool reexpansion_regular) {
  assert(m <= n);
  bool coeffs_regular;
  auto basis_func = optimet::bessel<Bessel>;
  auto re_basis_func = optimet::bessel<Bessel>;
  t_real theta_p;
  t_real r_p;
  t_real r_pq;
  t_real phi;
  t_real theta_q;
  t_real r_q;
  if(regular and reexpansion_regular) {
    theta_p = 0.0897864;
    r_p = 5.88056;
    r_pq = 2.4808;
    phi = 1.2131;
    theta_q = 0.154931;
    r_q = 3.41701;
    basis_func = optimet::bessel<Bessel>;
    re_basis_func = optimet::bessel<Bessel>;
    coeffs_regular = true;
  } else if(!regular and reexpansion_regular) {
    theta_p = 0.328313;
    r_p = 5.61927;
    r_pq = 4.25386;
    phi = 1.2131;
    theta_q = 1.0393;
    r_q = 2.10187;
    basis_func = optimet::bessel<Hankel1>;
    re_basis_func = optimet::bessel<Bessel>;
    coeffs_regular = false;
    // Convergence is problematic when r_q ~= r_pq such as below
    // t_real r_p = 5;
    // t_real r_q = 5;
    // t_real r_pq = 5;
    // t_real theta_p = consPi / 3.0;
    // t_real theta_q = 2 * consPi / 3.0;
    // t_real phi = 1.2131;
  } else if(!regular and !reexpansion_regular) {
    theta_p = 0.0897864;
    r_p = 5.88056;
    r_pq = 2.4808;
    phi = 1.2131;
    theta_q = 0.154931;
    r_q = 3.41701;
    basis_func = optimet::bessel<Hankel1>;
    re_basis_func = optimet::bessel<Hankel1>;
    coeffs_regular = true;
    // Convergence is problematic when r_q ~= r_pq such as below
    // t_real r_p = 5;
    // t_real r_q = 5;
    // t_real r_pq = 5;
    // t_real theta_p = consPi / 3.0;
    // t_real theta_q = 2 * consPi / 3.0;
    // t_real phi = 1.2131;
  } else {
    std::cout << "this should never happen" << std::endl;
  }
  CHECK(sin(theta_p) * r_p == Approx(sin(theta_q) * r_q));

  CachedCoAxialRecurrence tca(r_pq, waveK, coeffs_regular);
  t_complex translated = 0;
  for(t_int l = std::abs(m); l < std::abs(m) + 105; l++) {
    translated += static_cast<t_complex>(tca(n, m, l)) *
                  std::get<0>(re_basis_func(r_q * waveK, l)).back() *
                  boost::math::spherical_harmonic(l, m, theta_q, phi);
  }
  auto expected = std::get<0>(basis_func(r_p * waveK, n)).back() *
                  boost::math::spherical_harmonic(n, m, theta_p, phi);
  INFO("Testing translation for regular? "
       << regular << " reexpanded in regular? " << reexpansion_regular << " n: " << n << " m: " << m
       << " sph " << boost::math::spherical_harmonic(n, m, theta_p, phi) << " bessel "
       << std::get<0>(basis_func(r_p * waveK, n)).back() << r_p * waveK << " r_p " << r_p << " r_q "
       << r_q << " r_pq " << r_pq);
  CHECK(expected.real() == Approx(translated.real()));
  CHECK(expected.imag() == Approx(translated.imag()));
}

void check_coaxial_translation_off_axis_reexpand_iregular(t_int n, t_int m, t_complex waveK) {
  assert(m <= n);
  //! Inner floating point with higher precision
  typedef long double Real;
  //! Inner complex floating point with higher precision
  typedef std::complex<Real> Complex;
  bool coeffs_regular;
  std::uniform_real_distribution<> theta_dist(0, consPi / 2);
  std::uniform_real_distribution<> phi_dist(0, 2 * consPi);
  std::uniform_real_distribution<> r_p_dist(1, 10);
  t_real r_p = r_p_dist(*mersenne);
  t_real theta_p = theta_dist(*mersenne);
  t_real phi = phi_dist(*mersenne);
  t_real theta_p1 = theta_p;
  if(theta_p > consPi / 2)
    // internal angle in the triangle < pi/2
    theta_p1 = consPi - theta_p;

  t_real z_p = cos(theta_p1) * r_p;
  std::uniform_real_distribution<> r_pq_dist(-10, 10);
  t_real r_pq = 1;
  t_real r_q = 1;
  t_real rho_p;
  t_real rho_q;
  t_real z_q;
  t_real theta_q;
  t_int n1 = 0;
  r_pq = r_pq_dist(*mersenne);
  rho_p = sin(theta_p1) * r_p;
  rho_q = rho_p;
  z_q = 0;

  if((theta_p <= consPi / 2))
    z_q = std::abs(z_p - r_pq);
  else
    z_q = std::abs(z_p + r_pq);
  theta_q = atan(rho_q / z_q);
  r_q = z_q / cos(theta_q);
  n1 = n1 + 1;

  CHECK(sin(theta_p1) * r_p == Approx(sin(theta_q) * r_q));
  if((theta_p <= consPi / 2) and (r_pq > z_p))
    // swap to get the angle from positive z axis
    theta_q = consPi - theta_q;
  if((theta_p > consPi / 2) and (r_pq > -z_p))
    theta_q = consPi - theta_q;
  auto const basis_func = optimet::bessel<Hankel1>;
  auto re_basis_func = optimet::bessel<Bessel>;
  if(std::abs(r_q) <= std::abs(r_pq)) {
    coeffs_regular = false;
    re_basis_func = optimet::bessel<Bessel>;
  } else {
    coeffs_regular = true;
    re_basis_func = optimet::bessel<Hankel1>;
  }
  CachedCoAxialRecurrence tca(r_pq, waveK, coeffs_regular);
  Complex translated = 0;
  for(t_int l = std::abs(m); l < std::abs(m) + 105; l++) {
    translated += static_cast<t_complex>(tca(n, m, l)) *
                  std::get<0>(re_basis_func(r_q * waveK, l)).back() *
                  boost::math::spherical_harmonic(l, m, theta_q, phi);
  }
  auto expected = std::get<0>(basis_func(r_p * waveK, n)).back() *
                  boost::math::spherical_harmonic(n, m, theta_p, phi);
  INFO("Testing translation for regular ? 0 "
       << " reexpanded in regular? " << !coeffs_regular << " n: " << n << " m: " << m << " sph "
       << boost::math::spherical_harmonic(n, m, theta_p, phi) << " bessel "
       << std::get<0>(basis_func(r_p * waveK, n)).back() << r_p * waveK << " r_p " << r_p << " r_q "
       << r_q << " r_pq " << r_pq);
  CHECK(expected.real() == Approx(translated.real()));
  CHECK(expected.imag() == Approx(translated.imag()));
}

void check_coaxial_translation_off_axis_reexpand_regular(t_int n, t_int m, t_complex waveK) {
  assert(m <= n);
  //! Inner floating point with higher precision
  typedef long double Real;
  //! Inner complex floating point with higher precision
  typedef std::complex<Real> Complex;
  bool coeffs_regular;
  std::uniform_real_distribution<> theta_dist(0, consPi);
  std::uniform_real_distribution<> phi_dist(0, 2 * consPi);
  std::uniform_real_distribution<> r_p_dist(1, 10);
  t_real r_p = r_p_dist(*mersenne);
  t_real theta_p = theta_dist(*mersenne);
  t_real phi = phi_dist(*mersenne);
  t_real theta_p1 = theta_p;
  if(theta_p > consPi / 2)
    // internal angle in the triangle < pi/2
    theta_p1 = consPi - theta_p;

  t_real z_p = cos(theta_p1) * r_p;
  std::uniform_real_distribution<> r_pq_dist(-10, 10);
  t_real r_pq = r_pq_dist(*mersenne);
  t_real rho_p = sin(theta_p1) * r_p;
  t_real rho_q = rho_p;
  t_real z_q = 0;

  if((theta_p <= consPi / 2))
    z_q = std::abs(z_p - r_pq);
  else
    z_q = std::abs(z_p + r_pq);
  t_real theta_q = atan(rho_q / z_q);
  t_real r_q = z_q / cos(theta_q);

  CHECK(sin(theta_p1) * r_p == Approx(sin(theta_q) * r_q));
  if((theta_p <= consPi / 2) and (r_pq > z_p))
    // swap to get the angle from positive z axis
    theta_q = consPi - theta_q;
  if((theta_p > consPi / 2) and (r_pq > -z_p))
    theta_q = consPi - theta_q;
  auto const basis_func = optimet::bessel<Bessel>;
  auto const re_basis_func = optimet::bessel<Bessel>;
  coeffs_regular = true;

  CachedCoAxialRecurrence tca(r_pq, waveK, coeffs_regular);
  t_complex translated = 0;
  for(t_int l = std::abs(m); l < std::abs(m) + 105; l++) {
    translated += static_cast<t_complex>(tca(n, m, l)) *
                  std::get<0>(re_basis_func(r_q * waveK, l)).back() *
                  boost::math::spherical_harmonic(l, m, theta_q, phi);
  }
  auto expected = std::get<0>(basis_func(r_p * waveK, n)).back() *
                  boost::math::spherical_harmonic(n, m, theta_p, phi);
  INFO("Testing translation for regular 1 reexpanded in regular n: "
       << n << " m: " << m << " sph " << boost::math::spherical_harmonic(n, m, theta_p, phi)
       << " bessel " << std::get<0>(basis_func(r_p * waveK, n)).back() << r_p * waveK << " r_p "
       << r_p << " r_q " << r_q << " r_pq " << r_pq);
}

TEST_CASE("Coaxial translation") {
  std::uniform_real_distribution<> small_dist(0, 1);
  std::uniform_real_distribution<> large_dist(10, 50);
  std::uniform_real_distribution<> wave_dist(0.1, 1);
  t_real waver = wave_dist(*mersenne);
  t_real wavei = wave_dist(*mersenne);
  t_complex waveK(waver, wavei);
  for(t_int n = 0; n < 5; n++) {
    for(t_int m = -n; m <= n; m++) {
      // "Simple singular expanded in regular"
      t_real small = small_dist(*mersenne);
      t_real large = large_dist(*mersenne);
      t_real large_small_diff = large - small;

      check_coaxial_translation_zero(n, m, waveK, false);

      check_coaxial_translation_zero(n, m, waveK, true);

      check_coaxial_translation_onaxis(large, small, false, true, n, m, waveK);
      // "Simple singular expanded in singular"
      check_coaxial_translation_onaxis(large, large_small_diff, false, false, n, m, waveK);
      // "Simple regular expanded in regular"
      check_coaxial_translation_onaxis(large, large_small_diff, true, true, n, m, waveK);
      // "Simple regular expanded in regular"
      check_coaxial_translation_onaxis(large, small, true, true, n, m, waveK);
      // Zero translation
      check_coaxial_translation_onaxis(large, large, true, true, n, m, waveK);

      check_coaxial_translation_onaxis(large, large, false, false, n, m, waveK);

      check_coaxial_translation_off_axis_reexpand_fixed_r(n, m, waveK, true, true);

      check_coaxial_translation_off_axis_reexpand_fixed_r(n, m, waveK, false, true);

      check_coaxial_translation_off_axis_reexpand_fixed_r(n, m, waveK, false, false);

      check_coaxial_translation_off_axis_reexpand_regular(n, m, waveK);

      check_coaxial_translation_off_axis_reexpand_iregular(n, m, waveK);
    }
  }
}

Vector<t_real> to_spherical(Vector<t_real> const &x) {
  assert(x.size() == 3);
  auto const r = x.stableNorm();
  auto const theta = std::atan2(x[1], x[0]);
  auto const phi = r > 1e-12 ? std::acos(x[2] / r) : 0e0;
  Vector<t_real> result(3);
  result << r, phi, theta > 0 ? theta : theta + 2 * constant::pi;
  return result;
};

t_complex
basis_function(t_complex waveK, bool expansion, Vector<t_real> const &r, t_int n, t_int m) {
  auto const function = expansion ? optimet::bessel<Bessel> : optimet::bessel<Hankel1>;
  auto const spherical = to_spherical(r);
  return std::get<0>(function(spherical(0) * waveK, n)).back() *
         boost::math::spherical_harmonic(n, m, spherical(1), spherical(2));
};

std::function<t_complex(Vector<t_real> const &r)>
field(bool expansion, t_complex waveK, Vector<t_complex> const &coeffs) {

  t_int const nmax = std::lround(std::sqrt(coeffs.rows() + 1) - 1.0);
  assert(nmax * (nmax + 2) == coeffs.rows());
  assert(nmax > 0);

  return [nmax, waveK, expansion, coeffs](Vector<t_real> const &r) {
    t_complex result = 0;
    for(auto n = 1; n <= nmax; ++n)
      for(auto m = -n; m <= n; ++m) {
        auto const index = std::abs(m) > n ? 0 : n * n - 1 + n + m;
        result += coeffs(index) * basis_function(waveK, expansion, r, n, m);
      }
    return result;
  };
}

TEST_CASE("Matrix interface") {
  auto const nmax = 1;
  bool const in_regular = true;
  bool const out_regular = true;
  t_complex const waveK(1e0, 1.5e0);
  Vector<t_real> const r_p = Vector<t_real>::Zero(3);
  Vector<t_real> r_q = (Vector<t_real>::Random(3) * 10).array() - 5e0;
  r_q(0) = 0;
  r_q(1) = 0;
  r_q(2) = 0;
  CoAxialTranslationAdditionCoefficients tca((r_q - r_p)(2), waveK,
                                             in_regular == out_regular);

  Matrix<t_complex> const coeffs = Vector<t_complex>::Random(nmax * (nmax + 2), 1);
  auto const input_field = field(in_regular, waveK, coeffs);
  auto const output_field = field(out_regular, waveK, tca(coeffs));

  Vector<t_real> x = (Vector<t_real>::Random(3) * 10).array() - 5e0;
  x(0) = 0; x(1) = 0;
  CHECK(input_field(x - r_p).real() == Approx(output_field(x - r_q).real()));
  CHECK(input_field(x - r_p).imag() == Approx(output_field(x - r_q).imag()));
}
