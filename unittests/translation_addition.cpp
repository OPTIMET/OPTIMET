// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#include "catch.hpp"

#include "Bessel.h"
#include "TranslationAdditionCoefficients.h"
#include "constants.h"
#include <boost/math/special_functions/legendre.hpp>
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
