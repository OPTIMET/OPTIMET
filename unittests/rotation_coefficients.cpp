#include "catch.hpp"

#include "RotationCoefficients.h"
#include "Types.h"
#include "constants.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <memory>
#include <random>
#include <iostream>

extern std::unique_ptr<std::mt19937_64> mersenne;
using namespace optimet;

TEST_CASE("Check Initials") {
  std::uniform_real_distribution<> rdist(0, constant::pi);
  auto const theta = rdist(*mersenne);
  auto const phi = 2 * rdist(*mersenne);
  auto const chi = rdist(*mersenne);
  RotationCoefficients rot(theta, phi, chi);

  CHECK(std::abs(rot(0, 1, -1)) < 1e-12);
  for(t_uint n = 5; n > 0; --n) {
    auto const factor = std::sqrt(4 * constant::pi / static_cast<t_real>(2 * n + 1));
    for(t_int mu = -static_cast<t_int>(n); mu <= static_cast<t_int>(n); ++mu) {
      auto const expected =
          factor * boost::math::spherical_harmonic(n, mu, theta, phi) * std::pow(-1, -mu);
      auto const actual = rot(n, 0, -mu);
      CAPTURE(expected);
      CAPTURE(actual);
      CHECK(std::abs(actual - expected) < 1e-8);
    }
  }
}

TEST_CASE("Check recurrence") {
  std::uniform_real_distribution<> rdist(0, constant::pi);
  auto const theta = rdist(*mersenne);
  auto const phi = 2 * rdist(*mersenne);
  auto const chi = rdist(*mersenne);
  RotationCoefficients rot(theta, phi, chi);

  auto const bnm = [](t_uint n, t_int m) -> t_real {
    if(std::abs(m) > n)
      return 0;
    auto const sgn = m >= 0 ? 1 : -1;
    return sgn * std::sqrt(static_cast<t_real>((n - m - 1) * (n - m)) /
                           static_cast<t_real>((2 * n - 1) * (2 * n + 1)));
  };
  auto const anm = [](t_uint n, t_int m) -> t_real {
    if(std::abs(m) > n)
      return 0;
    auto const absm = std::abs(m);
    return std::sqrt(static_cast<t_real>((n + 1 + absm) * (n + 1 - absm)) /
                     static_cast<t_real>((2 * n + 1) * (2 * n + 3)));
  };

  // The right hand-side of the recursion formula itself
  auto const recurrence = [&rot, theta, phi, chi, &bnm, &anm](t_uint n, t_int m,
                                                              t_int mu) -> t_complex {
    auto const factor = std::exp(t_complex(0, chi)) / bnm(n + 1, m - 1);
    auto const first = 0.5 * bnm(n + 1, -mu - 1) * std::exp(t_complex(0, phi)) *
                       (1 - std::cos(theta)) * rot(n + 1, m - 1, mu + 1);
    auto const second = -0.5 * bnm(n + 1, mu - 1) * std::exp(-t_complex(0, phi)) *
                        (1 + std::cos(theta)) * rot(n + 1, m - 1, mu - 1);
    auto const third = -anm(n, mu) * std::sin(theta) * rot(n + 1, m - 1, mu);
    return factor * (first + second + third);
  };

  // Creates a set of n to look at
  std::uniform_int_distribution<> ndist(2, 50);
  std::set<t_uint> ns = {30};
  // for(auto i = 0; i < 10; ++i)
  //   ns.insert(ndist(*mersenne));

  for(auto const n : ns)
    for(t_int m(1); m <= static_cast<t_int>(n); ++m)
      for(t_int mu(-static_cast<t_int>(n)); mu <= static_cast<t_int>(n); ++mu) {
        std::cout << n << " " << m << " " << mu << ": " << rot(n, m, mu) << std::endl;
        REQUIRE(std::abs(rot(n, m, mu) - recurrence(n, m, mu)) < 1e-8);
        CHECK(std::abs(rot(n, m, mu) - std::conj(rot(n, -m, -mu))) < 1e-8);
      }
}

TEST_CASE("All coeffs") {
  std::uniform_real_distribution<> rdist(0, constant::pi);
  auto const theta = rdist(*mersenne);
  auto const phi = 2 * rdist(*mersenne);
  auto const chi = rdist(*mersenne);
  RotationCoefficients rot(theta, phi, chi);

  auto const through_coeffs = [&rot](t_uint n, t_int m, t_int mu) {
    auto const coeffs = rot.coefficients(n, m, mu);
    return std::accumulate(coeffs.begin(), coeffs.end(), t_complex(0));
  };
  CHECK(std::abs(rot(1, 1, 0) - through_coeffs(1, 1, 0)) < 1e-8);
  CHECK(std::abs(rot(2, 1, 0) - through_coeffs(2, 1, 0)) < 1e-8);
  CHECK(std::abs(rot(4, 4, 0) - through_coeffs(4, 4, 0)) < 1e-8);
  CHECK(std::abs(rot(4, -4, 0) - through_coeffs(4, -4, 0)) < 1e-8);
}
