#include "catch.hpp"

#include "RotationCoefficients.h"
#include "Types.h"
#include "constants.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>
#include <memory>
#include <random>

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
      auto const expected = ((-mu > 0 and -mu % 2 == 1) ? -1. : 1.) * factor *
                            boost::math::spherical_harmonic(n, -mu, theta, phi);
      auto const actual = rot(n, 0, mu);
      CHECK(std::abs(actual - expected) < 1e-8);
      CHECK(std::abs(rot(n, 0, -mu) - std::conj(rot(n, 0, mu))) < 1e-8);
    }
  }
}

typedef long double Real;
typedef std::complex<Real> Complex;
typedef std::tuple<t_uint, t_int, t_int> Triplet;
typedef std::map<Triplet, Complex> Cache;

Real bnm(t_uint n, t_int m) {
  if(std::abs(m) > n)
    return 0;
  auto const sgn = m >= 0 ? 1 : -1;
  return sgn * std::sqrt(static_cast<Real>((n - m - 1) * (n - m)) /
                         static_cast<Real>((2 * n - 1) * (2 * n + 1)));
}

Real anm(t_uint n, t_int m) {
  if(std::abs(m) > n)
    return 0;
  auto const absm = std::abs(m);
  return std::sqrt(static_cast<Real>((n + 1 + absm) * (n + 1 - absm)) /
                   static_cast<Real>((2 * n + 1) * (2 * n + 3)));
};

Complex recurrence(Real theta, Real phi, Real chi, t_uint n, t_int m, t_int mu, Cache &cache) {
  Real const pi = 3.14159265358979323846264338327950288419716939937510582097494459230781640628619;
  if(static_cast<t_uint>(std::abs(m)) > n or static_cast<t_uint>(std::abs(mu)) > n)
    return 0;
  if(m < 0)
    return std::conj(recurrence(theta, phi, chi, n, -m, -mu, cache));
  auto const i_found = cache.find(Triplet{n, m, mu});
  if(i_found != cache.end())
    return i_found->second;
  if(m == 0) {
    auto const gamma = -mu > 0 and -mu % 2 == 1 ?
                           -boost::math::spherical_harmonic(n, -mu, theta, phi) :
                           boost::math::spherical_harmonic(n, -mu, theta, phi);
    auto const result = std::sqrt(static_cast<Real>(4) * pi / static_cast<Real>(2 * n + 1)) * gamma;
    cache.emplace(Triplet{n, m, mu}, result);
    return result;
  }

  auto const factor = std::exp(Complex(0, chi)) / bnm(n + 1, m - 1);
  auto const first = factor * static_cast<Real>(0.5) * bnm(n + 1, -mu - 1) *
                     std::exp(Complex(0, phi)) * (static_cast<Real>(1) - std::cos(theta)) *
                     recurrence(theta, phi, chi, n + 1, m - 1, mu + 1, cache);
  auto const second = factor * static_cast<Real>(-0.5) * bnm(n + 1, mu - 1) *
                      std::exp(-Complex(0, phi)) * (static_cast<Real>(1) + std::cos(theta)) *
                      recurrence(theta, phi, chi, n + 1, m - 1, mu - 1, cache);
  auto const third =
      factor * -anm(n, mu) * std::sin(theta) * recurrence(theta, phi, chi, n + 1, m - 1, mu, cache);
  auto const result = first + second + third;
  cache.emplace(Triplet{n, m, mu}, result);
  return result;
};

TEST_CASE("Check recurrence") {
  std::uniform_real_distribution<> rdist(0, constant::pi);
  auto const theta = rdist(*mersenne);
  auto const phi = 2 * rdist(*mersenne);
  auto const chi = rdist(*mersenne);
  RotationCoefficients rot(theta, phi, chi);

  // Creates a set of n to look at
  Cache cache;
  std::uniform_int_distribution<> ndist(2, 50);
  std::set<Triplet> triplets = {
      Triplet{1, 0, 0},    Triplet{1, 0, 1},    Triplet{1, 1, 0},    Triplet{1, 1, 1},
      Triplet{1, 0, -1},   Triplet{1, -1, 0},   Triplet{1, -1, -1},  Triplet{2, 2, 2},
      Triplet{2, -2, 2},   Triplet{10, 10, 2},  Triplet{10, -1, 2},  Triplet{15, 15, 2},
      Triplet{20, 20, 20}, Triplet{21, 21, 21}, Triplet{22, 22, 22}, Triplet{23, 23, 23},
      Triplet{24, 24, 24}, Triplet{25, 25, 25}, Triplet{26, 26, 26}, Triplet{27, 27, 27},
      Triplet{28, 28, 28}, Triplet{29, 29, 29}, Triplet{30, 30, 30},
  };

  for(auto const triplet : triplets) {
    auto const n = std::get<0>(triplet);
    auto const m = std::get<1>(triplet);
    auto const mu = std::get<2>(triplet);
    Complex const expected = recurrence(theta, phi, chi, n, m, mu, cache);
    Complex const actual = rot(n, m, mu);
    INFO("n: " << n << ", m: " << m << ", mu: " << mu << ", actual: " << actual
               << ", expected: " << expected);
    auto const tolerance =
        std::max(static_cast<Real>(1e-12), std::abs(expected) * static_cast<Real>(1e-8));
    CHECK(std::abs(actual - expected) < tolerance);
    CHECK(std::abs(rot(n, m, mu) - std::conj(rot(n, -m, -mu))) < 1e-8);
  }
}

TEST_CASE("Coefficient sum") {
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

// TEST_CASE("All coefficients") {
//   std::uniform_real_distribution<> rdist(0, constant::pi);
//   auto const theta = rdist(*mersenne);
//   auto const phi = 2 * rdist(*mersenne);
//   auto const chi = rdist(*mersenne);
//   RotationCoefficients rot(theta, phi, chi);
//
//   auto const coeffs = rot.all_coefficients(10, 10, 0);
//   std::cout << std::setw(18) << std::setprecision(6) << std::setfill(' ') << std::right
//             << std::scientific << "\n";
//   for(auto const c : coeffs)
//     if(std::get<1>(c.first) >= 0)
//       std::cout << std::get<0>(c.first) << ", " << std::get<1>(c.first) << ", "
//                 << std::get<2>(c.first) << ": " << c.second.transpose() << std::endl;
// }
