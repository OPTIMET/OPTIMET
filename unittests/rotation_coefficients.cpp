#include "catch.hpp"

#include "RotationCoefficients.h"
#include "Types.h"
#include "constants.h"
#include <boost/math/special_functions/spherical_harmonic.hpp>
#include <iostream>
#include <memory>
#include <random>
#include <HarmonicsIterator.h>

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
  auto const first = static_cast<Real>(0.5) * bnm(n + 1, -mu - 1) * std::exp(Complex(0, phi)) *
                     (static_cast<Real>(1) - std::cos(theta)) *
                     recurrence(theta, phi, chi, n + 1, m - 1, mu + 1, cache);
  auto const second = static_cast<Real>(-0.5) * bnm(n + 1, mu - 1) * std::exp(-Complex(0, phi)) *
                      (static_cast<Real>(1) + std::cos(theta)) *
                      recurrence(theta, phi, chi, n + 1, m - 1, mu - 1, cache);
  auto const third =
      -anm(n, mu) * std::sin(theta) * recurrence(theta, phi, chi, n + 1, m - 1, mu, cache);
  auto const result = factor * (first + second + third);
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
  std::set<Triplet> triplets = {Triplet{1, 0, 0},    Triplet{1, 0, 1},   Triplet{1, 1, 0},
                                Triplet{1, 1, 1},    Triplet{1, 0, -1},  Triplet{1, -1, 0},
                                Triplet{1, -1, -1},  Triplet{2, 2, 2},   Triplet{2, -2, 2},
                                Triplet{10, 10, 2},  Triplet{10, -1, 2}, Triplet{15, 15, 2},
                                Triplet{10, 10, 10}, Triplet{15, 15, 5}};

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
    if(std::abs(expected) > tolerance)
      CHECK(std::abs(actual - expected) < tolerance);
    CHECK(std::abs(rot(n, m, mu) - std::conj(rot(n, -m, -mu))) < 1e-8);
  }
}

TEST_CASE("Distribution of rows (columns) across procs") {
  auto const world = mpi::Communicator();
  auto const nmax = 10;
  auto const nobjects = 25;
  auto const map = optimet::belos::map(nmax, nobjects, world);
  CHECK(map->isContiguous());
  CHECK(map->isDistributed() == (world.size() > 1));
  CHECK(map->isOneToOne());
  CHECK(map->getIndexBase() == 0);

  auto const min_local =
      (nobjects / world.size()) * world.rank() + std::min(nobjects % world.size(), world.rank());
  auto const local = (nobjects / world.size()) + (nobjects % world.size() > world.rank() ? 1 : 0);
  auto const all_harmonics = 2 * nmax * (nmax + 2);
  CHECK(map->getMinLocalIndex() == 0);
  CAPTURE(map->getMinLocalIndex());
  CAPTURE(map->getMaxLocalIndex());
  CAPTURE(map->getMinGlobalIndex());
  CAPTURE(map->getMaxGlobalIndex());
  CHECK(map->getNodeNumElements() == local * all_harmonics);
  CHECK(map->getMaxLocalIndex() == local * all_harmonics - 1);
  CHECK(map->getMinGlobalIndex() == min_local * all_harmonics);
  CHECK(map->getMaxGlobalIndex() == (min_local + local) * all_harmonics - 1);
}

TEST_CASE("Graph structure of the sparse rotation matrix") {
  auto const world = mpi::Communicator();
  auto const nmax = 10;
  auto const nobjects = 25;
  auto const graph = optimet::belos::crs_rotation_graph(nmax, nobjects, world);
  auto const sum_n = [](t_uint x) { return (x * x + x) / 2; };
  auto const sum_n2 = [](t_uint x) { return (2 * x * x * x + 3 * x * x + x) / 6; };
  auto const row_entries = nobjects * (4 * sum_n2(nmax) + 4 * sum_n(nmax) + 1);
  CHECK(graph.getGlobalNumEntries() == row_entries);
}
