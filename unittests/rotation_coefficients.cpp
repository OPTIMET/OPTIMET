#include "catch.hpp"

#include "RotationCoefficients.h"
#include "Types.h"
#include "constants.h"
#include <Eigen/Dense>
#include <Eigen/LU>
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
  RotationCoefficients rotation(theta, phi, chi);

  CHECK(std::abs(rotation(0, 1, -1)) < 1e-12);
  for(t_uint n = 5; n > 0; --n) {
    auto const factor = std::sqrt(4 * constant::pi / static_cast<t_real>(2 * n + 1));
    for(t_int mu = -static_cast<t_int>(n); mu <= static_cast<t_int>(n); ++mu) {
      auto const expected = ((-mu > 0 and -mu % 2 == 1) ? -1. : 1.) * factor *
                            boost::math::spherical_harmonic(n, -mu, theta, -phi);
      auto const actual = rotation(n, 0, mu);
      CHECK(std::abs(actual - expected) < 1e-8);
      CHECK(std::abs(rotation(n, 0, -mu) - std::conj(rotation(n, 0, mu))) < 1e-8);
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
                           -boost::math::spherical_harmonic(n, -mu, theta, -phi) :
                           boost::math::spherical_harmonic(n, -mu, theta, -phi);
    auto const result = std::sqrt(static_cast<Real>(4) * pi / static_cast<Real>(2 * n + 1)) * gamma;
    cache.emplace(Triplet{n, m, mu}, result);
    return result;
  }

  auto const factor = std::exp(Complex(0, chi)) / bnm(n + 1, m - 1);
  auto const first = static_cast<Real>(0.5) * bnm(n + 1, -mu - 1) * std::exp(Complex(0, -phi)) *
                     (static_cast<Real>(1) - std::cos(theta)) *
                     recurrence(theta, phi, chi, n + 1, m - 1, mu + 1, cache);
  auto const second = static_cast<Real>(-0.5) * bnm(n + 1, mu - 1) * std::exp(-Complex(0, -phi)) *
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
  RotationCoefficients rotation(theta, phi, chi);

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
    Complex const actual = rotation(n, m, mu);
    INFO("n: " << n << ", m: " << m << ", mu: " << mu << ", actual: " << actual
               << ", expected: " << expected);
    auto const tolerance =
        std::max(static_cast<Real>(1e-11), std::abs(expected) * static_cast<Real>(1e-8));
    if(std::abs(expected) > tolerance)
      CHECK(std::abs(actual - expected) < tolerance);
    CHECK(std::abs(rotation(n, m, mu) - std::conj(rotation(n, -m, -mu))) < 1e-8);
  }
}

TEST_CASE("Rotation by 0,pi,0 is identity") {
  std::uniform_real_distribution<> idist(0, 10);
  auto const n = 2; // std::uniform_int_distribution<>(1, 10)(*mersenne);

  CAPTURE(RotationCoefficients(0, constant::pi, 0).matrix(n));
  for(t_int i(-static_cast<t_int>(n)); i <= static_cast<t_int>(n); ++i)
    CAPTURE(boost::math::spherical_harmonic(n, i, 0, 0));
  CHECK(RotationCoefficients(0, constant::pi, 0)
            .matrix(n)
            .isApprox(Matrix<t_complex>::Identity(2 * n + 1, 2 * n + 1)));
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

TEST_CASE("Basis rotation from z to zp") {
  Eigen::Matrix<t_real, 3, 1> const z(0, 0, 1);
  Eigen::Matrix<t_real, 3, 1> zp = Eigen::Matrix<t_real, 3, 1>::Random().normalized();
  Eigen::Matrix<t_real, 3, 3> rotation(3, 3);
  rotation.col(2) = zp;
  rotation.col(1) = z.cross(zp).normalized();
  if(rotation.col(1).stableNorm() < 1e-12)
    rotation.col(1) = Eigen::Matrix<t_real, 3, 1>(0, 1, 0);
  rotation.col(0) = rotation.col(1).cross(rotation.col(2)).normalized();

  auto const theta = std::acos(rotation(2, 2));
  auto const phi = std::atan2(rotation(1, 2), rotation(0, 2));
  auto const chi = constant::pi;
  RotationCoefficients rotcoeffs(theta, phi, chi);

  CHECK((rotation * z).isApprox(zp));
  CHECK((rotation.transpose() * zp).isApprox(z));
  CHECK((rotcoeffs.basis_rotation() * z).isApprox(zp));
  CHECK(rotation.isApprox(rotcoeffs.basis_rotation()));

  Eigen::Matrix<t_real, 3, 1> const x0 = Eigen::Matrix<t_real, 3, 1>::Random().normalized();
  Eigen::Matrix<t_real, 3, 1> const x1 = rotation.transpose() * x0;
  auto const sphe0 = to_spherical(x0);
  auto const sphe1 = to_spherical(x1);

  SECTION("Rotation matrix for each n") {
    for(t_uint n(1); n < 10; ++n) {
      Vector<t_complex> rotated(2 * n + 1), original(2 * n + 1);
      for(t_int m(-static_cast<t_int>(n)); m <= static_cast<t_int>(n); ++m) {
        rotated(m + static_cast<t_int>(n)) = rotcoeffs.spherical_harmonic(n, m, sphe1[1], sphe1[2]);
        original(m + static_cast<t_int>(n)) =
            rotcoeffs.spherical_harmonic(n, m, sphe0[1], sphe0[2]);
      }
      CAPTURE(rotated.transpose());
      CAPTURE(original.transpose());
      CAPTURE((rotcoeffs.matrix(n) * rotated).transpose());
      CAPTURE(rotcoeffs.matrix(n));
      CHECK(original.isApprox(rotcoeffs.matrix(n) * rotated));
      CHECK(rotated.isApprox(rotcoeffs.matrix(n).adjoint() * original));
    }
  }

  SECTION("Rotation helper class") {
    auto const nmax = 5;
    Matrix<t_complex> const original = Matrix<t_complex>::Random(nmax * (nmax + 2), 2);
    Rotation const sphe_rot(theta, phi, chi, nmax);
    auto const rotated = sphe_rot.adjoint(original);
    for(t_uint n(1), i(0); n <= nmax; ++n) {
      auto const inc = 2 * n + 1;
      CAPTURE(n);
      CAPTURE(original.col(0).transpose());
      CAPTURE(rotated.col(0).transpose());
      CHECK(rotated.block(i, 0, inc, 2)
                .isApprox(rotcoeffs.matrix(n).adjoint() * original.block(i, 0, inc, 2)));
      CHECK(
          original.block(i, 0, inc, 2).isApprox(rotcoeffs.matrix(n) * rotated.block(i, 0, inc, 2)));
      i += inc;
    }
  }
}
