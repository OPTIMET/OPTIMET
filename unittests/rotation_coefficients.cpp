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
                            boost::math::spherical_harmonic(n, -mu, theta, phi);
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

TEST_CASE("Rotation matrix") {
  Eigen::Matrix<t_real, 3, 1> const axis = Eigen::Matrix<t_real, 3, 1>::Random();
  auto const basis = RotationCoefficients::basis_rotation(axis);

  SECTION("Axis goes to z") {
    CHECK((basis * basis.transpose()).isApprox(basis.Identity()));
    CHECK((basis.transpose() * basis).isApprox(basis.Identity()));
    CHECK((basis * axis.normalized()).isApprox(axis.Unit(2)));
    CHECK((basis * axis.Unit(2).cross(axis).normalized()).isApprox(axis.Unit(0)));
    CHECK((basis * axis.cross(axis.Unit(2).cross(axis)).normalized()).isApprox(axis.Unit(1)));
  }

  SECTION("Special axis --> z") {
    CHECK(RotationCoefficients::basis_rotation(axis.Unit(2)).isApprox(basis.Identity()));
    Eigen::Matrix<t_real, 3, 3> diag;
    diag.fill(0);
    diag.diagonal() = Eigen::Matrix<t_real, 3, 1>(1, -1, -1);
    CHECK(RotationCoefficients::basis_rotation(-axis.Unit(2)).isApprox(diag));
  }

  SECTION("ϑ vs 2ϑ") {
    auto const axis2 = basis.transpose() * axis;
    auto const basis2 = RotationCoefficients::basis_rotation(axis2);
    CHECK((basis2 * axis2.normalized()).isApprox(axis.Unit(2)));
    CHECK((basis * basis * axis2.normalized()).isApprox(axis.Unit(2)));
    CHECK((basis * basis).row(2).isApprox(basis2.row(2)));
  }
}

TEST_CASE("Rotation angles") {
  Eigen::Matrix<t_real, 3, 1> const axis = Eigen::Matrix<t_real, 3, 1>::Random();
  SECTION("angle <-> matrix coherence") {
    auto const angles = RotationCoefficients::rotation_angles(axis);
    CHECK(RotationCoefficients::basis_rotation(angles).isApprox(
        RotationCoefficients::basis_rotation(axis)));
  }

  SECTION("angles for z-axis") {
    auto const id = RotationCoefficients::rotation_angles(axis.Unit(2));
    CHECK(std::get<0>(id) == Approx(0));
    CHECK(std::get<1>(id) == Approx(0));
    CHECK(std::get<2>(id) == Approx(constant::pi));
    CHECK(RotationCoefficients::basis_rotation(id).isApprox(
        RotationCoefficients::basis_rotation(axis.Unit(2))));

    auto const mid = RotationCoefficients::rotation_angles(-axis.Unit(2));
    CHECK(std::get<0>(mid) == Approx(constant::pi));
    CHECK(std::get<1>(mid) == Approx(0));
    CHECK(std::get<2>(mid) == Approx(0));
    CHECK(RotationCoefficients::basis_rotation(mid).isApprox(
        RotationCoefficients::basis_rotation(-axis.Unit(2))));
  }
}

TEST_CASE("Rotation to axis z is identity") {
  auto const n = 5;
  CHECK(RotationCoefficients(0, constant::pi, 0)
            .matrix(n)
            .isApprox(Matrix<t_complex>::Identity(2 * n + 1, 2 * n + 1)));
  CHECK(RotationCoefficients(0, 0, constant::pi)
            .matrix(n)
            .isApprox(Matrix<t_complex>::Identity(2 * n + 1, 2 * n + 1)));
}

Eigen::Matrix<t_real, 3, 1> to_spherical(Eigen::Matrix<t_real, 3, 1> const &x) {
  auto const r = x.stableNorm();
  if(r < 1e-12)
    return Eigen::Matrix<t_real, 3, 1>::Zero();
  return Eigen::Matrix<t_real, 3, 1>(r, std::acos(x(2) / r), std::atan2(x(1), x(0)));
};

std::function<t_complex(Vector<t_real> const &r)> field(Vector<t_complex> const &coeffs) {
  t_int const nmax = std::lround(std::sqrt(coeffs.rows() + 1) - 1.0);
  assert(nmax * (nmax + 2) == coeffs.rows());
  assert(nmax > 0);

  return [nmax, coeffs](Vector<t_real> const &r) {
    auto const spherical = to_spherical(r);
    t_complex result = 0;
    for(auto n = 1, i = 0; n <= nmax; ++n)
      for(auto m = -n; m <= n; ++m, ++i)
        result +=
            coeffs(i) * RotationCoefficients::spherical_harmonic(n, m, spherical(1), spherical(2));
    return result;
  };
};

template <class T>
typename std::enable_if<std::is_fundamental<T>::value, T>::type pretty(T const &x) {
  return std::abs(x) < 1e-12 ? 0 : x;
}
template <class T>
typename std::enable_if<std::is_fundamental<T>::value, std::complex<T>>::type
pretty(std::complex<T> const &x) {
  return std::complex<T>(pretty(std::real(x)), pretty(std::imag(x)));
}
template <class T> Matrix<typename T::Scalar> pretty(Eigen::MatrixBase<T> const &x) {
  Matrix<typename T::Scalar> result = x;
  for(t_int i(0); i < x.rows(); ++i)
    for(t_int j(0); j < x.rows(); ++j)
      result(i, j) = pretty(result(i, j));
  return result;
}

TEST_CASE("Basis rotation from z to zp") {
  auto const N = 10;
  Eigen::Matrix<t_real, 3, 1> const z = Eigen::Matrix<t_real, 3, 1>::Unit(2);
  Eigen::Matrix<t_real, 3, 1> const zp = Eigen::Matrix<t_real, 3, 1>::Random().normalized();

  RotationCoefficients rotcoeffs(zp);
  auto const basis = rotcoeffs.basis_rotation();

  CHECK((basis * zp).isApprox(z));
  CHECK((basis.transpose() * z).isApprox(zp));

  Eigen::Matrix<t_real, 3, 1> const x0 = Eigen::Matrix<t_real, 3, 1>::Random().normalized();
  auto const sphe0 = to_spherical(x0);
  auto const sphe1 = to_spherical(basis * x0);

  SECTION("Rotation matrix for each n") {
    for(t_uint n(0); n <= N; ++n) {
      Vector<t_complex> rotated(2 * n + 1), original(2 * n + 1);
      for(t_int m(-static_cast<t_int>(n)); m <= static_cast<t_int>(n); ++m) {
        rotated(m + static_cast<t_int>(n)) = rotcoeffs.spherical_harmonic(n, m, sphe1(1), sphe1(2));
        original(m + static_cast<t_int>(n)) =
            rotcoeffs.spherical_harmonic(n, m, sphe0(1), sphe0(2));
      }
      CHECK(original.isApprox(rotcoeffs.matrix(n) * rotated));
      CHECK(rotated.isApprox(rotcoeffs.matrix(n).adjoint() * original));
    }
  }

  Rotation const sphe_rot(zp, N);

  SECTION("Rotation helper class") {
    auto const size = N * (N + 2);
    Matrix<t_complex> const original = Matrix<t_complex>::Random(size, 2);
    auto const rotated = sphe_rot(original);
    SECTION("Check coefficients") {
      for(t_uint n(1), i(0); n <= N; ++n) {
        auto const inc = 2 * n + 1;
        REQUIRE(rotated.block(i, 0, inc, 2)
                    .isApprox(rotcoeffs.matrix(n) * original.block(i, 0, inc, 2)));
        CHECK(original.block(i, 0, inc, 2)
                  .isApprox(rotcoeffs.matrix(n).adjoint() * rotated.block(i, 0, inc, 2)));
        i += inc;
      }
    }

    SECTION("Field rotation") {
      Vector<t_complex> const original = Vector<t_complex>::Ones(size);
      auto const in_field = field(original);
      auto const out_field = field(sphe_rot.transpose(original));
      CHECK(in_field(x0).real() == Approx(out_field(basis * x0).real()));
      CHECK(in_field(x0).imag() == Approx(out_field(basis * x0).imag()));
    }
  }
}

TEST_CASE("Direct vs transpose, conjugate, adjoint") {
  auto const theta = 2e0 * std::uniform_real_distribution<>(0, constant::pi)(*mersenne);
  auto const phi = std::uniform_real_distribution<>(0, constant::pi)(*mersenne);
  auto const chi = std::uniform_real_distribution<>(0, constant::pi)(*mersenne);
  auto const N = 10;
  auto const size = N * (N + 2);
  Rotation const sphe_rot(theta, phi, chi, N);

  Matrix<t_complex> direct(size, size), transpose(size, size), conjugate(size, size),
      adjoint(size, size);
  for(t_int i(0); i < size; ++i) {
    direct.col(i) = sphe_rot(Vector<t_complex>::Unit(size, i));
    transpose.col(i) = sphe_rot.transpose(Vector<t_complex>::Unit(size, i));
    conjugate.col(i) = sphe_rot.conjugate(Vector<t_complex>::Unit(size, i));
    adjoint.col(i) = sphe_rot.adjoint(Vector<t_complex>::Unit(size, i));
  }
  CHECK(transpose.isApprox(direct.transpose()));
  CHECK(conjugate.isApprox(direct.conjugate()));
  CHECK(adjoint.isApprox(direct.adjoint()));
  CHECK((adjoint * direct).isApprox(Matrix<t_complex>::Identity(size, size)));
  CHECK((transpose * conjugate).isApprox(Matrix<t_complex>::Identity(size, size)));
}

TEST_CASE("Additivity, e.g. ϑ vs 2ϑ vs 3ϑ") {
  auto const N = 2;
  auto const size = N * (N + 2);
  SECTION("Three fold rotation along 111") {
    Eigen::Matrix<t_real, 3, 1> const axis(1, 0, 0);
    auto const basis = RotationCoefficients::basis_rotation(axis);
    Rotation const rotation(basis, N);
    Rotation const rotation2(basis * basis, N);

    Matrix<t_complex> direct(size, size), twice(size, size);
    for(t_int i(0); i < size; ++i) {
      twice.col(i) = rotation2(Vector<t_complex>::Unit(size, i));
      direct.col(i) = rotation(Vector<t_complex>::Unit(size, i));
    }
    CHECK((basis * basis * basis).isApprox(basis.Identity()));
    CHECK((direct * direct * direct).isApprox(direct.Identity(size, size)));
    CHECK(twice.isApprox(direct * direct));
  }

  SECTION("Random rotation") {
    Eigen::Matrix<t_real, 3, 1> const axis = Eigen::Matrix<t_real, 3, 1>::Random().normalized();
    auto const basis = RotationCoefficients::basis_rotation(axis);
    Rotation const rotation(basis, N);
    Rotation const rotation2(basis * basis, N);

    Matrix<t_complex> direct(size, size), twice(size, size);
    for(t_int i(0); i < size; ++i) {
      twice.col(i) = rotation2(Vector<t_complex>::Unit(size, i));
      direct.col(i) = rotation(Vector<t_complex>::Unit(size, i));
    }
    CHECK(twice.isApprox(direct * direct));
  }
}
