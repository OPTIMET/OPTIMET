#include "catch.hpp"

#include "Coefficients.h"
#include "RotationCoaxialDecomposition.h"
#include "Types.h"
#include "constants.h"

using namespace optimet;

TEST_CASE("Brute check") {
  using optimet::coefficient::a;
  auto const nmax = 4;
  Matrix<t_complex> input(nmax * (nmax + 2), 2);
  input.fill(0);
  SECTION("Nothing in n±1") {
    SECTION("n=1") {
      input(0, 0) = 2;
      input(0, 1) = t_complex(0, 1);

      for(auto const k : {1.0, 1.5}) {
        auto out = rotation_coaxial_decomposition(k, 2, input);
        CHECK(out(0, 0).real() == Approx(2 + k));
        CHECK(out(0, 0).imag() == Approx(0));
        CHECK(out(0, 1).real() == Approx(0));
        CHECK(out(0, 1).imag() == Approx(1 - 2 * k));
      }
    }
    SECTION("n=2") {
      input(3, 0) = 2;
      input(3, 1) = t_complex(0, 1);

      auto out = rotation_coaxial_decomposition(1, 6, input);
      CHECK(out(3, 0).real() == Approx(4));
      CHECK(out(3, 0).imag() == Approx(0));
      CHECK(out(3, 1).real() == Approx(0));
      CHECK(out(3, 1).imag() == Approx(-3));
    }
  }

  SECTION("Something in n - 1") {
    input(0, 0) = 2;
    input(0, 1) = t_complex(0, 1);

    auto out = rotation_coaxial_decomposition(1, 6, input);
    CHECK(out(4, 0).real() == Approx(3 * a<>(1, -1) * 2));
    CHECK(out(4, 0).imag() == Approx(0));
    CHECK(out(4, 1).real() == Approx(0));
    CHECK(out(4, 1).imag() == Approx(3 * a<>(1, -1)));
  }

  SECTION("Something in n + 1") {
    input(4, 0) = 2;
    input(4, 1) = t_complex(0, 1);

    for(auto const k : {1.0, 1.5}) {
      auto out = rotation_coaxial_decomposition(k, 2, input);
      CHECK(out(0, 0).real() == Approx(a<>(1, -1) * 2 * k));
      CHECK(out(0, 0).imag() == Approx(0));
      CHECK(out(0, 1).real() == Approx(0));
      CHECK(out(0, 1).imag() == Approx(a<>(1, -1) * k));
    }
  }
}

TEST_CASE("Transpose operation") {
  auto const wavenumber = 10;
  auto const tz = 7;
  auto const N = 10;
  auto const size = N * (N + 2);

  Matrix<t_complex> actual(2 * size, 2 * size), expected(2 * size, 2 * size);
  for(t_int i(0); i < size; ++i) {
    Matrix<t_complex> input = Matrix<t_complex>::Zero(size, 2);
    input.col(0) = Vector<t_complex>::Unit(size, i);

    auto const phi = rotation_coaxial_decomposition(wavenumber, tz, input);
    expected.col(i).head(size) = phi.col(0);
    expected.col(i).tail(size) = phi.col(1);

    auto const phiT = rotation_coaxial_decomposition_transpose(wavenumber, tz, input);
    actual.col(i).head(size) = phiT.col(0);
    actual.col(i).tail(size) = phiT.col(1);

    input.col(1) = Vector<t_complex>::Unit(size, i);
    input.col(0) = Vector<t_complex>::Zero(size);

    auto const psi = rotation_coaxial_decomposition(wavenumber, tz, input);
    expected.col(i + size).head(size) = psi.col(0);
    expected.col(i + size).tail(size) = psi.col(1);

    auto const psiT = rotation_coaxial_decomposition_transpose(wavenumber, tz, input);
    actual.col(i + size).head(size) = psiT.col(0);
    actual.col(i + size).tail(size) = psiT.col(1);
  }
  CHECK(actual.transpose().isApprox(expected));
}
