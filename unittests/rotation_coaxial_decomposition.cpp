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

// returns unit matrix
Matrix<t_complex> unit_matrix(t_uint nx, t_uint ny, t_uint i, t_uint j) {
  Matrix<t_complex> result(nx, ny);
  result.fill(0);
  result(i, j) = 1;
  return result;
}

Vector<t_complex> rotcoax(t_real wavenumber, t_real tz, Matrix<t_complex> const &input) {
  Vector<t_complex> output(input.size());
  Eigen::Map<Matrix<t_complex>>(output.data(), input.rows(), input.cols()) =
      rotation_coaxial_decomposition(wavenumber, tz, input);
  return output;
};

Vector<t_complex> rotcoax_transpose(t_real wavenumber, t_real tz, Matrix<t_complex> const &input) {
  Vector<t_complex> output(input.size());
  Eigen::Map<Matrix<t_complex>>(output.data(), input.rows(), input.cols()) =
      rotation_coaxial_decomposition_transpose(wavenumber, tz, input);
  return output;
};

TEST_CASE("Transpose operation") {
  auto const wavenumber = 10;
  auto const tz = 7;
  auto const N = 10;
  auto const size = N * (N + 2) + 1;

  Matrix<t_complex> actual(2 * size, 2 * size), expected(2 * size, 2 * size);
  for(t_int i(0); i < size; ++i) {
    expected.col(i) = rotcoax(wavenumber, tz, unit_matrix(size, 2, i, 0));
    expected.col(i + size) = rotcoax(wavenumber, tz, unit_matrix(size, 2, i, 1));
    actual.col(i) = rotcoax_transpose(wavenumber, tz, unit_matrix(size, 2, i, 0));
    actual.col(i + size) = rotcoax_transpose(wavenumber, tz, unit_matrix(size, 2, i, 1));
  }
  CAPTURE(actual.row(0));
  CAPTURE(expected.col(0).transpose());
  CHECK(actual.transpose().isApprox(expected));
}
