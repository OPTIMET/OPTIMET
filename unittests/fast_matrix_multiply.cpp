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

#include "Aliases.h"
#include "FastMatrixMultiply.h"
#include "Geometry.h"
#include "HarmonicsIterator.h"
#include "Result.h"
#include "PreconditionedMatrix.h"
#include "Tools.h"
#include "catch.hpp"
#include <iostream>

ElectroMagnetic const silicon{13.1, 1.0};
auto const wavenumber = 2 * optimet::constant::pi / (1200 * 1e-9);
auto const nHarmonics = 5;
auto const radius = 500e-9;

TEST_CASE("No objects") {
  using namespace optimet;
  FastMatrixMultiply fmm{wavenumber, std::vector<Scatterer>()};
  CHECK(fmm.couplings().size() == 0);
  CHECK(fmm.rows() == fmm.cols());
  CHECK(fmm.rows() == 0);
  CHECK((fmm * Vector<t_complex>::Random(0)).size() == 0);
}

TEST_CASE("Single object") {
  using namespace optimet;
  std::vector<Scatterer> const scatterers{{{0, 0, 0}, silicon, radius, nHarmonics}};

  FastMatrixMultiply fmm{wavenumber, scatterers};
  CHECK(fmm.couplings().size() == 1);
  CHECK(fmm.couplings().back().first == 0);
  CHECK(fmm.couplings().back().second == 0);
  CHECK(fmm.rows() == fmm.cols());
  CHECK(fmm.rows() == 2 * nHarmonics * (nHarmonics + 2));

  Vector<t_complex> const input = Vector<t_complex>::Random(2 * nHarmonics * (nHarmonics + 2));
  auto const result = fmm(input);
  CHECK(result.isApprox(input));
}

TEST_CASE("Transparent objects") {
  // For transparent objects, the coaxial transpose rotation mess should be zero
  // Only the left hand side of Eq 106 in Gumerov, Duraiswami (2007) is left (e.g identity part).
  using namespace optimet;
  Eigen::Matrix<t_real, 3, 1> const direction = Vector<t_real>::Random(3).normalized();
  ElectroMagnetic const bground;
  std::vector<Scatterer> scatterers{
      {direction.Zero(3), bground, radius, nHarmonics},
      {direction * 3 * radius * 1.00001, bground, 2 * radius, nHarmonics},
  };

  optimet::FastMatrixMultiply fmm(bground, wavenumber, scatterers);
  CHECK(fmm.couplings().size() == 4);
  CHECK(fmm.rows() == fmm.cols());
  CHECK(fmm.rows() == scatterers.size() * 2 * nHarmonics * (nHarmonics + 2));

  Vector<t_complex> const input = Vector<t_complex>::Random(fmm.cols());
  auto const actual = fmm * input;
  CHECK(actual.size() == input.size());
  CHECK(actual.isApprox(input));
  CHECK(fmm.transpose(input).isApprox(input));
  CHECK(fmm.adjoint(input).isApprox(input));
  CHECK(fmm.conjugate(input).isApprox(input));
}

TEST_CASE("Transpose/conjugate/adjoint of the fast matrix multiply") {
  using namespace optimet;
  auto const radius = 500.0e-9;
  Eigen::Matrix<t_real, 3, 1> const direction(2, 1, 1); // = Vector<t_real>::Random(3).normalized();
  auto geometry = std::make_shared<Geometry>();
  geometry->pushObject({{0, 0, 0}, silicon, radius, nHarmonics});
  geometry->pushObject({direction * 3 * radius * 1.500001, silicon, 2 * radius, nHarmonics});

  auto const wavelength = 1490.0e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry->update(excitation);

  optimet::FastMatrixMultiply fmm(geometry->bground, excitation->omega() / constant::c,
                                  geometry->objects);

  auto const size = geometry->objects.size() * nHarmonics * (nHarmonics + 2) * 2;
  Matrix<t_complex> expected(size, size), transpose(size, size), conjugate(size, size),
      adjoint(size, size);
  for(t_int i(0); i < size; ++i) {
    expected.col(i) = fmm(Vector<t_complex>::Unit(size, i));
    transpose.col(i) = fmm.transpose(Vector<t_complex>::Unit(size, i));
    conjugate.col(i) = fmm.conjugate(Vector<t_complex>::Unit(size, i));
    adjoint.col(i) = fmm.adjoint(Vector<t_complex>::Unit(size, i));
  }
  CAPTURE(transpose.row(0));
  CAPTURE(expected.col(0).transpose());
  CHECK(transpose.isApprox(expected.transpose()));
  CHECK(conjugate.isApprox(expected.conjugate()));
  CHECK(adjoint.isApprox(expected.adjoint()));
}

TEST_CASE("Standard vs Fast matrix multiply") {
  using namespace optimet;
  auto const radius = 500.0e-9;
  Eigen::Matrix<t_real, 3, 1> const direction = Vector<t_real>::Random(3).normalized();
  auto geometry = std::make_shared<Geometry>();
  geometry->pushObject({{0, 0, 0}, silicon, radius, nHarmonics});
  geometry->pushObject({direction * 3 * radius * 1.500001, silicon, 2 * radius, nHarmonics});

  auto const wavelength = 1490.0e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry->update(excitation);

  SECTION("Two particles") {
    optimet::FastMatrixMultiply fmm(geometry->bground, excitation->omega() / constant::c,
                                    geometry->objects);

    auto const S = preconditioned_scattering_matrix(*geometry, excitation);
    auto const size = S.cols();
    for(t_int i(0); i < size; ++i) {
      auto const input = Vector<t_complex>::Unit(size, i);
      Vector<t_complex> const expected = S * input;
      Vector<t_complex> const actual = fmm * input;
      CAPTURE(actual.transpose());
      REQUIRE(expected.isApprox(actual));
    }
  }

  SECTION("Three particles") {
    ElectroMagnetic const other{9, 2.0};
    auto const x = Eigen::Matrix<t_real, 3, 1>::Unit(0).eval();
    geometry->pushObject(
        {direction * 1.5 * radius * 1.500001 + x * radius * 8, other, 0.5 * radius, nHarmonics});
    optimet::FastMatrixMultiply fmm(geometry->bground, excitation->omega() / constant::c,
                                    geometry->objects);
    auto const S = preconditioned_scattering_matrix(*geometry, excitation);
    auto const size = S.cols();

    for(t_int i(0); i < size; ++i) {
      auto const input = Vector<t_complex>::Unit(size, i);
      Vector<t_complex> const expected = S * input;
      Vector<t_complex> const actual = fmm * input;
      CHECK(expected.isApprox(actual));
    }
  }
}

TEST_CASE("Ranged matrices") {
  using namespace optimet;
  auto const wavelength = 1490.0e-9;
  auto const radius = 500.0e-9;
  Eigen::Matrix<t_real, 3, 1> const direction = Vector<t_real>::Random(3).normalized();
  ElectroMagnetic const other{9, 2.0};

  std::vector<Scatterer> scatterers;
  scatterers.emplace_back(Vector<t_real>::Zero(3), silicon, radius, nHarmonics);
  scatterers.emplace_back(direction * 3 * radius * 1.500001, silicon, 2 * radius, nHarmonics);
  auto const x = Eigen::Matrix<t_real, 3, 1>::Unit(0).eval();
  scatterers.emplace_back(direction * 1.5 * radius * 1.500001 + x * radius * 8, other, 0.5 * radius,
                          nHarmonics);

  auto const omega = 2 * constant::pi / wavelength;
  auto const wavenumber = omega / constant::c;
  ElectroMagnetic const bground;

  auto const N = nHarmonics * (nHarmonics + 2);
  Matrix<t_int> splitting = -Matrix<t_int>::Ones(scatterers.size(), scatterers.size());
  splitting(0, 0) = 1;
  splitting.col(0).tail(scatterers.size() - 1).fill(2);
  splitting.rightCols(scatterers.size() - 1).fill(3);
  CHECK(not(splitting.array() == -1).any());
  optimet::FastMatrixMultiply whole(wavenumber, scatterers);
  optimet::FastMatrixMultiply partA(wavenumber, scatterers, splitting.array() == 1);
  optimet::FastMatrixMultiply partB(wavenumber, scatterers, splitting.array() == 2);
  optimet::FastMatrixMultiply partC(wavenumber, scatterers, splitting.array() == 3);

  CHECK(whole.rows() == whole.cols());
  CHECK(whole.rows() == 2 * N * scatterers.size());
  CHECK(partA.rows() == partA.cols());
  CHECK(partA.rows() == 2 * N);
  CHECK(partB.cols() == 2 * N);
  CHECK(partB.rows() == 2 * N * 2);
  CHECK(partC.cols() == 2 * N * 2);
  CHECK(partC.rows() == 2 * N * 3);

  SECTION("Transpose, conjugate, adjoint") {
    for(auto const &part : {partA, partB, partC}) {
      Matrix<t_complex> direct(part.rows(), part.cols()), transpose(part.cols(), part.rows()),
          adjoint(part.cols(), part.rows()), conjugate(part.rows(), part.cols());
      for(t_int i(0); i < part.cols(); ++i) {
        direct.col(i) = part(Vector<t_complex>::Unit(part.cols(), i));
        conjugate.col(i) = part.conjugate(Vector<t_complex>::Unit(part.cols(), i));
      }
      for(t_int i(0); i < part.rows(); ++i) {
        transpose.col(i) = part.transpose(Vector<t_complex>::Unit(part.rows(), i));
        adjoint.col(i) = part.adjoint(Vector<t_complex>::Unit(part.rows(), i));
      }
      CHECK(transpose.isApprox(direct.transpose()));
      CHECK(conjugate.isApprox(direct.conjugate()));
      CHECK(adjoint.isApprox(direct.adjoint()));
    }
  }

  SECTION("Whole against parts") {
    Matrix<t_complex> matrix_whole(whole.rows(), whole.cols()),
        matrix_A(partA.rows(), partA.cols()), matrix_B(partB.rows(), partB.cols()),
        matrix_C(partC.rows(), partC.cols());
    for(t_uint i(0); i < whole.rows(); ++i) {
      matrix_whole.col(i) = whole(Vector<t_complex>::Unit(whole.cols(), i));
      if(i < matrix_A.cols())
        matrix_A.col(i) = partA(Vector<t_complex>::Unit(partA.cols(), i));
      if(i < matrix_B.cols())
        matrix_B.col(i) = partB(Vector<t_complex>::Unit(partB.cols(), i));
      if(i < matrix_C.cols())
        matrix_C.col(i) = partC(Vector<t_complex>::Unit(partC.cols(), i));
    }

    Matrix<t_complex> reconstructed = Matrix<t_complex>::Zero(whole.rows(), whole.cols());
    reconstructed.block(0, 0, partA.rows(), partA.cols()) = matrix_A;
    reconstructed.block(partA.rows(), 0, partB.rows(), partB.cols()) = matrix_B;
    reconstructed.block(0, partA.cols(), partC.rows(), partC.cols()) = matrix_C;
    CHECK(reconstructed.isApprox(matrix_whole));
  }
}
