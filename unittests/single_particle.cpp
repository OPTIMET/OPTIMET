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
#include <iostream>

#include "Aliases.h"
#include "Geometry.h"
#include "Scatterer.h"
#include "PreconditionedMatrix.h"
#include "Tools.h"
#include "Types.h"
#include "constants.h"

using namespace optimet;

TEST_CASE("Add scatterers to geometry") {
  Geometry geometry;
  // spherical coords, ε, μ, radius, nmax
  geometry.pushObject({{1, 0, 0}, {1.1e0, 1.2e0}, 0.5, 2});

  SECTION("Checks first sphere") {
    CHECK(geometry.objects.size() == 1);
    CHECK(geometry.objects.front().vR.rrr == Approx(1));
    CHECK(geometry.objects.front().vR.the == Approx(0));
    CHECK(geometry.objects.front().vR.phi == Approx(0));
    CHECK(geometry.objects.front().elmag.epsilon_r.real() == Approx(1.1));
    CHECK(geometry.objects.front().elmag.epsilon_r.imag() == Approx(0));
    CHECK(geometry.objects.front().elmag.epsilon.real() == Approx(1.1 * consEpsilon0));
    CHECK(geometry.objects.front().elmag.epsilon.imag() == Approx(0));
    CHECK(geometry.objects.front().elmag.mu_r.real() == Approx(1.2));
    CHECK(geometry.objects.front().elmag.mu_r.imag() == Approx(0));
    CHECK(geometry.objects.front().elmag.mu.real() == Approx(1.2 * consMu0));
    CHECK(geometry.objects.front().elmag.mu.imag() == Approx(0));
    CHECK(geometry.objects.front().radius == Approx(0.5));
    CHECK(geometry.objects.front().nMax == 2);
    CHECK(geometry.objects.front().sourceCoef.size() == 2 * Tools::iteratorMax(2));
  }

  SECTION("Add second scatterer") {
    SECTION("No overlap -- far") {
      CHECK_NOTHROW(geometry.pushObject({{3, 0, 0}, {1.1e0, 1.2e0}, 0.5, 2}));
      CHECK(geometry.objects.size() == 2);
      CHECK(geometry.objects.back().vR.rrr == Approx(3.0));
      CHECK(geometry.objects.back().vR.the == Approx(0.0));
      CHECK(geometry.objects.back().vR.phi == Approx(0.0));
    }
    SECTION("Overlap -- touching")
    CHECK_THROWS_AS(geometry.pushObject({{2, 0, 0}, {1.1e0, 1.2e0}, 0.5, 2}), std::runtime_error);
  }
}

TEST_CASE("Two spheres") {
  auto geometry = std::make_shared<Geometry>();
  // spherical coords, ε, μ, radius, nmax
  auto const nHarmonics = 5;
  geometry->pushObject({{0, 0, 0}, {1.0e0, 1.0e0}, 0.5, nHarmonics});
  geometry->pushObject({{1.5, 0, 0}, {1.0e0, 1.0e0}, 0.5, nHarmonics});
  CHECK(geometry->objects.size() == 2);

  // Create excitation
  auto const wavelength = 14960e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto const excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry->update(excitation);

#ifdef OPTIMET_SCALAPACK
  auto const context = scalapack::Context().split(1, scalapack::global_size());
#else
  scalapack::Context const context;
#endif
  auto const S = preconditioned_scattering_matrix(*geometry, excitation);
  auto const Q = source_vector(*geometry, excitation);

  auto const nb = 2 * nHarmonics * (nHarmonics + 2);
  CHECK(S.rows() == S.cols());
  CHECK(S.rows() == nb * geometry->objects.size());
  CHECK(Q.size() == S.cols());

  SECTION("Check transparent <==> identity") {
    CHECK(S.isApprox(Matrix<>::Identity(S.rows(), S.cols())));
  }
  SECTION("Check structure for only one transparent sphere") {
    geometry->objects.front() = {{0, 0, 0}, {10.0e0, 1.0e0}, 0.5, nHarmonics};
    geometry->update(excitation);
    auto const S = preconditioned_scattering_matrix(*geometry, excitation);
    CHECK(S.topLeftCorner(nb, nb).isIdentity());
    CHECK(S.bottomRightCorner(nb, nb).isIdentity());
    CHECK(S.topRightCorner(nb, nb).isZero());
    CHECK(not S.bottomLeftCorner(nb, nb).isZero());
  }
  SECTION("Check structure for two identical spheres") {
    geometry->objects.front() = {{-1, 0, 0}, {10.0e0, 1.0e0}, 0.5, nHarmonics};
    geometry->objects.back() = {{1, 0, 0}, {10.0e0, 1.0e0}, 0.5, nHarmonics};
    geometry->update(excitation);
    auto const S = preconditioned_scattering_matrix(*geometry, excitation);
    CHECK(S.topLeftCorner(nb, nb).isIdentity());
    CHECK(S.bottomRightCorner(nb, nb).isIdentity());
    auto const AB = S.topRightCorner(nb, nb);
    auto const BA = S.bottomLeftCorner(nb, nb);
    CHECK(AB.diagonal().isApprox(BA.diagonal()));
  }
}
