#include <iostream>
#include "catch.hpp"

#include "types.h"
#include "Geometry.h"
#include "Scatterer.h"
#include "constants.h"
#include "Tools.h"
#include "Solver.h"
#include "Aliases.h"

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
    CHECK(geometry.objects.front().elmag.epsilon.real() ==
          Approx(1.1 * consEpsilon0));
    CHECK(geometry.objects.front().elmag.epsilon.imag() == Approx(0));
    CHECK(geometry.objects.front().elmag.mu_r.real() == Approx(1.2));
    CHECK(geometry.objects.front().elmag.mu_r.imag() == Approx(0));
    CHECK(geometry.objects.front().elmag.mu.real() == Approx(1.2 * consMu0));
    CHECK(geometry.objects.front().elmag.mu.imag() == Approx(0));
    CHECK(geometry.objects.front().radius == Approx(0.5));
    CHECK(geometry.objects.front().nMax == 2);
    CHECK(geometry.objects.front().sourceCoef.size() ==
          2 * Tools::iteratorMax(2));
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
    CHECK_THROWS_AS(geometry.pushObject({{2, 0, 0}, {1.1e0, 1.2e0}, 0.5, 2}),
                    std::runtime_error);
  }
}

TEST_CASE("Two spheres") {
  Geometry geometry;
  // spherical coords, ε, μ, radius, nmax
  auto const nHarmonics = 5;
  geometry.pushObject({{0, 0, 0}, {1.0e0, 1.0e0}, 0.5, nHarmonics});
  geometry.pushObject({{1.5, 0, 0}, {1.0e0, 1.0e0}, 0.5, nHarmonics});
  CHECK(geometry.objects.size() == 2);

  // Create excitation
  auto const wavelength = 14960e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0,
                                90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  Excitation excitation{0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics};
  excitation.populate();
  geometry.update(&excitation);

  Solver solver(&geometry, &excitation, O3DSolverIndirect, nHarmonics);

  auto const nb = 2 * nHarmonics * (nHarmonics + 2);
  CHECK(solver.S.rows() == solver.S.cols());
  CHECK(solver.S.rows() == nb * geometry.objects.size());
  CHECK(solver.Q.size() == solver.S.cols());

  SECTION("Check transparent <==> identity") {
    solver.populate();
    CHECK(solver.S.isApprox(
        Matrix<>::Identity(solver.S.rows(), solver.S.cols())));
  }
  SECTION("Check structure for only one transparent sphere") {
    geometry.objects.front() = {{0, 0, 0}, {10.0e0, 1.0e0}, 0.5, nHarmonics};
    geometry.update(&excitation);
    solver.populate();
    CHECK(solver.S.topLeftCorner(nb, nb).isIdentity());
    CHECK(solver.S.bottomRightCorner(nb, nb).isIdentity());
    CHECK(solver.S.topRightCorner(nb, nb).isZero());
    CHECK(not solver.S.bottomLeftCorner(nb, nb).isZero());
  }
  SECTION("Check structure for two identical spheres") {
    geometry.objects.front() = {{-1, 0, 0}, {10.0e0, 1.0e0}, 0.5, nHarmonics};
    geometry.objects.back() = {{1, 0, 0}, {10.0e0, 1.0e0}, 0.5, nHarmonics};
    geometry.update(&excitation);
    solver.populate();
    CHECK(solver.S.topLeftCorner(nb, nb).isIdentity());
    CHECK(solver.S.bottomRightCorner(nb, nb).isIdentity());
    auto const AB = solver.S.topRightCorner(nb, nb);
    auto const BA = solver.S.bottomLeftCorner(nb, nb);
    CHECK(AB.diagonal().isApprox(BA.diagonal()));
  }
}
