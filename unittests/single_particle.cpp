#include "catch.hpp"

#include "Types.h"
#include "Geometry.h"
#include "Scatterer.h"
#include "constants.h"
#include "Tools.h"

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
