#include <iostream>
#include "catch.hpp"

#include "Types.h"
#include "Geometry.h"
#include "Scatterer.h"
#include "constants.h"
#include "Tools.h"
#include "Solver.h"
#include "Aliases.h"

using namespace optimet;

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

  optimet::Result parallel(&geometry, &excitation, nHarmonics);
  Solver solver(&geometry, &excitation, O3DSolverIndirect, nHarmonics);
  solver.solve(parallel.scatter_coef, parallel.internal_coef);

  optimet::Result serial(&geometry, &excitation, nHarmonics);
  solver.context({1, 1}).solve(parallel.scatter_coef, parallel.internal_coef);
}
