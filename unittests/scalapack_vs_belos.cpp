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

TEST_CASE("Scalapack vs Belos") {
  Geometry geometry;
  // spherical coords, ε, μ, radius, nmax
  auto const nSpheres = 10;
  auto const nHarmonics = 10;
  for(t_uint i(0); i < 10; ++i)
    geometry.pushObject({{static_cast<t_real>(i) * 1.5, 0, 0}, {0.45e0, 1.1e0}, 0.5, nHarmonics});

  // Create excitation
  auto const wavelength = 14960e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0,
                                90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  Excitation excitation{0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics};
  excitation.populate();
  geometry.update(&excitation);

  Solver solver(&geometry, &excitation, O3DSolverIndirect, nHarmonics);

  optimet::Result scalapack(&geometry, &excitation, nHarmonics);
  solver.belos_parameters()->set("Solver", "scalapack");
  solver.solve(scalapack.scatter_coef, scalapack.internal_coef);

  solver.belos_parameters()->set("Solver", "GMRES");
  solver.belos_parameters()->set("Num Blocks", 1000);
  solver.belos_parameters()->set("Maximum Iterations", 4000);
  solver.belos_parameters()->set("Convergence Tolerance", 1.0e-10);
  optimet::Result belos(&geometry, &excitation, nHarmonics);
  solver.solve(belos.scatter_coef, belos.internal_coef);

  CHECK(belos.scatter_coef.isApprox(belos.scatter_coef));
  CHECK(belos.internal_coef.isApprox(belos.internal_coef));
}
