#include "catch.hpp"
#include <iostream>

#include "Aliases.h"
#include "Geometry.h"
#include "Scatterer.h"
#include "Solver.h"
#include "Tools.h"
#include "Types.h"
#include "constants.h"

using namespace optimet;

TEST_CASE("Scalapack vs Belos") {
  Geometry geometry;
  // spherical coords, ε, μ, radius, nmax
  auto const nSpheres = 10;
  auto const nHarmonics = 10;
  for(t_uint i(0); i < 10; ++i) {
    t_real const ii(i);
    geometry.pushObject(
        {{ii * 1.5, 0, 0}, {0.45e0 + 0.1 * ii, 1.1e0}, 0.5 + 0.01 * ii, nHarmonics});
  }

  // Create excitation
  auto const wavelength = 14960e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry.update(excitation);

  Solver solver(&geometry, excitation, O3DSolverIndirect, nHarmonics);

  optimet::Result scalapack(&geometry, excitation, nHarmonics);
  solver.belos_parameters()->set("Solver", "scalapack");
  solver.solve(scalapack.scatter_coef, scalapack.internal_coef);

  // known to fail: "CGPOLY", "FLEXIBLE GMRES", "RECYCLING CG", "RCG", "PCPG", "MINRES", "LSQR",
  // "SEED CG",
  auto const names = {
      "BICGSTAB",
      "BLOCK GMRES",
      "CG",                 // "PSEUDO BLOCK CG",
      "GMRES",              // "PSEUDO BLOCK GMRES",
      "GMRESPOLY",          // "HYBRID BLOCK GMRES", "SEED GMRES"
      "PSEUDO BLOCK TFQMR", // "PSEUDO BLOCK TRANSPOSE-FREE QMR",
      "GCRODR",             // "RECYCLING GMRES",
      "STOCHASTIC CG",      // "PSEUDO BLOCK STOCHASTIC CG"
      "TFQMR",              // "TRANSPOSE-FREE QMR",
      "BLOCK CG",
      "BLOCK GMRES",
      "FIXED POINT",
  };
  for(auto const name : names) {
    SECTION(name) {
      solver.belos_parameters()->set("Solver", name);
      solver.belos_parameters()->set("Num Blocks", 1000);
      solver.belos_parameters()->set("Maximum Iterations", 4000);
      solver.belos_parameters()->set("Convergence Tolerance", 1.0e-10);
      optimet::Result belos(&geometry, excitation, nHarmonics);
      solver.solve(belos.scatter_coef, belos.internal_coef);

      CHECK(belos.scatter_coef.isApprox(belos.scatter_coef));
      CHECK(belos.internal_coef.isApprox(belos.internal_coef));
    }
  }
}
