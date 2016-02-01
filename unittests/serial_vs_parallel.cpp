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

  optimet::Result parallel(&geometry, &excitation, nHarmonics);
  Solver solver(&geometry, &excitation, O3DSolverIndirect, nHarmonics);
  solver.solve(parallel.scatter_coef, parallel.internal_coef);

  optimet::Result serial(&geometry, &excitation, nHarmonics);
  auto const serial_context = solver.context().serial();
  if(serial_context.is_valid()) {
    solver.context(serial_context).solve(serial.scatter_coef, serial.internal_coef);

    CHECK(parallel.scatter_coef.isApprox(serial.scatter_coef));
    CHECK(parallel.internal_coef.isApprox(serial.internal_coef));
  }
}

TEST_CASE("Simultaneous") {
  if(scalapack::global_size() < 3) {
    WARN("Need at least 3 processes to run test");
    return;
  }

  Geometry geometry;
  // spherical coords, ε, μ, radius, nmax
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

  optimet::Result parallel(&geometry, &excitation, nHarmonics);


  auto const gridsize = scalapack::squarest_largest_grid(scalapack::global_size() - 1);
  Matrix<t_uint> grid_map(gridsize.rows, gridsize.cols);
  for(t_uint i(0), k(1); i < gridsize.rows; ++i)
    for(t_uint j(0); j < gridsize.cols; ++j, ++k)
      grid_map(i, j) = k;
  scalapack::Context const world_context;
  auto const parallel_context = world_context.subcontext(grid_map);
  auto const serial_context = world_context.serial();

  CHECK((
      (serial_context.is_valid() xor parallel_context.is_valid()) or
      (not (parallel_context.is_valid() and serial_context.is_valid()))
  ));

  Solver solver(&geometry, &excitation, O3DSolverIndirect, nHarmonics);
  optimet::Result result(&geometry, &excitation, nHarmonics);

  if(parallel_context.is_valid())
    solver.context(parallel_context).solve(result.scatter_coef, result.internal_coef);
  else if(serial_context.is_valid())
    solver.context(serial_context).solve(result.scatter_coef, result.internal_coef);

  auto const proot = world_context.process_coordinates(1);
  auto const scatter_serial = world_context.broadcast(result.scatter_coef, 0, 0);
  auto const scatter_parallel = world_context.broadcast(result.scatter_coef, proot.row, proot.col);

  auto const internal_serial = world_context.broadcast(result.internal_coef, 0, 0);
  auto const internal_parallel = world_context.broadcast(result.internal_coef, proot.row, proot.col);

  if(parallel_context.is_valid()) {
    CHECK(scatter_serial.isApprox(result.scatter_coef));
    CHECK(internal_serial.isApprox(result.internal_coef));
  }
  CHECK(scatter_serial.isApprox(scatter_parallel));
  CHECK(internal_serial.isApprox(internal_parallel));
}
