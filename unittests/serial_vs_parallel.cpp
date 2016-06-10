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

TEST_CASE("N spheres") {
  Geometry geometry;
  // spherical coords, ε, μ, radius, nmax
  auto const nHarmonics = 10;
  auto const nSpheres = 10;
  for(t_uint i(0); i < nSpheres; ++i)
    geometry.pushObject({{static_cast<t_real>(i) * 1.5, 0, 0}, {0.45e0, 1.1e0}, 0.5, nHarmonics});

  // Create excitation
  auto const wavelength = 14960e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto const excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry.update(excitation);

  optimet::mpi::Communicator world;
  optimet::Result parallel(&geometry, excitation, nHarmonics);
  Solver solver(&geometry, excitation, O3DSolverIndirect, nHarmonics);
  // called with world communicator: world has solution
  solver.solve(parallel.scatter_coef, parallel.internal_coef, world);

  optimet::Result serial(&geometry, excitation, nHarmonics);
  auto const serial_context = solver.context().serial();
  if(serial_context.is_valid()) {
    auto const serial_solver =
        Solver(&geometry, excitation, O3DSolverIndirect, nHarmonics, serial_context);
    // called without world communicator: only serial context has solution
    serial_solver.solve(serial.scatter_coef, serial.internal_coef);
  }
  // sending serial solution to rest of the world
  broadcast_to_out_of_context(serial.scatter_coef, serial_context, world);
  broadcast_to_out_of_context(serial.internal_coef, serial_context, world);

  REQUIRE(parallel.scatter_coef.rows() == serial.scatter_coef.rows());
  REQUIRE(parallel.scatter_coef.cols() == serial.scatter_coef.cols());
  CHECK(parallel.scatter_coef.isApprox(serial.scatter_coef));
  REQUIRE(parallel.internal_coef.rows() == serial.internal_coef.rows());
  REQUIRE(parallel.internal_coef.cols() == serial.internal_coef.cols());
  CHECK(parallel.internal_coef.isApprox(serial.internal_coef));
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
    geometry.pushObject({{static_cast<t_real>(i) * 1.5, 0, 0}, {0.45e0, 1.1e0}, 0.5,
    nHarmonics});

  // Create excitation
  auto const wavelength = 14960e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi /
  180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto const excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry.update(excitation);

  auto const gridsize = scalapack::squarest_largest_grid(scalapack::global_size() - 1);
  Matrix<t_uint> grid_map(gridsize.rows, gridsize.cols);
  for(t_uint i(0), k(1); i < gridsize.rows; ++i)
    for(t_uint j(0); j < gridsize.cols; ++j, ++k)
      grid_map(i, j) = k;

  scalapack::Context const world_context;
  auto const parallel_context = world_context.subcontext(grid_map);
  auto const serial_context = world_context.serial();

  CHECK(((serial_context.is_valid() xor parallel_context.is_valid()) or
         (not(parallel_context.is_valid() and serial_context.is_valid()))));

  Solver solver(&geometry, excitation, O3DSolverIndirect, nHarmonics);
  optimet::Result serial(&geometry, excitation, nHarmonics);
  optimet::Result parallel(&geometry, excitation, nHarmonics);

  if(parallel_context.is_valid()) {
    Solver solver(&geometry, excitation, O3DSolverIndirect, nHarmonics, parallel_context);
    solver.solve(parallel.scatter_coef, parallel.internal_coef);
  } else if(serial_context.is_valid()) {
    Solver solver(&geometry, excitation, O3DSolverIndirect, nHarmonics, serial_context);
    solver.solve(serial.scatter_coef, serial.internal_coef);
  }

  mpi::Communicator world;
  broadcast_to_out_of_context(serial.scatter_coef, serial_context, world);
  broadcast_to_out_of_context(serial.internal_coef, serial_context, world);
  broadcast_to_out_of_context(parallel.scatter_coef, parallel_context, world);
  broadcast_to_out_of_context(parallel.internal_coef, parallel_context, world);

  CHECK(serial.scatter_coef.isApprox(parallel.scatter_coef));
  CHECK(serial.internal_coef.isApprox(parallel.internal_coef));
}
