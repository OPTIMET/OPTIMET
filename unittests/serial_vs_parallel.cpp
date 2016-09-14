#include "catch.hpp"
#include <iostream>

#include "Aliases.h"
#include "Geometry.h"
#include "Scatterer.h"
#include "Solver.h"
#include "Tools.h"
#include "Types.h"
#include "constants.h"
#include <BelosTypes.hpp>

#ifndef OPTIMET_SOLVER
#define OPTIMET_SOLVER "scalapack"
#endif

using namespace optimet;
#ifndef NDEBUG
auto const nHarmonics = 5;
auto const nSpheres = 5;
#else
auto const nHarmonics = 10;
auto const nSpheres = 10;
#endif

TEST_CASE("N spheres") {
  Geometry geometry;
  // spherical coords, ε, μ, radius, nmax
  for(t_uint i(0); i < nSpheres; ++i)
    geometry.pushObject(
        {{static_cast<t_real>(i) * 1.5 * 2e-6, 0, 0}, {0.45e0, 1.1e0}, 0.5 * 2e-6, nHarmonics});

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
#ifdef OPTIMET_BELOS
  solver.belos_parameters()->set("Solver", OPTIMET_SOLVER);
  solver.belos_parameters()->set<int>("Num Blocks", 100);
  solver.belos_parameters()->set("Maximum Iterations", 4000);
  solver.belos_parameters()->set("Convergence Tolerance", 1.0e-14);
#endif
  solver.solve(parallel.scatter_coef, parallel.internal_coef, world);

  optimet::Result serial(&geometry, excitation, nHarmonics);
  auto const serial_context = solver.context().serial();
  auto const comm = world.split(serial_context.is_valid());
  if(serial_context.is_valid()) {
    auto const serial_solver =
        Solver(&geometry, excitation, O3DSolverIndirect, nHarmonics, serial_context);
#ifdef OPTIMET_BELOS
    serial_solver.belos_parameters()->set("Solver", "eigen");
#endif
    serial_solver.solve(serial.scatter_coef, serial.internal_coef, comm);
  }
  // sending serial solution to rest of the world
  broadcast_to_out_of_context(serial.scatter_coef, serial_context, world);
  broadcast_to_out_of_context(serial.internal_coef, serial_context, world);

  REQUIRE(parallel.scatter_coef.rows() == serial.scatter_coef.rows());
  REQUIRE(parallel.scatter_coef.cols() == serial.scatter_coef.cols());
  auto const scatter_tol = 1e-6 * std::max(1., serial.scatter_coef.array().abs().maxCoeff());
  CHECK(parallel.scatter_coef.isApprox(serial.scatter_coef, scatter_tol));
  REQUIRE(parallel.internal_coef.rows() == serial.internal_coef.rows());
  REQUIRE(parallel.internal_coef.cols() == serial.internal_coef.cols());
  auto const internal_tol = 1e-6 * std::max(1., serial.internal_coef.array().abs().maxCoeff());
  CHECK(parallel.internal_coef.isApprox(serial.internal_coef, internal_tol));
}

TEST_CASE("Simultaneous") {
  if(scalapack::global_size() < 3) {
    WARN("Need at least 3 processes to run test");
    return;
  }

  Geometry geometry;
  // spherical coords, ε, μ, radius, nmax
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

  auto const gridsize = scalapack::squarest_largest_grid(scalapack::global_size() - 1);
  Matrix<t_uint> grid_map(gridsize.rows, gridsize.cols);
  for(t_uint i(0), k(1); i < gridsize.rows; ++i)
    for(t_uint j(0); j < gridsize.cols; ++j, ++k)
      grid_map(i, j) = k;

  scalapack::Context const world_context;
  auto const parallel_context = world_context.subcontext(grid_map);
  auto const serial_context = world_context.serial();
  CHECK(parallel_context.is_valid() != serial_context.is_valid());

  CHECK(((serial_context.is_valid() xor parallel_context.is_valid()) or
         (not(parallel_context.is_valid() and serial_context.is_valid()))));

  optimet::Result serial(&geometry, excitation, nHarmonics);
  optimet::Result parallel(&geometry, excitation, nHarmonics);
  auto const comm = mpi::Communicator().split(serial_context.is_valid());

  if(parallel_context.is_valid()) {
    Solver solver(&geometry, excitation, O3DSolverIndirect, nHarmonics, parallel_context);
#ifdef OPTIMET_BELOS
    solver.belos_parameters()->set("Solver", OPTIMET_SOLVER);
    solver.belos_parameters()->set<int>("Num Blocks", 100);
    solver.belos_parameters()->set("Maximum Iterations", 4000);
    solver.belos_parameters()->set("Convergence Tolerance", 1.0e-14);
#endif
    solver.solve(parallel.scatter_coef, parallel.internal_coef, comm);
  } else if(serial_context.is_valid()) {
    Solver solver(&geometry, excitation, O3DSolverIndirect, nHarmonics, serial_context);
#ifdef OPTIMET_BELOS
    solver.belos_parameters()->set("Solver", "eigen");
#endif
    solver.solve(serial.scatter_coef, serial.internal_coef, comm);
  }

  mpi::Communicator world;
  broadcast_to_out_of_context(serial.scatter_coef, serial_context, world);
  broadcast_to_out_of_context(serial.internal_coef, serial_context, world);
  broadcast_to_out_of_context(parallel.scatter_coef, parallel_context, world);
  broadcast_to_out_of_context(parallel.internal_coef, parallel_context, world);

  auto const scatter_tol = 1e-6 * std::max(1., serial.scatter_coef.array().abs().maxCoeff());
  CHECK(serial.scatter_coef.isApprox(parallel.scatter_coef, scatter_tol));
  auto const internal_tol = 1e-6 * std::max(1., serial.internal_coef.array().abs().maxCoeff());
  CHECK(serial.internal_coef.isApprox(parallel.internal_coef, internal_tol));
}
