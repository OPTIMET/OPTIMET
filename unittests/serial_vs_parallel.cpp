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
#include "MatrixBelosSolver.h"
#include "ScalapackSolver.h"
#include "Scatterer.h"
#include "Tools.h"
#include "Types.h"
#include "constants.h"
#include "scalapack/BroadcastToOutOfContext.h"
#ifdef OPTIMET_BELOS
#include <BelosTypes.hpp>
#endif

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
  auto geometry = std::make_shared<Geometry>();
  // spherical coords, ε, μ, radius, nmax
  for(t_uint i(0); i < nSpheres; ++i)
    geometry->pushObject(
        {{static_cast<t_real>(i) * 1.5 * 2e-6, 0, 0}, {0.45e0, 1.1e0}, 0.5 * 2e-6, nHarmonics});

  // Create excitation
  auto const wavelength = 14960e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto const excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry->update(excitation);

  optimet::mpi::Communicator world;
  optimet::Result parallel(geometry, excitation);
#ifdef OPTIMET_BELOS
  optimet::solver::MatrixBelos solver(geometry, excitation, world);
  solver.belos_parameters()->set("Solver", OPTIMET_SOLVER);
  solver.belos_parameters()->set<int>("Num Blocks", 500);
  solver.belos_parameters()->set("Maximum Iterations", 4000);
  solver.belos_parameters()->set("Convergence Tolerance", 1.0e-10);
#else
  optimet::solver::Scalapack const solver(geometry, excitation, world);
#endif
  solver.solve(parallel.scatter_coef, parallel.internal_coef);

  optimet::Result serial(geometry, excitation);
  auto const serial_context = solver.context().serial();
  auto const comm = world.split(serial_context.is_valid());
  if(serial_context.is_valid()) {
    optimet::solver::PreconditionedMatrix const serial_solver(geometry, excitation, comm);
    serial_solver.solve(serial.scatter_coef, serial.internal_coef);
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

  auto geometry = std::make_shared<Geometry>();
  // spherical coords, ε, μ, radius, nmax
  for(t_uint i(0); i < nSpheres; ++i)
    geometry->pushObject({{static_cast<t_real>(i) * 1.5, 0, 0}, {0.45e0, 1.1e0}, 0.5, nHarmonics});

  // Create excitation
  auto const wavelength = 14960e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto const excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, nHarmonics);
  excitation->populate();
  geometry->update(excitation);

  auto const gridsize = scalapack::squarest_largest_grid(scalapack::global_size() - 1);
  Matrix<t_uint> grid_map(gridsize.rows, gridsize.cols);
  for(t_uint i(0), k(1); i < gridsize.rows; ++i)
    for(t_uint j(0); j < gridsize.cols; ++j, ++k)
      grid_map(i, j) = k;

  scalapack::Context const world_context;
  auto const parallel_context = world_context.subcontext(grid_map);
  auto const serial_context = world_context.serial(0);
  if(parallel_context.is_valid())
    REQUIRE(not serial_context.is_valid());

  REQUIRE(((serial_context.is_valid() xor parallel_context.is_valid()) or
           (not(parallel_context.is_valid() and serial_context.is_valid()))));

  optimet::Result serial(geometry, excitation);
  optimet::Result parallel(geometry, excitation);
  auto const serial_comm = mpi::Communicator().split(serial_context.is_valid());
  auto const parallel_comm = mpi::Communicator().split(parallel_context.is_valid());

  if(parallel_context.is_valid()) {
#ifdef OPTIMET_BELOS
    optimet::solver::MatrixBelos solver(geometry, excitation, parallel_comm, parallel_context);
    solver.belos_parameters()->set("Solver", OPTIMET_SOLVER);
    solver.belos_parameters()->set<int>("Num Blocks", 500);
    solver.belos_parameters()->set("Maximum Iterations", 4000);
    solver.belos_parameters()->set("Convergence Tolerance", 1.0e-10);
#else
    optimet::solver::Scalapack const solver(geometry, excitation, parallel_comm, parallel_context);
#endif
    solver.solve(parallel.scatter_coef, parallel.internal_coef);
  } else if(serial_context.is_valid()) {
    optimet::solver::PreconditionedMatrix const solver(geometry, excitation, serial_comm);
    solver.solve(serial.scatter_coef, serial.internal_coef);
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
