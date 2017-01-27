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
#include <numeric>

#include "Types.h"
#include "mpi/Collectives.h"
#include "mpi/Communicator.h"
#include "scalapack/BroadcastToOutOfContext.h"
#include "scalapack/InitExit.h"
#include "scalapack/LinearSystemSolver.h"

using namespace optimet;

template <class SCALAR> void check(scalapack::Sizes const &grid, t_uint n, t_uint nb) {
  scalapack::Sizes const size = {n, n};
  scalapack::Sizes const blocks = {nb, nb};
  if(mpi::Communicator().size() < grid.rows * grid.cols) {
    WARN("Not enough processes to run test: " << grid.rows << "x" << grid.cols << " < "
                                              << mpi::Communicator().size());
    return;
  }
  scalapack::Context const parallel_context(grid.rows, grid.cols);
  scalapack::Context const serial_context(1, 1);
  scalapack::Matrix<SCALAR> Aserial(serial_context, size, blocks);
  scalapack::Matrix<SCALAR> bserial(serial_context, {size.cols, 1}, blocks);
  if(serial_context.is_valid()) {
    Aserial.local() = Matrix<SCALAR>::Random(size.rows, size.cols) -
                      0.5 * Matrix<SCALAR>::Ones(size.rows, size.cols) +
                      10e0 * Matrix<SCALAR>::Identity(size.rows, size.cols);
    bserial.local() = Matrix<SCALAR>::Random(size.rows, 1).array() - 0.5;
  }

  auto const Aparallel = Aserial.transfer_to(parallel_context);
  auto const bparallel = bserial.transfer_to(parallel_context);

  auto const result = general_linear_system(Aparallel, bparallel);
  REQUIRE(std::get<1>(result) == 0);
  auto const x = std::get<0>(result).transfer_to(serial_context);
  if(serial_context.is_valid())
    CHECK((Aserial.local() * x.local()).isApprox(bserial.local(), 1e-8));
}

TEST_CASE("Compute solver in parallel, check result in serial") {
  SECTION("double") {
    SECTION("1x1") { check<double>({1, 1}, 1024, 64); }
    SECTION("1x3") { check<double>({1, 3}, 1024, 64); }
    SECTION("3x1") { check<double>({3, 1}, 1024, 64); }
    SECTION("3x2") { check<double>({3, 2}, 1024, 64); }
    SECTION("matrix small than block size") { check<double>({2, 2}, 36, 64); }
    SECTION("matrix size equal to block size") { check<double>({2, 2}, 64, 64); }
  }
  SECTION("complex double") {
    SECTION("1x1") { check<std::complex<double>>({1, 1}, 1024, 64); }
    SECTION("1x3") { check<std::complex<double>>({1, 3}, 1024, 64); }
    SECTION("3x1") { check<std::complex<double>>({3, 1}, 1024, 64); }
    SECTION("3x2") { check<std::complex<double>>({3, 2}, 1024, 64); }
  }
}

TEST_CASE("Empty matrix") {
  scalapack::Sizes const grid = {2, 2};
  t_uint const n = 0;
  t_uint const nb = 64;
  scalapack::Sizes const size = {n, n};
  scalapack::Sizes const blocks = {nb, nb};
  if(mpi::Communicator().size() < grid.rows * grid.cols) {
    WARN("Not enough processes to run test: " << grid.rows << "x" << grid.cols << " < "
                                              << mpi::Communicator().size());
    return;
  }
  scalapack::Context const parallel_context(grid.rows, grid.cols);
  scalapack::Context const serial_context(1, 1);
  scalapack::Matrix<double> Aserial(serial_context, size, blocks);
  scalapack::Matrix<double> bserial(serial_context, {size.cols, 1}, blocks);

  auto const Aparallel = Aserial.transfer_to(parallel_context);
  auto const bparallel = bserial.transfer_to(parallel_context);

  auto const result = general_linear_system(Aparallel, bparallel);
  REQUIRE(std::get<1>(result) == 0);
  auto const x = std::get<0>(result).transfer_to(serial_context);
  CHECK(x.size() == 0);
  CHECK(x.local().size() == 0);
}

TEST_CASE("Shuttle data to out-of-context procs") {
  mpi::Communicator const world;
  auto const input = world.broadcast(Vector<t_real>::Random(6).eval());
  SECTION("All procs in context") {
    scalapack::Context const context(1, world.size());
    auto const local = Vector<t_real>::Random(input.size()).eval();
    auto inout = local;
    broadcast_to_out_of_context(inout, context, world);
    CHECK(inout.isApprox(local));
  }

  for(size_t i(1); i < world.size() - 1; ++i) {
    std::ostringstream sstr;
    sstr << i << " procs are missing from context";
    SECTION(sstr.str()) {
      scalapack::Context const context(1, world.size() - i);
      auto const is_context_root = context.is_valid() and context.row() == 0 and context.col() == 0;
      auto const root_rank = world.all_reduce(is_context_root ? world.rank() : 0, MPI_SUM);
      CHECK(i == world.all_reduce<int>(context.is_valid() ? 0 : 1, MPI_SUM));
      auto const local = Vector<t_real>::Random(input.size()).eval();
      auto inout = world.rank() == root_rank ? input : local;
      broadcast_to_out_of_context(inout, context, world);
      CHECK(inout.isApprox(world.rank() == root_rank or not context.is_valid() ? input : local));
    }
  }
}
