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
#include "scalapack/Context.h"
#include "scalapack/InitExit.h"
#include "scalapack/Matrix.h"

using namespace optimet;

TEST_CASE("Creates a matrix in 1x1 context") {
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(1, 1);
  REQUIRE(context.is_valid() == (rank == 0));
  scalapack::Matrix<> matrix(context, {64, 64}, {20, 10});

  CHECK(matrix.local().rows() == (context.is_valid() ? 64 : 0));
  CHECK(matrix.local().cols() == (context.is_valid() ? 64 : 0));
  CHECK(matrix.blacs()[0] == 1);
  if(context.is_valid())
    CHECK(matrix.blacs()[1] == *context);
  CHECK(matrix.blacs()[2] == 64);
  CHECK(matrix.blacs()[3] == 64);
  CHECK(matrix.blacs()[4] == 20);
  CHECK(matrix.blacs()[5] == 10);
  CHECK(matrix.blacs()[6] == 0);
  CHECK(matrix.blacs()[7] == 0);
  auto const leading = matrix.local().IsRowMajor ? matrix.local().cols() : matrix.local().rows();
  CHECK(matrix.blacs()[8] == std::max(leading, static_cast<decltype(leading)>(1)));
}

TEST_CASE("Creates a matrix in 1x2 context") {
  if(scalapack::global_size() < 2)
    return;
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(1, 2);
  scalapack::Matrix<> matrix(context, {64, 64}, {16, 8});

  CHECK(matrix.local().rows() == (context.is_valid() ? 64 : 0));
  CHECK(matrix.local().cols() == (context.is_valid() ? 32 : 0));
  CHECK(matrix.blacs()[0] == 1);
  CHECK(matrix.blacs()[1] == (context.is_valid() ? *context : -1));
  CHECK(matrix.blacs()[2] == 64);
  CHECK(matrix.blacs()[3] == 64);
  CHECK(matrix.blacs()[4] == 16);
  CHECK(matrix.blacs()[5] == 8);
  CHECK(matrix.blacs()[6] == 0);
  CHECK(matrix.blacs()[7] == 0);
  auto const leading = matrix.local().IsRowMajor ? matrix.local().cols() : matrix.local().rows();
  CHECK(matrix.blacs()[8] == std::max(leading, static_cast<decltype(leading)>(1)));
}

void check_creation(t_uint n, t_uint m) {
  if(scalapack::global_size() < n * m)
    return;
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(n, m);
  scalapack::Sizes const size = {1024, 1024};
  scalapack::Sizes const blocks = {31, 65};
  scalapack::Matrix<> matrix(context, size, blocks);

  if(context.is_valid()) {
    REQUIRE(matrix.local().rows() > 0);
    REQUIRE(matrix.local().cols() > 0);
  } else {
    REQUIRE(matrix.local().rows() == 0);
    REQUIRE(matrix.local().cols() == 0);
  }

  auto const split = mpi::Communicator().split(context.is_valid());
  if(not context.is_valid())
    return;

  REQUIRE(split.size() == n * m);
  auto const rows = split.all_gather(matrix.local().rows());
  auto const cols = split.all_gather(matrix.local().cols());
  auto n_elements = 0;
  for(std::size_t i(0); i < split.size(); ++i)
    n_elements += rows[i] * cols[i];
  CHECK(static_cast<t_uint>(n_elements) == size.rows * size.cols);

  CHECK(matrix.blacs()[0] == 1);
  CHECK(matrix.blacs()[1] == *context);
  CHECK(static_cast<t_uint>(matrix.blacs()[2]) == size.rows);
  CHECK(static_cast<t_uint>(matrix.blacs()[3]) == size.cols);
  CHECK(static_cast<t_uint>(matrix.blacs()[4]) == blocks.rows);
  CHECK(static_cast<t_uint>(matrix.blacs()[5]) == blocks.cols);
  CHECK(matrix.blacs()[6] == 0);
  CHECK(matrix.blacs()[7] == 0);
  auto const leading = matrix.local().IsRowMajor ? matrix.local().cols() : matrix.local().rows();
  CHECK(matrix.blacs()[8] == std::max(leading, static_cast<decltype(leading)>(1)));
}

TEST_CASE("Creates a matrice in nxm context") {
  SECTION("1x1") { check_creation(1, 1); }
  SECTION("2x1") { check_creation(2, 1); }
  SECTION("1x2") { check_creation(1, 2); }
  SECTION("2x2") { check_creation(2, 2); }
  SECTION("3x1") { check_creation(3, 1); }
  SECTION("3x2") { check_creation(3, 2); }
}

TEST_CASE("Transfer from 1x1 to 2x1") {
  if(scalapack::global_size() < 2)
    return;
  mpi::Communicator const world;
  auto const ranks = world.all_gather(scalapack::global_rank());
  auto const root = std::find(ranks.begin(), ranks.end(), 0) - ranks.begin();
  scalapack::Context const single(1, 1);
  scalapack::Context const parallel(2, 1);
  scalapack::Context const all(world.size(), 1);

  auto const split = mpi::Communicator().split(single.is_valid() or parallel.is_valid());
  scalapack::Sizes const size = {64, 128};
  scalapack::Sizes const blocks = {16, 32};
  scalapack::Matrix<> input(single, size, blocks);
  if(single.is_valid())
    input.local() = optimet::Matrix<t_real>::Random(size.rows, size.cols);
  auto const input_matrix = split.broadcast(input.local(), root);

  // perform transfer
  auto const scattered = input.transfer_to(all, parallel);
  // Global size remain the same
  CHECK(scattered.rows() == input.rows());
  CHECK(scattered.cols() == input.cols());
  if(parallel.is_valid()) {
    // Check local sizes
    CHECK(static_cast<t_uint>(scattered.local().cols()) == size.cols);
    CHECK(static_cast<t_uint>(scattered.local().rows()) == size.rows / 2);
    // Check upper left block remains the same
    // Other blocks are too painful to check
    if(parallel.row() == 0 and parallel.col() == 0)
      CHECK(input_matrix.topLeftCorner(blocks.rows, blocks.cols)
                .isApprox(scattered.local().topLeftCorner(blocks.rows, blocks.cols)));
  }

  // gather matrix back to origin
  auto const gathered = scattered.transfer_to(single);
  // Global size remain the same
  CHECK(gathered.rows() == input.rows());
  CHECK(gathered.cols() == input.cols());
  if(single.is_valid()) {
    // Check local sizes
    CHECK(static_cast<t_uint>(gathered.local().cols()) == size.cols);
    CHECK(static_cast<t_uint>(gathered.local().rows()) == size.rows);
    // Check input matrix was recovered
    CHECK(input_matrix.isApprox(gathered.local(), 1e-12));
  }
}

void check_distribute(optimet::scalapack::Sizes const &grid, optimet::scalapack::Sizes const &size,
                      optimet::scalapack::Sizes const &blocks) {
  mpi::Communicator const world;
  if(world.size() < grid.rows * grid.cols) {
    WARN("Not enough processes to run test: " << grid.rows << "x" << grid.cols << " < "
                                              << world.size());
    return;
  }

  auto const ranks = world.all_gather(scalapack::global_rank());
  auto const root = std::find(ranks.begin(), ranks.end(), 0) - ranks.begin();
  scalapack::Context const single(1, 1);
  scalapack::Context const parallel(grid.rows, grid.cols);
  scalapack::Context const transpose(grid.cols, grid.rows);
  scalapack::Context const all(world.size(), 1);

  auto const split = mpi::Communicator().split(single.is_valid() or parallel.is_valid());
  scalapack::Matrix<> input(single, size, blocks);
  if(single.is_valid())
    input.local() = optimet::Matrix<t_real>::Random(size.rows, size.cols);
  auto const input_matrix = split.broadcast(input.local(), root);

  // perform transfer
  auto scattered = input.transfer_to(all, parallel);
  scattered.local() *= 2;
  // go to transpose
  auto transposed = scattered.transfer_to(transpose);
  transposed.local() *= 1.5;
  // return to single
  auto const gathered = transposed.transfer_to(single);
  // Global size remain the same
  CHECK(gathered.rows() == input.rows());
  CHECK(gathered.cols() == input.cols());
  if(single.is_valid())
    CHECK((3e0 * input_matrix).isApprox(gathered.local(), 1e-12));
}

TEST_CASE("Transfer from 1x1 to nxm to mxn to 1x1") {
  SECTION("1x2") { check_distribute({1, 2}, {1024, 2048}, {32, 65}); }
  SECTION("2x1") { check_distribute({1, 2}, {1024, 2048}, {32, 65}); }
  SECTION("3x1") { check_distribute({3, 1}, {1024, 2048}, {32, 65}); }
  SECTION("2x2") { check_distribute({2, 2}, {1024, 2048}, {32, 65}); }
  SECTION("3x2") { check_distribute({3, 2}, {1024, 2048}, {32, 65}); }
}

void check_assignement(optimet::scalapack::Sizes const &grid, optimet::scalapack::Sizes const &size,
                       optimet::scalapack::Sizes const &blocks) {
  mpi::Communicator const world;
  if(world.size() < grid.rows * grid.cols) {
    WARN("Not enough processes to run test: " << grid.rows << "x" << grid.cols << " < "
                                              << world.size());
    return;
  }

  auto const ranks = world.all_gather(scalapack::global_rank());
  auto const root = std::find(ranks.begin(), ranks.end(), 0) - ranks.begin();
  scalapack::Context const single(1, 1);
  scalapack::Context const parallel(grid.rows, grid.cols);

  auto const split = mpi::Communicator().split(single.is_valid() or parallel.is_valid());
  scalapack::Matrix<> input(single, size, blocks);
  if(single.is_valid())
    input.local() = optimet::Matrix<t_real>::Random(size.rows, size.cols);
  auto const input_matrix = split.broadcast(input.local(), root);

  // Perform transfer to parallel matrix
  scalapack::Matrix<> mid(parallel, size, blocks);
  scalapack::Matrix<> output(single, size, blocks);
  mid = input;
  output = mid;

  if(single.is_valid())
    CHECK(output.local().isApprox(input.local()));

  output.local().fill(0);
  output = input;
  if(single.is_valid())
    CHECK(output.local().isApprox(input.local()));
}

TEST_CASE("Assignement") {
  SECTION("1x1") { check_assignement({1, 1}, {1024, 2048}, {32, 65}); }
  SECTION("1x2") { check_assignement({1, 2}, {1024, 2048}, {32, 65}); }
  SECTION("2x1") { check_assignement({2, 1}, {1024, 2048}, {32, 65}); }
  SECTION("3x1") { check_assignement({3, 1}, {1024, 2048}, {32, 65}); }
  SECTION("2x2") { check_assignement({2, 2}, {1024, 2048}, {32, 65}); }
  SECTION("3x2") { check_assignement({3, 2}, {1024, 2048}, {32, 65}); }
}

void check_multiply(optimet::scalapack::Sizes const &grid, optimet::scalapack::Sizes const &size,
                    optimet::scalapack::Sizes const &blocks) {
  mpi::Communicator const world;
  if(world.size() < grid.rows * grid.cols) {
    WARN("Not enough processes to run test: " << grid.rows << "x" << grid.cols << " < "
                                              << world.size());
    return;
  }

  auto const ranks = world.all_gather(scalapack::global_rank());
  auto const root = std::find(ranks.begin(), ranks.end(), 0) - ranks.begin();
  scalapack::Context const single(1, 1);
  scalapack::Context const matrix_context(grid.rows, grid.cols);
  scalapack::Context const vector_context(grid.rows, 1);

  auto const condition =
      single.is_valid() or matrix_context.is_valid() or vector_context.is_valid();
  auto const split = mpi::Communicator().split(condition);
  if(not condition)
    return;
  scalapack::Matrix<> Aserial(single, size, blocks);
  scalapack::Matrix<> xserial(single, {Aserial.cols(), static_cast<t_uint>(1)}, blocks);
  if(single.is_valid()) {
    Aserial.local() = optimet::Matrix<t_real>::Random(Aserial.rows(), Aserial.cols());
    xserial.local() = optimet::Matrix<t_real>::Random(xserial.rows(), xserial.cols());
  }

  // perform transfer
  auto Aparallel = Aserial.transfer_to(matrix_context);
  auto xparallel = xserial.transfer_to(vector_context);
  // Do multiplication
  scalapack::Matrix<> result(single, {Aserial.rows(), 1}, blocks);
  result.local().fill(0);
  pdgemm(1.5e0, Aparallel, xparallel, 0.0, result);

  if(single.is_valid()) {
    auto const expected = (1.5 * Aserial.local() * xserial.local()).eval();
    CHECK(expected.isApprox(result.local()));
  }
}

TEST_CASE("Matrix multiplication") {
  SECTION("1x1") { check_multiply({1, 1}, {1024, 2048}, {32, 65}); }
  SECTION("1x2") { check_multiply({1, 2}, {1024, 2048}, {32, 65}); }
  SECTION("2x1") { check_multiply({2, 1}, {1024, 2048}, {32, 65}); }
  SECTION("2x3") { check_multiply({2, 3}, {1024, 2048}, {32, 65}); }
  SECTION("matrix smaller than block size in one dimensions") {
    check_multiply({2, 2}, {64, 32}, {64, 64});
  }
  SECTION("matrix smaller than block size in both dimensions") {
    check_multiply({2, 2}, {32, 32}, {64, 64});
  }
}

TEST_CASE("Empty matrix") {
  scalapack::Matrix<> A(scalapack::Context(), {0, 0}, {64, 64});
  CHECK(A.rows() == 0);
  CHECK(A.cols() == 0);
}
