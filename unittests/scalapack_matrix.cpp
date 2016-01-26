#include <iostream>
#include <numeric>
#include "catch.hpp"

#include "Types.h"
#include "scalapack/Context.h"
#include "scalapack/InitExit.h"
#include "scalapack/Matrix.h"
#include "mpi/Communicator.h"
#include "mpi/Collectives.h"

using namespace optimet;

TEST_CASE("Creates a matrix in 1x1 context") {
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(1, 1);
  CAPTURE(context.is_valid());
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
  CHECK(matrix.blacs()[8] == leading);
}

TEST_CASE("Creates a matrix in 1x2 context") {
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(1, 2);
  CAPTURE(context.is_valid());
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
  CHECK(matrix.blacs()[8] == leading);
}

void check_creation(t_uint n, t_uint m) {
  if(scalapack::global_size() < n * m)
    return;
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(n, m);
  CAPTURE(context.is_valid());
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
  CHECK(matrix.blacs()[8] == leading);
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

void check_distribute(optimet::scalapack::Sizes const &grid,
                      optimet::scalapack::Sizes const &size,
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
  SECTION("3x1") { check_distribute({3, 1}, {1024, 2048}, {32, 65}); }
  SECTION("2x2") { check_distribute({2, 2}, {1024, 2048}, {32, 65}); }
  SECTION("3x2") { check_distribute({3, 2}, {1024, 2048}, {32, 65}); }
}
