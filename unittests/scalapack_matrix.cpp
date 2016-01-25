#include <iostream>
#include <numeric>
#include "catch.hpp"

#include "Types.h"
#include "scalapack/Context.h"
#include "scalapack/InitExit.h"
#include "scalapack/Matrix.h"
#include "MpiCommunicator.h"

using namespace optimet;

TEST_CASE("Creates a matrix in 1x1 context") {
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(1, 1);
  CAPTURE(context.is_valid());
  REQUIRE(context.is_valid() == (rank == 0));
  scalapack::Matrix matrix(context, {64, 64}, {20, 10});

  CHECK(matrix.eigen().rows() == (context.is_valid() ? 64 : 0));
  CHECK(matrix.eigen().cols() == (context.is_valid() ? 64 : 0));
  CHECK(matrix.blacs()[0] == 1);
  if(context.is_valid())
    CHECK(matrix.blacs()[1] == *context);
  CHECK(matrix.blacs()[2] == 64);
  CHECK(matrix.blacs()[3] == 64);
  CHECK(matrix.blacs()[4] == 20);
  CHECK(matrix.blacs()[5] == 10);
  CHECK(matrix.blacs()[6] == 0);
  CHECK(matrix.blacs()[7] == 0);
  auto const leading = matrix.eigen().IsRowMajor ? matrix.eigen().rows() : matrix.eigen().cols();
  CHECK(matrix.blacs()[8] == leading);
}

void check(t_uint n, t_uint m) {
  if(scalapack::global_size() < n * m)
    return;
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(n, m);
  CAPTURE(context.is_valid());
  scalapack::Matrix::Sizes const size = {1024, 1024};
  scalapack::Matrix::Sizes const blocks = {31, 65};
  scalapack::Matrix matrix(context, size, blocks);

  if(context.is_valid()) {
    REQUIRE(matrix.eigen().rows() > 0);
    REQUIRE(matrix.eigen().cols() > 0);
  } else {
    REQUIRE(matrix.eigen().rows() == 0);
    REQUIRE(matrix.eigen().cols() == 0);
  }

  auto const split = mpi::Communicator().split(context.is_valid());
  if(not context.is_valid())
    return;

  REQUIRE(split.size() == n * m);
  auto const rows = split.all_gather(matrix.eigen().rows());
  auto const cols = split.all_gather(matrix.eigen().cols());
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
  auto const leading = matrix.eigen().IsRowMajor ? matrix.eigen().rows() : matrix.eigen().cols();
  CHECK(matrix.blacs()[8] == leading);
}

TEST_CASE("Creates a matrice in nxm context") {
  check(1, 1);
  check(2, 1);
  check(1, 2);
  check(2, 2);
  check(3, 1);
  check(3, 2);
}
