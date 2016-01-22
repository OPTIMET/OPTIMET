#include <iostream>
#include "catch.hpp"

#include "Types.h"
#include "ScalapackContext.h"
#include "ScalapackExit.h"
#include "MpiCommunicator.h"

using namespace optimet;

TEST_CASE("Check pinfo") {
  mpi::Communicator world;
  CHECK(world.size() == scalapack::global_size());

  auto ranks = world.gather(scalapack::global_rank());
  if(not world.is_root())
    return;

  REQUIRE(ranks.size() == world.size());
  std::sort(ranks.begin(), ranks.end(), std::less<t_uint>());
  for(std::size_t i(0); i < ranks.size(); ++i)
    CHECK(i == ranks[i]);
}

TEST_CASE("Creates a blacs context 1x1") {
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(1, 1);
  REQUIRE(context.is_valid() == (rank == 0));
  if(rank == 0) {
    CHECK(context.rows() == 1);
    CHECK(context.cols() == 1);
    CHECK(context.row() == 0);
    CHECK(context.col() == 0);
  }
}

void check_nxm(t_uint n, t_uint m) {
  if(mpi::Communicator().size() < m * n)
    return;
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(n, m);
  REQUIRE(context.is_valid() == (rank < n * m));
  if(rank < n * m) {
    CHECK(context.rows() == n);
    CHECK(context.cols() == m);
    // Fortran is column major
    CHECK(context.col() == rank / n);
    CHECK(context.row() == rank % n);
  }
}

TEST_CASE("Create different blacs context") {
  check_nxm(1, 1);
  check_nxm(2, 1);
  check_nxm(1, 2);
  check_nxm(2, 2);
  check_nxm(3, 2);
}
