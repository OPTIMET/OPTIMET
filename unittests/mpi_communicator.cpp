#include <iostream>
#include "catch.hpp"

#include "Types.h"
#include "MpiCommunicator.h"
#include "MpiExit.h"

using namespace optimet;

TEST_CASE("Creates an mpi communicator") {
  CHECK(mpi::initialized());
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  mpi::Communicator const world;

  SECTION("General stuff") {
    REQUIRE(*world == MPI_COMM_WORLD);
    REQUIRE(static_cast<t_int>(world.rank()) == rank);
    REQUIRE(static_cast<t_int>(world.size()) == size);

    mpi::Communicator shallow = world;
    CHECK(*shallow == *world);
  }

  SECTION("Split") {
    auto const split = world.split(world.is_root() ? 0: 1);
    if(world.is_root())
      CHECK(split.size() == 1);
    else {
      CHECK(split.size() == world.size() - 1);
      CHECK(split.rank() == world.rank() - 1);
    }
  }

  SECTION("Duplicate") {
    mpi::Communicator dup = world.duplicate();
    CHECK(*dup != *world);
  }
}
