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

TEST_CASE("Broadcasting") {
  mpi::Communicator const world;

  SECTION("From root") {
    CHECK(world.broadcast(world.rank() * 2) == 0);
    CHECK(world.broadcast(world.rank() * 2 + 1) == 1);
    CHECK(world.broadcast(static_cast<double>(world.rank() * 2) + 1.5) == 1.5);

    auto const value = world.is_root() ?
      world.broadcast('c'): world.broadcast<char>();
    CHECK(value == 'c');
  }

  SECTION("From other") {
    if(world.size() == 1)
      return;
    // Test written expecting root is zero;
    REQUIRE(world.root_id() == 0);
    auto const root = 1u;

    CHECK(world.broadcast(world.rank() * 2u, root) == 2u);
    CHECK(world.broadcast(world.rank() * 2u + 1u, root) == 3u);
    CHECK(world.broadcast(static_cast<double>(world.rank() * 2u) + 1.5, root) == 3.5);

    auto const value = world.is_root() ?
      world.broadcast('c', root):
      world.rank() == 1 ?
      world.broadcast('d', root): world.broadcast<char>(root);
    CHECK(value == 'd');
  }
}

TEST_CASE("Gathering") {
  mpi::Communicator const world;

  SECTION("From root") {
    auto const result = world.gather(world.rank() * 2);
    for(std::size_t i(0); i < result.size(); ++i)
      CHECK(result[i] == i * 2);
  }
}
