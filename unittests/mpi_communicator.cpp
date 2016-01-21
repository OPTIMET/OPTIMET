#include <iostream>
#include "catch.hpp"

#include "Types.h"
#include "MpiCommunicator.h"
#include "MpiExit.h"

using namespace optimet;

TEST_CASE("Creates an mpi communicator") {
  CHECK(mpi_initialized());
  int rank, size;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_rank(MPI_COMM_WORLD, &size);

  MpiCommunicator const world;

  CHECK(*world == MPI_COMM_WORLD);
  CHECK(world.rank() == rank);
  CHECK(world.size() == size);
}
