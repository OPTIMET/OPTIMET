#include "catch.hpp"
#include "mpi/FastMatrixMultiply.h"
#include <iostream>

TEST_CASE("ReduceComputation") {
  using namespace optimet;
  mpi::Communicator const world;
  if(world.size() < 3)
    return;

  ElectroMagnetic const silicon{13.1, 1.0};
  auto const nprocs = 3;
  auto const nHarmonics = 5;
  auto const radius = 500e-9;
  std::vector<Scatterer> const scatterers = {
      {{0, 0, 0}, silicon, radius, nHarmonics},     {{0, 0, 0}, silicon, radius, nHarmonics + 1},
      {{0, 0, 0}, silicon, radius, nHarmonics + 2}, {{0, 0, 0}, silicon, radius, nHarmonics + 3},
      {{0, 0, 0}, silicon, radius, nHarmonics + 4}, {{0, 0, 0}, silicon, radius, nHarmonics + 5},
  };

  auto const distribution = mpi::details::vector_distribution(scatterers.size(), nprocs);
  Matrix<bool> locals = Matrix<bool>::Zero(scatterers.size(), scatterers.size());
  locals.topLeftCorner(3, 3).fill(true);
  locals.bottomRightCorner(3, 3).fill(true);
  mpi::GraphCommunicator const graph_comm(
      world, mpi::GraphCommunicator::symmetrize(
                 mpi::details::non_local_graph_edges(locals.array() == false, distribution)),
      false);

  auto const nfunctions = [&scatterers](t_int i) {
    return 2 * scatterers[i].nMax * (scatterers[i].nMax + 2);
  };
  int const sizes[] = {nfunctions(0) + nfunctions(1) + nfunctions(2),
                       nfunctions(3) + nfunctions(4) + nfunctions(5)};
  Vector<int> const messages[] = {
      Vector<int>::LinSpaced(sizes[1], sizes[0], sizes[1]).eval(),
      10 * Vector<int>::LinSpaced(sizes[0] + sizes[1], 0, sizes[0] + sizes[1]).eval(),
      100 * Vector<int>::LinSpaced(sizes[0], 0, sizes[0]).eval(), Vector<int>::Zero(0).eval()};

  mpi::FastMatrixMultiply::ReduceComputation reduction(graph_comm, locals.array() == false,
                                                       distribution, scatterers);
  std::vector<int> const buffer_sizes{2 * (nfunctions(0) + nfunctions(1)),
                                      2 * (nfunctions(2) + nfunctions(3)),
                                      2 * (nfunctions(4) + nfunctions(5)), 0};
  SECTION("Send message") {
    Vector<int> buffer;
    if(auto const request = reduction.send(messages[std::min<int>(world.rank(), nprocs)], buffer))
      REQUIRE(request);
    else
      REQUIRE(false);

    switch(world.rank()) {
    case 0:
      REQUIRE(buffer.size() == 2 * (nfunctions(0) + nfunctions(1)));
      CHECK(buffer.head(nfunctions(0) + nfunctions(1)) ==
            messages[1].head(nfunctions(0) + nfunctions(1)));
      CHECK(buffer.tail(nfunctions(0) + nfunctions(1)) ==
            messages[2].head(nfunctions(0) + nfunctions(1)));
      break;
    case 1:
      REQUIRE(buffer.size() == 2 * (nfunctions(2) + nfunctions(3)));
      CHECK(buffer.head(nfunctions(3)) == messages[0].head(nfunctions(3)));
      CHECK(buffer.segment(nfunctions(3), nfunctions(2) + nfunctions(3)) ==
            messages[1].segment(nfunctions(0) + nfunctions(1), nfunctions(2) + nfunctions(3)));
      CHECK(buffer.tail(nfunctions(2)) == messages[2].tail(nfunctions(2)));
      break;
    case 2:
      REQUIRE(buffer.size() == 2 * (nfunctions(4) + nfunctions(5)));
      CHECK(buffer.head(nfunctions(4) + nfunctions(5)) ==
            messages[0].tail(nfunctions(4) + nfunctions(5)));
      CHECK(buffer.tail(nfunctions(4) + nfunctions(5)) ==
            messages[1].tail(nfunctions(4) + nfunctions(5)));
      break;

    default:
      REQUIRE(buffer.size() == 0);
    }
  }

  SECTION("Do Reduction") {
    std::vector<int> const result_sizes{nfunctions(0) + nfunctions(1),
                                        nfunctions(2) + nfunctions(3),
                                        nfunctions(4) + nfunctions(5), 0};
    auto const rank = std::min<int>(world.rank(), nprocs);
    Vector<int> const buffer =
        (world.rank() + 2) * Vector<int>::LinSpaced(buffer_sizes[rank], 0, buffer_sizes[rank]);
    Vector<int> result(result_sizes[rank]);
    result.fill(world.rank() + 1);
    reduction.reduce(result, buffer);

    switch(world.rank()) {
    case 0:
      CHECK(result ==
            result.Constant(result.size(), world.rank() + 1) +
                buffer.head(nfunctions(0) + nfunctions(1)) +
                buffer.tail(nfunctions(0) + nfunctions(1)));
      break;
    case 1:
      CHECK(result.head(nfunctions(2)) ==
            result.Constant(nfunctions(2), world.rank() + 1) +
                buffer.segment(nfunctions(3), nfunctions(2)) + buffer.tail(nfunctions(2)));
      CHECK(result.tail(nfunctions(3)) ==
            result.Constant(nfunctions(3), world.rank() + 1) +
                buffer.segment(nfunctions(3) + nfunctions(2), nfunctions(3)) +
                buffer.head(nfunctions(3)));
      break;
    case 2:
      CHECK(result ==
            result.Constant(result.size(), world.rank() + 1) +
                buffer.head(nfunctions(4) + nfunctions(5)) +
                buffer.tail(nfunctions(4) + nfunctions(5)));
      break;
    default:
      REQUIRE(result.size() == 0);
    }
  }
}

TEST_CASE("DistributeInput") {
  using namespace optimet;
  mpi::Communicator const world;
  if(world.size() < 3)
    return;

  ElectroMagnetic const silicon{13.1, 1.0};
  auto const nprocs = 3;
  auto const nHarmonics = 5;
  auto const radius = 500e-9;
  std::vector<Scatterer> const scatterers = {
      {{0, 0, 0}, silicon, radius, nHarmonics},     {{0, 0, 0}, silicon, radius, nHarmonics + 1},
      {{0, 0, 0}, silicon, radius, nHarmonics + 2}, {{0, 0, 0}, silicon, radius, nHarmonics + 3},
      {{0, 0, 0}, silicon, radius, nHarmonics + 4}, {{0, 0, 0}, silicon, radius, nHarmonics + 5},
  };

  auto const vector_distribution = mpi::details::vector_distribution(scatterers.size(), nprocs);
  Matrix<bool> locals = Matrix<bool>::Zero(scatterers.size(), scatterers.size());
  locals.topLeftCorner(3, 3).fill(true);
  locals.bottomRightCorner(3, 3).fill(true);
  mpi::GraphCommunicator const graph_comm(
      world, mpi::GraphCommunicator::symmetrize(
                 mpi::details::non_local_graph_edges(locals.array() == false, vector_distribution)),
      false);

  auto const nfunctions = [&scatterers](t_int i) {
    return 2 * scatterers[i].nMax * (scatterers[i].nMax + 2);
  };
  int const sizes[] = {nfunctions(0) + nfunctions(1), nfunctions(2) + nfunctions(3),
                       nfunctions(4) + nfunctions(5), 0};
  Vector<int> const messages[] = {Vector<int>::LinSpaced(sizes[0], 0, sizes[0]).eval(),
                                  10 * Vector<int>::LinSpaced(sizes[1], 0, sizes[1]).eval(),
                                  100 * Vector<int>::LinSpaced(sizes[2], 0, sizes[2]).eval(),
                                  Vector<int>::Zero(0).eval()};

  mpi::FastMatrixMultiply::DistributeInput distribution(graph_comm, locals.array() == false,
                                                        vector_distribution, scatterers);
  std::vector<int> const buffer_sizes{sizes[1] + sizes[2], sizes[0] + sizes[1] + sizes[2],
                                      sizes[0] + sizes[1], 0};
  auto const rank = std::min<int>(world.rank(), nprocs);
  SECTION("Send message") {
    Vector<int> buffer;
    if(auto const request = distribution.send(messages[rank], buffer))
      REQUIRE(request);
    else
      REQUIRE(false);

    REQUIRE(buffer.size() == buffer_sizes[rank]);
    switch(world.rank()) {
    case 0:
      CHECK(buffer.head(sizes[1]) == messages[1]);
      CHECK(buffer.tail(sizes[2]) == messages[2]);
      break;
    case 1:
      CHECK(buffer.head(sizes[0]) == messages[0]);
      CHECK(buffer.tail(sizes[2]) == messages[2]);
      break;
    case 2:
      CHECK(buffer.head(sizes[0]) == messages[0]);
      CHECK(buffer.tail(sizes[1]) == messages[1]);
      break;
    default:
      break;
    }
  }

  SECTION("Synthesize output") {
    Vector<int> const buffer =
        Vector<int>::LinSpaced(buffer_sizes[rank], 1, buffer_sizes[rank] + 1);
    Vector<int> const result;
    distribution.synthesize(buffer, result);
    switch(world.rank()) {
    case 0:
      REQUIRE(result.size() == nfunctions(3) + nfunctions(4) + nfunctions(5));
      CHECK(result == buffer.tail(result.size()));
      break;
    case 1:
      REQUIRE(result.size() ==
              nfunctions(0) + nfunctions(1) + nfunctions(2) + nfunctions(3) + nfunctions(4) +
                  nfunctions(5));
      CHECK(result == buffer);
      break;
    case 2:
      REQUIRE(result.size() == nfunctions(0) + nfunctions(1) + nfunctions(2));
      CHECK(result == buffer.head(nfunctions(0) + nfunctions(1) + nfunctions(2)));
      break;
    default:
      REQUIRE(result.size() == 0);
      break;
    }
  }
}
