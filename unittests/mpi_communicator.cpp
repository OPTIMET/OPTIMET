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

#include "Types.h"
#include "mpi/Collectives.h"
#include "mpi/Communicator.h"
#include "mpi/GraphCommunicator.h"
#include "mpi/Session.h"

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
    auto const split = world.split(world.is_root() ? 0 : 1);
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
  if(world.size() == 1)
    return;

  SECTION("From root") {
    CHECK(world.broadcast(world.rank() * 2) == 0);
    CHECK(world.broadcast(world.rank() * 2 + 1) == 1);
    CHECK(world.broadcast(static_cast<double>(world.rank() * 2) + 1.5) == 1.5);

    auto const value = world.is_root() ? world.broadcast('c') : world.broadcast<char>();
    CHECK(value == 'c');
  }

  SECTION("From other") {
    // Test written expecting root is zero;
    REQUIRE(world.root_id() == 0);
    auto const root = 1u;

    CHECK(world.broadcast(world.rank() * 2u, root) == 2u);
    CHECK(world.broadcast(world.rank() * 2u + 1u, root) == 3u);
    CHECK(world.broadcast(static_cast<double>(world.rank() * 2u) + 1.5, root) == 3.5);

    auto const value = world.is_root() ? world.broadcast('c', root) : world.rank() == 1 ?
                                         world.broadcast('d', root) :
                                         world.broadcast<char>(root);
    CHECK(value == 'd');
  }

  SECTION("Matrix") {
    Matrix<t_real> input(2, 3);
    for(Matrix<>::Index i(0); i < input.rows(); ++i)
      for(Matrix<>::Index j(0); j < input.cols(); ++j)
        input(i, j) = 2 * i + j;
    Matrix<t_real> matrix =
        world.is_root() ? input : Matrix<t_real>::Zero(input.rows(), input.cols());
    auto const result = world.broadcast(matrix);
    CHECK(result.isApprox(input, 1e-12));
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

TEST_CASE("Symmetric graph communicators") {
  mpi::Communicator world;
  if(world.size() < 3)
    return;

  std::vector<std::set<t_uint>> const comms =
      mpi::GraphCommunicator::symmetrize({{1, 2}, {0, 2}, {0, 1}, {}});

  auto const rank = std::min<t_uint>(world.rank(), 3);
  mpi::GraphCommunicator graph(world, comms);
  CHECK(graph.neighborhood_size(rank) == comms[rank].size());
  for(decltype(comms)::size_type i(0); i < comms.size(); ++i) {
    REQUIRE(graph.neighborhood(i).size() == comms[i].size());
    for(decltype(comms)::value_type::size_type j(0); j < comms[i].size(); ++j) {
      auto iter = comms[i].begin();
      std::advance(iter, j);
      CHECK(graph.neighborhood(i)[j] == *iter);
    }
  }
}

TEST_CASE("Blocking gather of scalar on graph") {
  mpi::Communicator world;
  if(world.size() < 3)
    return;

  std::vector<std::set<t_uint>> const comms =
      mpi::GraphCommunicator::symmetrize({{2}, {2}, {1, 0}, {}});
  std::vector<int> const values = {2, 4, 1, 3};

  auto const rank = std::min<t_uint>(world.rank(), 3);
  mpi::GraphCommunicator graph(world, comms);
  auto const actual = graph.allgather(values[rank]);

  CHECK(actual.size() == comms[rank].size());
  for(decltype(actual)::size_type i(0); i < actual.size(); ++i) {
    auto iter = comms[rank].begin();
    std::advance(iter, i);
    CHECK(actual[i] == values[*iter]);
  }
}

TEST_CASE("Additive symmetrization") {
  std::vector<std::set<t_uint>> const expected = {{2}, {2}, {1, 0}, {}};
  std::vector<std::set<t_uint>> const actual =
      mpi::GraphCommunicator::symmetrize({{}, {2}, {0}, {}});

  CHECK(expected.size() == actual.size());
  for(decltype(actual)::size_type i(0); i < actual.size(); ++i) {
    CHECK(actual[i].size() == expected[i].size());
    for(auto const node : expected[i])
      CHECK(actual[i].count(node) == 1);
  }
}

TEST_CASE("Non-blocking gather of Eigen vectors on graph") {
  mpi::Communicator world;
  if(world.size() < 3)
    return;

  auto const size = [](t_int rank) { return 3 * (rank + 1); };
  std::vector<std::set<t_uint>> const comms =
      mpi::GraphCommunicator::symmetrize({{2}, {2}, {1, 0}, {}});
  std::vector<int> const values = {3, 5, 1, 0};

  auto const rank = std::min<t_uint>(world.rank(), 3);
  mpi::GraphCommunicator graph(world, comms);

  auto const receive_count = graph.allgather(size(rank));
  Vector<int> input = Vector<int>::Constant(size(rank), values[rank]);
  Vector<int> result;
  if(auto request = graph.iallgather(input, result, receive_count)) {
    auto const N = std::accumulate(receive_count.begin(), receive_count.end(), 0u);
    CHECK(result.size() == N);
    CHECK(static_cast<bool>(request));
  }

  for(t_uint i(0), j(0); i < receive_count.size(); j += receive_count[i++]) {
    auto iter = comms[rank].begin();
    std::advance(iter, i);
    CHECK((result.segment(j, receive_count[i]).array() == values[*iter]).all());
  }
}

TEST_CASE("Non-blocking all-2-all of Eigen vectors on graph") {
  mpi::Communicator world;
  if(world.size() < 3)
    return;

  auto const vecsize = [](int a, int b) { return a == b ? 0 : 3 * a + b; };
  auto const values = [](int a, int b) { return 2 * a + 3 * b + 1; };
  std::vector<std::set<t_uint>> const comms =
      mpi::GraphCommunicator::symmetrize({{2}, {2}, {1, 0}, {}});

  auto const rank = std::min<t_uint>(world.rank(), 3);
  mpi::GraphCommunicator graph(world, comms);

  std::vector<int> receive_counts;
  std::transform(comms[rank].begin(), comms[rank].end(), std::back_inserter(receive_counts),
                 [vecsize, rank](int a) { return vecsize(rank, a); });
  std::vector<int> send_counts;
  std::transform(comms[rank].begin(), comms[rank].end(), std::back_inserter(send_counts),
                 [vecsize, rank](int a) { return vecsize(a, rank); });

  // construct message to send from this proc
  auto const Nsend = std::accumulate(send_counts.begin(), send_counts.end(), 0);
  Vector<int> message(Nsend);
  for(int i(0), j(0); i < comms[rank].size(); j += send_counts[i++]) {
    auto iter = comms[rank].begin();
    std::advance(iter, i);
    message.segment(j, send_counts[i]).fill(values(*iter, rank));
  }

  // construct message we expect to *receive* on this proc
  auto const Nreceive = std::accumulate(receive_counts.begin(), receive_counts.end(), 0u);
  Vector<int> expected(Nreceive);
  for(int i(0), j(0); i < comms[rank].size(); j += receive_counts[i++]) {
    auto iter = comms[rank].begin();
    std::advance(iter, i);
    expected.segment(j, receive_counts[i]).fill(values(rank, *iter));
  }

  // Now does all to all over graph
  Vector<int> actual;
  if(auto request = graph.ialltoall(message, actual, send_counts, receive_counts)) {
    CHECK(actual.size() == Nreceive);
    CHECK(static_cast<bool>(request));
  }

  CHECK(actual.transpose() == expected.transpose());
}
