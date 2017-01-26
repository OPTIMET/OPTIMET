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

#include "mpi/GraphCommunicator.h"
#include <algorithm>
#include <exception>
#include <iterator>
#include <mpi.h>

namespace optimet {
namespace mpi {

GraphCommunicator::GraphCommunicator(Communicator const &comm,
                                     std::vector<std::set<t_uint>> const &graph, bool reorder) {

  Communicator::reset(nullptr);
  std::vector<int> degrees;
  std::vector<int> edges;
  for(t_uint i(0); i < std::min<t_uint>(comm.size(), graph.size()); ++i) {
    degrees.push_back((degrees.size() > 0 ? degrees.back() : 0) + graph[i].size());
    std::copy(graph[i].begin(), graph[i].end(), std::back_inserter(edges));
  }
  if(edges.size() == 0)
    return;
  if(*std::max_element(edges.begin(), edges.end()) >= comm.size())
    throw std::out_of_range("Edge index in graph communicator is larger than parent comm size");

  if(degrees.back() == 0)
    return;
  degrees.resize(comm.size(), degrees.back());
  MPI_Comm result;
  auto const error =
      MPI_Graph_create(*comm, degrees.size(), degrees.data(), edges.data(), reorder, &result);
  if(error != MPI_SUCCESS)
    throw std::runtime_error("Could not create graph communicator");
  Communicator::reset(&result);
}

t_uint GraphCommunicator::neighborhood_size(int rank) const {
  if(not is_valid())
    return 0;

  int result;
  auto const error = MPI_Graph_neighbors_count(**this, rank, &result);
  if(error != MPI_SUCCESS)
    throw std::runtime_error("Could not retrieve the number of input edges");
  return static_cast<t_uint>(result);
}

std::vector<t_uint> GraphCommunicator::neighborhood(int rank) const {
  if(not is_valid())
    return std::vector<t_uint>();

  int const N = neighborhood_size(rank);
  if(N == 0)
    return std::vector<t_uint>();
  std::vector<int> result(N);
  auto const error = MPI_Graph_neighbors(**this, rank, N, result.data());
  if(error != MPI_SUCCESS)
    throw std::runtime_error("Could not retrieve the number of input edges");
  return std::vector<t_uint>(result.begin(), result.end());
}

namespace {
void wait_on_delete(MPI_Request *request) {
  if(request == nullptr)
    return;

  MPI_Status result;
  auto const error = MPI_Wait(request, &result);
  if(error != MPI_SUCCESS)
    throw std::runtime_error("Got an error when waiting for request to complete");
  delete request;
}
}
Request mpi_request_wait_on_delete(MPI_Request *const request) {
  return Request(request, wait_on_delete);
}

std::vector<std::set<t_uint>>
GraphCommunicator::symmetrize(std::vector<std::set<t_uint>> const &graph) {
  std::vector<std::set<t_uint>> result(graph);
  for(decltype(result)::size_type i(0); i < result.size(); ++i)
    for(auto const node : result[i]) {
      if(node > result.size())
        result.resize(node + 1);
      result[node].insert(i);
    }
  return result;
}

} /* optimet::mpi */
} /* optimet  */
