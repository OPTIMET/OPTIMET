#include "mpi/GraphCommunicator.h"
#include <algorithm>
#include <exception>
#include <iterator>
#include <mpi.h>

namespace optimet {
namespace mpi {

GraphCommunicator::GraphCommunicator(Communicator const &comm,
                                     std::vector<std::set<t_uint>> const &graph, bool reorder) {

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

namespace {
void wait_on_delete(MPI_Request *request) {
  if(request == nullptr)
    return;

  MPI_Status result;
  auto const error = MPI_Wait(request, &result);
  if(error != MPI_SUCCESS)
    throw std::runtime_error("Got an error when waiting for request to complete");
}
}
Request make_wait_on_delete(MPI_Request *const request) {
  return Request(request, &wait_on_delete);
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
