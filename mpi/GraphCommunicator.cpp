#include "mpi/GraphCommunicator.h"
#include <algorithm>
#include <exception>
#include <iterator>
#include <mpi.h>

namespace optimet {
namespace mpi {

DistGraphCommunicator::DistGraphCommunicator(Communicator const &comm,
                                             std::vector<int> const &sources,
                                             std::vector<int> const &destinations, bool reorder) {
  if(sources.size() != 0) {
    if(*std::max_element(sources.begin(), sources.end()) >= comm.size())
      throw std::out_of_range("Sources in graph communicator is larger than parent comm size");
    if(sources.size() == 0 || *std::min_element(sources.begin(), sources.end()) < 0)
      throw std::out_of_range("Sources in graph communicator is smaller than 0");
  }
  if(destinations.size() != 0) {
    if(*std::max_element(destinations.begin(), destinations.end()) >= comm.size())
      throw std::out_of_range("Destination in graph communicator is larger than parent comm size");
    if(*std::min_element(destinations.begin(), destinations.end()) < 0)
      throw std::out_of_range("Destination in graph communicator is smaller than 0");
  }

  MPI_Comm result;
  auto const error = MPI_Dist_graph_create_adjacent(
      *comm, sources.size(), sources.data(), MPI_UNWEIGHTED, destinations.size(),
      destinations.data(), MPI_UNWEIGHTED, MPI_INFO_NULL, reorder, &result);
  if(error != MPI_SUCCESS)
    throw std::runtime_error("Could not create graph communicator");
  Communicator::reset(&result);
}

std::tuple<t_uint, t_uint, bool> DistGraphCommunicator::nedges() const {
  if(not is_valid())
    return 0;

  int indeg, outdeg, weighted;
  auto const error = MPI_Dist_graph_neighbors_count(**this, &indeg, &outdeg, &weighted);
  if(error != MPI_SUCCESS)
    throw std::runtime_error("Could not retrieve the number of input edges");
  return {static_cast<t_uint>(indeg), static_cast<t_uint>(outdeg), weighted != 0};
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

} /* optimet::mpi */
} /* optimet  */
