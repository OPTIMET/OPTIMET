#ifndef OPTIMET_MPI_GRAPH_COMMUNICATOR_H
#define OPTIMET_MPI_GRAPH_COMMUNICATOR_H

#include "Types.h"

#ifdef OPTIMET_MPI

#include "mpi/Collectives.h"
#include "mpi/Communicator.h"
#include "mpi/RegisteredTypes.h"
#include <iostream>
#include <memory>
#include <numeric>
#include <set>
#include <vector>

namespace optimet {
namespace mpi {

//! Type of a request
typedef std::shared_ptr<MPI_Request const> Request;
//! A request that calls wait when going out of scope
Request make_wait_on_delete(MPI_Request *const request);

class GraphCommunicator : public Communicator {
public:
  //! Creates a graph communicator
  GraphCommunicator(mpi::Communicator const &comm, std::vector<std::set<t_uint>> const &graph,
                    bool reorder = false);
  GraphCommunicator(std::vector<std::set<t_uint>> const &graph, bool reorder = false)
      : GraphCommunicator(Communicator(), graph, reorder) {}

  //! Number of neighbors
  t_uint neighborhood_size(int rank) const;
  t_uint neighborhood_size() const { return neighborhood_size(rank()); }
  //! Non-blocking Gather over graph
  template <class T0, class T1>
  typename std::enable_if<is_registered_type<typename T0::Scalar>::value and
                              is_registered_type<typename T1::Scalar>::value,
                          Request>::type
  iallgather(Eigen::PlainObjectBase<T0> const &input, Eigen::PlainObjectBase<T1> const &out,
             std::vector<int> rcvcounts) const;

  //! Blocking Send/Receive single scalar to neighbor
  template <class T>
  typename std::enable_if<is_registered_type<T>::value, std::vector<T>>::type
  allgather(T const &input) const;

  //! Additive symmetrization of graph
  static std::vector<std::set<t_uint>> symmetrize(std::vector<std::set<t_uint>> const &graph);
};

// Send/Receive single scalar to neighbor
template <class T>
typename std::enable_if<is_registered_type<T>::value, std::vector<T>>::type
GraphCommunicator::allgather(T const &input) const {
  if(not is_valid())
    return std::vector<T>();

  // number of sources
  auto const N = neighborhood_size();
  std::vector<T> result(N);
  auto const error = MPI_Neighbor_allgather(
      &input, std::min<int>(N, 1), mpi::registered_type(input), result.data(), std::min<int>(N, 1),
      mpi::registered_type(input), **this);
  if(error != MPI_SUCCESS)
    throw std::runtime_error("Gathering scalar data over graph communicator failed");

  return result;
}

template <class T0, class T1>
typename std::enable_if<is_registered_type<typename T0::Scalar>::value and
                            is_registered_type<typename T1::Scalar>::value,
                        Request>::type
GraphCommunicator::iallgather(Eigen::PlainObjectBase<T0> const &input,
                              Eigen::PlainObjectBase<T1> const &out,
                              std::vector<int> rcvcounts) const {
  if(not is_valid()) {
    const_cast<Eigen::PlainObjectBase<T1> &>(out).resize(0);
    return Request(nullptr);
  }

  std::vector<int> out_disps(rcvcounts.size());
  if(out_disps.size() != 0) {
    out_disps[0] = 0;
    for(std::vector<int>::size_type i(1); i < rcvcounts.size(); ++i)
      out_disps[i] = rcvcounts[i - 1] + out_disps[i - 1];

    const_cast<Eigen::PlainObjectBase<T1> &>(out).resize(out_disps.back() + rcvcounts.back());
  } else
    const_cast<Eigen::PlainObjectBase<T1> &>(out).resize(0);

  std::unique_ptr<MPI_Request> request(new MPI_Request);
  typename T1::Scalar dummy(0);
  typename T1::Scalar *const out_buffer =
      out.size() == 0 ? &dummy : const_cast<Eigen::PlainObjectBase<T1> &>(out).data();
  int const idummy(0);
  int const *const rcvcount_ptr = rcvcounts.size() == 0 ? &idummy : rcvcounts.data();
  int const *const outdisp_ptr = rcvcounts.size() == 0 ? &idummy : out_disps.data();
  typename T0::Scalar const input_dummy = 0;
  typename T0::Scalar const *const input_ptr = input.size() == 0 ? &input_dummy : input.data();

  auto const N = neighborhood_size();
  auto const error = MPI_Ineighbor_allgatherv(
      input_ptr, N == 0 ? 0 : input.size(), mpi::Type<typename T0::Scalar>::value, out_buffer,
      rcvcount_ptr, outdisp_ptr, mpi::Type<typename T1::Scalar>::value, **this, request.get());

  if(error != MPI_SUCCESS)
    throw std::runtime_error("Gathering data over graph communicator failed");

  return make_wait_on_delete(request.release());
}
}
}
#endif
#endif
