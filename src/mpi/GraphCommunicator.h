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
typedef std::unique_ptr<MPI_Request, void (*)(MPI_Request *const)> Request;
//! A request that calls wait when going out of scope
Request mpi_request_wait_on_delete(MPI_Request *const request);
//! Wait for request to complete
inline void wait(Request &&request) { Request const locally_scoped(std::move(request)); }

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

  //! Rank of neighbors
  std::vector<t_uint> neighborhood(int rank) const;
  std::vector<t_uint> neighborhood() const { return neighborhood(rank()); }
  //! Non-blocking Gather over graph
  template <class T0, class T1>
  typename std::enable_if<is_registered_type<typename T0::Scalar>::value and
                              is_registered_type<typename T1::Scalar>::value,
                          Request>::type
  iallgather(Eigen::PlainObjectBase<T0> const &input, Eigen::PlainObjectBase<T1> const &out,
             std::vector<int> const &rcvcounts) const;

  //! Non-blocking all-to-all over graph
  template <class T0, class T1>
  typename std::enable_if<is_registered_type<typename T0::Scalar>::value and
                              is_registered_type<typename T1::Scalar>::value,
                          Request>::type
  ialltoall(Eigen::PlainObjectBase<T0> const &input, Eigen::PlainObjectBase<T1> const &out,
            std::vector<int> const &send_counts, std::vector<int> const &receive_counts) const;

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
                              std::vector<int> const &rcvcounts) const {
  if(not is_valid()) {
    const_cast<Eigen::PlainObjectBase<T1> &>(out).resize(0);
    return mpi_request_wait_on_delete(nullptr);
  }

  std::vector<int> receive_disps(rcvcounts.size());
  if(receive_disps.size() != 0) {
    receive_disps[0] = 0;
    for(std::vector<int>::size_type i(1); i < rcvcounts.size(); ++i)
      receive_disps[i] = rcvcounts[i - 1] + receive_disps[i - 1];

    const_cast<Eigen::PlainObjectBase<T1> &>(out).resize(receive_disps.back() + rcvcounts.back());
  } else
    const_cast<Eigen::PlainObjectBase<T1> &>(out).resize(0);

  std::unique_ptr<MPI_Request> request(new MPI_Request);
  typename T1::Scalar dummy(0);
  typename T1::Scalar *const out_buffer =
      out.size() == 0 ? &dummy : const_cast<Eigen::PlainObjectBase<T1> &>(out).data();
  int const idummy(0);
  int const *const rcvcount_ptr = rcvcounts.size() == 0 ? &idummy : rcvcounts.data();
  int const *const outdisp_ptr = rcvcounts.size() == 0 ? &idummy : receive_disps.data();
  typename T0::Scalar const input_dummy = 0;
  typename T0::Scalar const *const input_ptr = input.size() == 0 ? &input_dummy : input.data();

  auto const N = neighborhood_size();
  auto const error = MPI_Ineighbor_allgatherv(
      input_ptr, N == 0 ? 0 : input.size(), mpi::Type<typename T0::Scalar>::value, out_buffer,
      rcvcount_ptr, outdisp_ptr, mpi::Type<typename T1::Scalar>::value, **this, request.get());

  if(error != MPI_SUCCESS)
    throw std::runtime_error("Gathering data over graph communicator failed");

  return mpi_request_wait_on_delete(request.release());
}

template <class T0, class T1>
typename std::enable_if<is_registered_type<typename T0::Scalar>::value and
                            is_registered_type<typename T1::Scalar>::value,
                        Request>::type
GraphCommunicator::ialltoall(Eigen::PlainObjectBase<T0> const &input,
                             Eigen::PlainObjectBase<T1> const &out,
                             std::vector<int> const &send_counts,
                             std::vector<int> const &receive_counts) const {
  if(not is_valid()) {
    const_cast<Eigen::PlainObjectBase<T1> &>(out).resize(0);
    return mpi_request_wait_on_delete(nullptr);
  }

  std::vector<int> receive_disps(receive_counts.size());
  if(receive_disps.size() != 0) {
    receive_disps[0] = 0;
    for(std::vector<int>::size_type i(1); i < receive_counts.size(); ++i)
      receive_disps[i] = receive_counts[i - 1] + receive_disps[i - 1];

    const_cast<Eigen::PlainObjectBase<T1> &>(out).resize(receive_disps.back() +
                                                         receive_counts.back());
  } else
    const_cast<Eigen::PlainObjectBase<T1> &>(out).resize(0);

  std::vector<int> send_disps(send_counts.size());
  if(send_disps.size() != 0) {
    send_disps[0] = 0;
    for(std::vector<int>::size_type i(1); i < send_counts.size(); ++i)
      send_disps[i] = send_counts[i - 1] + send_disps[i - 1];

    if(input.size() != send_disps.back() + send_counts.back())
      throw std::out_of_range("Input vector does not match send_counts");
  }

  std::unique_ptr<MPI_Request> request(new MPI_Request);
  typename T1::Scalar receive_dummy(0);
  typename T1::Scalar *const receive_buffer =
      out.size() == 0 ? &receive_dummy : const_cast<Eigen::PlainObjectBase<T1> &>(out).data();
  int const idummy(0);
  int const *const receive_count_ptr = receive_counts.size() == 0 ? &idummy : receive_counts.data();
  int const *const receive_disp_ptr = receive_counts.size() == 0 ? &idummy : receive_disps.data();

  typename T0::Scalar const send_dummy = 0;
  typename T0::Scalar const *const send_ptr = input.size() == 0 ? &send_dummy : input.data();
  int const *const send_count_ptr = send_counts.size() == 0 ? &idummy : send_counts.data();
  int const *const send_disp_ptr = send_counts.size() == 0 ? &idummy : send_disps.data();

  auto const error = MPI_Ineighbor_alltoallv(
      send_ptr, send_count_ptr, send_disp_ptr, mpi::registered_type(send_dummy), receive_buffer,
      receive_count_ptr, receive_disp_ptr, mpi::registered_type(receive_dummy), **this,
      request.get());

  if(error != MPI_SUCCESS)
    throw std::runtime_error("All-to-all data over graph communicator failed");

  return mpi_request_wait_on_delete(request.release());
}
}
}
#endif
#endif
