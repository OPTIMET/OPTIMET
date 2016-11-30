#ifndef OPTIMET_MPI_GRAPH_COMMUNICATOR_H
#define OPTIMET_MPI_GRAPH_COMMUNICATOR_H

#include "Types.h"

#ifdef OPTIMET_MPI

#include "mpi/Collectives.h"
#include "mpi/Communicator.h"
#include "mpi/RegisteredTypes.h"
#include <memory>
#include <numeric>
#include <set>
#include <vector>
#include <iostream>

namespace optimet {
namespace mpi {

//! Type of a request
typedef std::shared_ptr<MPI_Request const> Request;
//! A request that calls wait when going out of scope
Request make_wait_on_delete(MPI_Request *const request);

class DistGraphCommunicator : public Communicator {
public:
  //! Creates a graph communicator
  DistGraphCommunicator(mpi::Communicator const &comm, std::vector<int> const &sources,
                        std::vector<int> const &destinations, bool reorder = false);
  DistGraphCommunicator(std::vector<int> const &sources, std::vector<int> const &destination,
                        bool reorder = false)
      : DistGraphCommunicator(Communicator(), sources, destination, reorder) {}

  //! Number of input/output edges and wether the graph is weighted
  std::tuple<t_uint, t_uint, bool> nedges() const;
  //! Broadcast whole input vector to neighbors
  template <class T0, class T1>
  typename std::enable_if<is_registered_type<typename T0::Scalar>::value and
                              is_registered_type<typename T1::Scalar>::value,
                          Request>::type
  iallgatherv(Eigen::PlainObjectBase<T0> const &input, Eigen::PlainObjectBase<T1> const &out,
              std::vector<int> rcvcounts) const;

  //! Blocking Send/Receive single scalar to neighbor
  template <class T>
  typename std::enable_if<is_registered_type<T>::value, std::vector<T>>::type
  allgather(T const &input) const;
};

// Send/Receive single scalar to neighbor
template <class T>
typename std::enable_if<is_registered_type<T>::value, std::vector<T>>::type
DistGraphCommunicator::allgather(T const &input) const {
  if(not is_valid())
    return std::vector<T>();

  // number of sources
  auto const Nsources = std::get<0>(nedges());
  auto const Ndestinations = std::get<1>(nedges());
  std::vector<T> result(Nsources + 1);
  std::fill(result.begin(), result.end(), 666);
  auto const error = MPI_Neighbor_allgather(
      &input, std::min<int>(Ndestinations, 1), mpi::registered_type(input), result.data(),
      std::min<int>(Nsources, 1), mpi::registered_type(input), **this);
  if(error != MPI_SUCCESS)
    throw std::runtime_error("Gathering scalar data over graph communicator failed");

  auto const back = result.back();
  result.pop_back();
  return result;
}

template <class T0, class T1>
typename std::enable_if<is_registered_type<typename T0::Scalar>::value and
                            is_registered_type<typename T1::Scalar>::value,
                        Request>::type
DistGraphCommunicator::iallgatherv(Eigen::PlainObjectBase<T0> const &input,
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

  auto const Nin = std::get<0>(nedges());

  std::unique_ptr<MPI_Request> request(new MPI_Request);
  typename T1::Scalar dummy(0);
  typename T1::Scalar *const out_buffer =
      out.size() == 0 ? &dummy : const_cast<Eigen::PlainObjectBase<T1> &>(out).data();
  int const idummy(0);
  int const *const rcvcount_ptr = rcvcounts.size() == 0 ? &idummy : rcvcounts.data();
  int const *const outdisp_ptr = rcvcounts.size() == 0 ? &idummy : out_disps.data();
  typename T0::Scalar const input_dummy = 0;
  typename T0::Scalar const *const input_ptr = input.size() == 0 ? &input_dummy : input.data();

  auto const error = MPI_Ineighbor_allgatherv(
      input_ptr, Nin == 0 ? 0 : input.size(), mpi::Type<typename T0::Scalar>::value, out_buffer,
      rcvcount_ptr, outdisp_ptr, mpi::Type<typename T1::Scalar>::value, **this, request.get());

  if(error != MPI_SUCCESS)
    throw std::runtime_error("Gathering data over graph communicator failed");

  return make_wait_on_delete(request.release());
}
}
}
#endif
#endif
