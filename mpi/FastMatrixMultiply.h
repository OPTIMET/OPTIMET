#ifndef OPTIMET_MPI_FAST_MATRIX_MULTIPLY_H

#include "../FastMatrixMultiply.h"
#include "mpi/GraphCommunicator.h"
#include <tuple>
#include <utility>
#include <vector>

namespace optimet {
namespace mpi {
namespace details {
//! \brief Splits interactions into local and non-local processes
//! \details If true, then that element will be computed locally.
Matrix<bool> local_interactions(t_int nscatterers);
//! \brief Distribution of input/output vector
Vector<t_int> vector_distribution(t_int nscatterers, t_int nprocs);
//! \brief Splits interactions into local and non-local processes
inline Matrix<t_int> matrix_distribution(Matrix<bool> const &local, Vector<t_int> const &columns) {
  return local.select(columns.transpose().replicate(columns.size(), 1),
                      columns.replicate(1, columns.size()));
}
//! Combines splitting alongst local-non-local and alongs columns
inline Matrix<t_int> matrix_distribution(t_int nscatterers, t_int nprocs) {
  return matrix_distribution(local_interactions(nscatterers),
                             vector_distribution(nscatterers, nprocs));
}

//! Figures out graph connectivity for given distribution
std::vector<std::set<t_uint>>
local_graph_edges(Matrix<bool> const &locals, Vector<t_int> const &vector_distribution);
//! Figures out graph connectivity for given distribution
std::vector<std::set<t_uint>>
non_local_graph_edges(Matrix<bool> const &nonlocals, Vector<t_int> const &vector_distribution);

//! Amount of input data to receive from neighboring processes
std::vector<int> neighborhood_input_counts(Matrix<bool> const &nonlocals,
                                           Vector<t_int> const &vector_distribution,
                                           std::vector<Scatterer> const &scatterers, t_uint rank);
}

//! \brief MPI version of the Fast-Matrix-Multiply
//! \details This operation distributes a Matrix-vector multiplication. The only data of import
//! outside this operation are the input/output vectors. The matrix/operation can be distributed as
//! we best see fit.
//!
//! The vectors are distributed across procs homogeneously: each proc owns all the coefficients
//! associated with a set of particles, contiguous in the array of input particles.
//!
//! To overlap calculations and communications, we perform several steps:
//!
//! 1. distribute to the processes that require it the input coefficients owned by this process.
//! 2. compute parts of the vector-matrix calculation that require only knowledge of the owned
//! coefficients. The output of this calculation affects the coefficients owned by other procs.
//! 3. distribute the computations from the previous step to the relevant procs
//! 4. receive the input data from step 1
//! 5. perform parts of the matrix-vector multiplication using the data received in 4 that affects
//! output coefficients owned by this process
//! 6. receive data from 3.
//! 7. performs a reductions over data received in 6 and computed in 5
//!
//! The crux is to separate those calculations that require data from other processes but outputs
//! only this to process, from calculations requiring only data from this process but outputs to any
//! process. There diagonal part of the matrix-vector multiplication can be computed eitehr at step
//! 2 or 5.
class FastMatrixMultiply {

public:
  //! Helper class to perform steps 1 and 4
  class DistributeInput;
  //! Helper class to perform steps 3, 6, and 7
  class ReduceComputation;

  //! Creates an MPI fast-matrix-multiply
  //! \param[in] em_background: Electromagnetic properties of the background medium
  //! \param[in] wavenumber: of the impinging wave
  //! \param[in] scatterers: array of scatterers
  //! \param[in] local_nonlocal: Each element refers to the coupling of a pair of particles. Those
  //!                            elements that are true will be computed in step 2 of the algorithm,
  //!                            and those that are false in step 5. This process will compute only
  //!                            the relevant elements that are along it's own columns (for step 2)
  //!                            and rows (for step 5).
  //! \param[in] vector_distribution: defines the rank of that own each element in the input
  //!                                 vector. It should define *contiguous* ranges.
  //! \param[in] comm: Communicator from which to create the graph communicators for the two
  //!                  communication steps. This communicator should hold all and only those
  //!                  processes involved in the matrix-vector multiplication.
  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers, Matrix<bool> const &local_nonlocal,
                     Vector<t_int> const &vector_distribution,
                     mpi::Communicator const &comm = mpi::Communicator())
      : FastMatrixMultiply(em_background, wavenumber, scatterers, local_nonlocal,
                           local_dist(local_nonlocal, vector_distribution, comm.rank()),
                           nonlocal_dist(local_nonlocal, vector_distribution, comm.rank()),
                           vector_distribution, comm) {}
  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers,
                     mpi::Communicator const &comm = mpi::Communicator())
      : FastMatrixMultiply(em_background, wavenumber, scatterers,
                           details::local_interactions(scatterers.size()),
                           details::vector_distribution(scatterers.size(), comm.size()), comm) {}
  FastMatrixMultiply(t_real wavenumber, std::vector<Scatterer> const &scatterers,
                     mpi::Communicator const &comm = mpi::Communicator())
      : FastMatrixMultiply(ElectroMagnetic(), wavenumber, scatterers, comm) {}

  //! \brief Applies fast matrix multiplication to effective incident field
  void operator()(Vector<t_complex> const &in, Vector<t_complex> &out) const;
  //! \brief Applies fast matrix multiplication to effective incident field
  Vector<t_complex> operator()(Vector<t_complex> const &in) const;
  //! \brief Applies fast matrix multiplication to effective incident field
  Vector<t_complex> operator*(Vector<t_complex> const &in) const { return operator()(in); }

  //! Local rows
  t_uint rows() const { return nonlocal_fmm_.rows(); }
  //! Local cols
  t_uint cols() const { return local_fmm_.cols(); }

private:
  //! Computed with local input vector
  optimet::FastMatrixMultiply local_fmm_;
  //! Computed with non-local input vector
  optimet::FastMatrixMultiply nonlocal_fmm_;
  //! MPI communicator over which to distribute input data
  mpi::GraphCommunicator distribute_comm_;
  //! MPI communicator over which to reduce output data
  mpi::GraphCommunicator reduce_comm_;

  template <class T0, class T1>
  Matrix<bool> local_dist(Eigen::MatrixBase<T0> const &locals, Eigen::MatrixBase<T1> const &vecdist,
                          t_uint rank) {
    return locals.array() && (vecdist.transpose().array() == rank).replicate(vecdist.size(), 1);
  }
  template <class T0, class T1>
  Matrix<bool> nonlocal_dist(Eigen::MatrixBase<T0> const &locals,
                             Eigen::MatrixBase<T1> const &vecdist, t_uint rank) {
    return (locals.array() == false) && (vecdist.array() == rank).replicate(1, vecdist.size());
  }

  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers, Matrix<bool> const &lnl,
                     Matrix<bool> const &local_dist, Matrix<bool> const &nonlocal_dist,
                     Vector<t_int> const &vector_distribution,
                     mpi::Communicator const &comm = mpi::Communicator());
};

class FastMatrixMultiply::DistributeInput {
public:
  DistributeInput(GraphCommunicator const &comm, Matrix<bool> const &allowed,
                  Vector<int> const &distribution, std::vector<Scatterer> const &scatterers);
  DistributeInput(GraphCommunicator const &comm, Matrix<bool> const &allowed,
                  Vector<int> const &distribution, Vector<int> const &sizes);

  //! Performs input distribution request
  template <class T0, class T1>
  Request
  send(Eigen::PlainObjectBase<T0> const &input, Eigen::PlainObjectBase<T1> const &receiving) const {
    return comm.iallgather(input, receiving, receive_counts);
  }

  template <class T0, class T1>
  void synthesize(Eigen::MatrixBase<T0> const &received,
                  Eigen::PlainObjectBase<T1> const &synthesis) const;

protected:
  GraphCommunicator comm;
  std::vector<std::tuple<t_uint, t_uint, t_uint>> relocate_receive;
  std::vector<int> receive_counts;
  t_uint synthesis_size;
};

class FastMatrixMultiply::ReduceComputation {
public:
  ReduceComputation(GraphCommunicator const &comm, Matrix<bool> const &allowed,
                    Vector<int> const &distribution, std::vector<Scatterer> const &scatterers);
  ReduceComputation(GraphCommunicator const &comm, Matrix<bool> const &allowed,
                    Vector<int> const &distribution, Vector<int> const &sizes);

  //! Performs reduction request over computed data
  template <class T0, class T1>
  Request
  send(Eigen::MatrixBase<T0> const &input, Eigen::PlainObjectBase<T1> const &receiving) const;
  //! Performs reduction over received data
  template <class T0, class T1>
  void reduce(Eigen::MatrixBase<T0> const &inout, Eigen::MatrixBase<T1> const &received) const;

  template <class T0, class T1>
  void reduce(Request &&request, Eigen::MatrixBase<T0> const &inout,
              Eigen::MatrixBase<T1> const &received) const;

protected:
  GraphCommunicator comm;
  std::vector<std::tuple<t_uint, t_uint, t_uint>> relocate_receive;
  std::vector<std::tuple<t_uint, t_uint, t_uint>> relocate_send;
  std::vector<int> receive_counts, send_counts;
  std::unique_ptr<Vector<int>> message;
};

template <class T0, class T1>
void FastMatrixMultiply::DistributeInput::synthesize(
    Eigen::MatrixBase<T0> const &received, Eigen::PlainObjectBase<T1> const &synthesis) const {
  const_cast<Eigen::PlainObjectBase<T1> &>(synthesis).resize(synthesis_size, 1);
  const_cast<Eigen::PlainObjectBase<T1> &>(synthesis).fill(0);
  for(auto const &location : relocate_receive) {
    auto const &size = std::get<0>(location);
    auto const &rcv_index = std::get<1>(location);
    auto const &in_index = std::get<2>(location);
    assert(in_index + size <= synthesis.size());
    assert(rcv_index + size <= received.size());
    const_cast<Eigen::PlainObjectBase<T1> &>(synthesis).segment(in_index, size) =
        received.segment(rcv_index, size);
  }
}

template <class T0, class T1>
Request
FastMatrixMultiply::ReduceComputation::send(Eigen::MatrixBase<T0> const &input,
                                            Eigen::PlainObjectBase<T1> const &receiving) const {
  for(auto const &location : relocate_send) {
    auto const &size = std::get<0>(location);
    auto const &message_index = std::get<1>(location);
    auto const &send_index = std::get<2>(location);
    assert(send_index + size <= input.size());
    assert(message_index + size <= message->size());
    message->segment(message_index, size) = input.segment(send_index, size);
  }

  return comm.ialltoall(*message, receiving, send_counts, receive_counts);
}

template <class T0, class T1>
void FastMatrixMultiply::ReduceComputation::reduce(Eigen::MatrixBase<T0> const &inout,
                                                   Eigen::MatrixBase<T1> const &received) const {
  for(auto const &location : relocate_receive) {
    auto const &size = std::get<0>(location);
    auto const &rcv_index = std::get<1>(location);
    auto const &in_index = std::get<2>(location);
    assert(in_index + size <= inout.size());
    assert(rcv_index + size <= received.size());
    const_cast<Eigen::MatrixBase<T0> &>(inout).segment(in_index, size) +=
        received.segment(rcv_index, size);
  }
}

template <class T0, class T1>
void FastMatrixMultiply::ReduceComputation::reduce(Request &&request,
                                                   Eigen::MatrixBase<T0> const &inout,
                                                   Eigen::MatrixBase<T1> const &received) const {
  assert(request);
  auto const deleter = request.get_deleter();
  deleter(request.release());
  return reduce(inout, received);
}
}
}

#endif
