#ifndef OPTIMET_MPI_FAST_MATRIX_MULTIPLY_H
#define OPTIMET_MPI_FAST_MATRIX_MULTIPLY_H

#include "../FastMatrixMultiply.h"
#include "mpi/GraphCommunicator.h"
#include <array>
#include <utility>
#include <vector>

namespace optimet {
namespace mpi {
namespace details {
//! \brief Splits interactions into local and non-local processes
//! \details True elements are computed using local input data. This function creates a banded
//! matrix, with a given number of diagonal and subdiagonal set to true. The return is a symmetric
//! matrix.
Matrix<bool> local_interactions(t_int nscatterers, t_int diagonal);
//! \brief Splits interactions into local and non-local processes
//! \details Creates a banded matrix with roughly as many local as non-local interactions.
inline Matrix<bool> local_interactions(t_int nscatterers) {
  return local_interactions(nscatterers, std::max<t_int>(1, (nscatterers - 1) / 2));
};
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
graph_edges(Matrix<bool> const &locals, Vector<t_int> const &vector_distribution);
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
  //! \param[in] diagonal: In some constructors, the `locals` matrix is constructed as a diagonal
  //!                      banded matrix with this number of subdiagonals set to local (computations
  //!                      from locally available input data).
  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers, Matrix<bool> const &locals,
                     Vector<t_int> const &vector_distribution,
                     Communicator const &comm = Communicator())
      : FastMatrixMultiply(
            em_background, wavenumber, scatterers, locals,
            // reordering in graph communicators would require re-mapping vector_distribution
            GraphCommunicator(
                comm, details::graph_edges(locals.array() == false, vector_distribution), false),
            // reordering in graph communicators would require re-mapping vector_distribution
            GraphCommunicator(comm, details::graph_edges(locals, vector_distribution), false),
            vector_distribution, comm) {}
  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers, t_int diagonal,
                     Vector<t_int> const &vector_distribution,
                     Communicator const &comm = Communicator())
      : FastMatrixMultiply(em_background, wavenumber, scatterers,
                           details::local_interactions(scatterers.size(), diagonal),
                           vector_distribution, comm) {}
  FastMatrixMultiply(t_real wavenumber, std::vector<Scatterer> const &scatterers, t_int diagonal,
                     Vector<t_int> const &vector_distribution,
                     Communicator const &comm = Communicator())
      : FastMatrixMultiply(ElectroMagnetic(), wavenumber, scatterers, diagonal, vector_distribution,
                           comm) {}
  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers,
                     Communicator const &comm = Communicator())
      : FastMatrixMultiply(em_background, wavenumber, scatterers,
                           details::local_interactions(scatterers.size()),
                           details::vector_distribution(scatterers.size(), comm.size()), comm) {}
  FastMatrixMultiply(t_real wavenumber, std::vector<Scatterer> const &scatterers,
                     Communicator const &comm = Communicator())
      : FastMatrixMultiply(ElectroMagnetic(), wavenumber, scatterers, comm) {}
  FastMatrixMultiply(t_real wavenumber, std::vector<Scatterer> const &scatterers, t_int diagonal,
                     Communicator const &comm = Communicator())
      : FastMatrixMultiply(wavenumber, scatterers, diagonal,
                           details::vector_distribution(scatterers.size(), comm.size()), comm) {}
  FastMatrixMultiply(t_real wavenumber, std::vector<Scatterer> const &scatterers,
                     Matrix<bool> const &local_nonlocal, Vector<t_int> const &vector_distribution,
                     Communicator const &comm = Communicator())
      : FastMatrixMultiply(ElectroMagnetic(), wavenumber, scatterers, local_nonlocal,
                           vector_distribution, comm) {}
  FastMatrixMultiply(t_real wavenumber, std::vector<Scatterer> const &scatterers,
                     Matrix<bool> const &local_nonlocal, Communicator const &comm = Communicator())
      : FastMatrixMultiply(wavenumber, scatterers, local_nonlocal,
                           details::vector_distribution(scatterers.size(), comm.size()), comm) {}

  //! \brief Applies fast matrix multiplication to effective incident field
  void operator()(Vector<t_complex> const &in, Vector<t_complex> &out) const;
  //! \brief Applies fast matrix multiplication to effective incident field
  Vector<t_complex> operator()(Vector<t_complex> const &in) const;
  //! \brief Applies fast matrix multiplication to effective incident field
  Vector<t_complex> operator*(Vector<t_complex> const &in) const { return operator()(in); }
  //! \brief Applies transpose fast matrix multiplication to effective incident field
  void transpose(Vector<t_complex> const &in, Vector<t_complex> &out) const;
  //! \brief Applies transpose fast matrix multiplication to effective incident field
  Vector<t_complex> transpose(Vector<t_complex> const &in) const;
  //! \brief Applies conjugate fast matrix multiplication to effective incident field
  void conjugate(Vector<t_complex> const &in, Vector<t_complex> &out) const {
    operator()(in.conjugate(), out);
    out = out.conjugate();
  }
  //! \brief Applies conjugate fast matrix multiplication to effective incident field
  Vector<t_complex> conjugate(Vector<t_complex> const &in) const {
    return operator()(in.conjugate()).conjugate();
  }
  //! \brief Applies adjoint fast matrix multiplication to effective incident field
  void adjoint(Vector<t_complex> const &in, Vector<t_complex> &out) const {
    transpose(in.conjugate(), out);
    out = out.conjugate();
  }
  //! \brief Applies conjugate fast matrix multiplication to effective incident field
  Vector<t_complex> adjoint(Vector<t_complex> const &in) const {
    return transpose(in.conjugate()).conjugate();
  }

  //! Local rows
  t_uint rows() const { return nonlocal_fmm_.rows(); }
  //! Local cols
  t_uint cols() const { return local_fmm_.cols(); }

  //! Reconstructs output according to argument indices
  template <class T0, class T1>
  static void reconstruct(std::vector<std::array<t_uint, 3>> const &indices,
                          Eigen::MatrixBase<T0> const &input,
                          Eigen::PlainObjectBase<T1> const &output, bool sum = false);
  //! Reconstructs output according to argument indices
  template <class T0>
  static Vector<typename T0::Scalar>
      reconstruct(std::vector<std::array<t_uint, 3>> const &indices,
                  Eigen::MatrixBase<T0> const &input, bool sum = false);

  //! Creates indices needed to reconstruct ouputs
  template <class T0, class T1, class FUNCTOR>
  static std::vector<std::array<t_uint, 3>>
  reconstruct(std::vector<t_uint> const &neighborhood, Eigen::DenseBase<T0> const &allowed,
              Eigen::DenseBase<T1> const &sizes, FUNCTOR const &func);

  class DistributeInput {
  public:
    DistributeInput(GraphCommunicator const &comm, Matrix<bool> const &allowed,
                    Vector<int> const &distribution, std::vector<Scatterer> const &scatterers);
    DistributeInput(GraphCommunicator const &comm, Matrix<bool> const &allowed,
                    Vector<int> const &distribution, Vector<int> const &sizes);

    //! Performs input distribution request
    template <class T0, class T1>
    Request send(Eigen::PlainObjectBase<T0> const &input,
                 Eigen::PlainObjectBase<T1> const &receiving) const {
      return comm.iallgather(input, receiving, receive_counts);
    }

    //! Creates an input vector for the matrix-multiplication terms from un-owned input
    template <class T0, class T1>
    void synthesize(Eigen::MatrixBase<T0> const &received,
                    Eigen::PlainObjectBase<T1> const &synthesis) const {
      FastMatrixMultiply::reconstruct(relocate_receive, received, synthesis);
    }
    template <class T0>
    Vector<typename T0::Scalar> synthesize(Eigen::MatrixBase<T0> const &received) const {
      Vector<typename T0::Scalar> result;
      synthesize(received, result);
      return result;
    }

  protected:
    GraphCommunicator comm;
    std::vector<std::array<t_uint, 3>> relocate_receive;
    std::vector<int> receive_counts;
    t_uint synthesis_size;
  };

  class ReduceComputation {
  public:
    ReduceComputation(GraphCommunicator const &comm, Matrix<bool> const &allowed,
                      Vector<int> const &distribution, std::vector<Scatterer> const &scatterers);
    ReduceComputation(GraphCommunicator const &comm, Matrix<bool> const &allowed,
                      Vector<int> const &distribution, Vector<int> const &sizes);

    //! Performs reduction request over computed data
    template <class T0, class T1, class T2>
    Request send(Eigen::MatrixBase<T0> const &input, Eigen::PlainObjectBase<T1> const &send_buffer,
                 Eigen::PlainObjectBase<T2> const &receiving) const;
    //! \brief Performs reduction over received data
    //! \details The operation is equivalent to reconstructing the output vectors such and doing a
    //! reduction. In practice, this operation performs the sum over the different
    //! matrix-multiplication terms computed by the different processes.
    template <class T0, class T1>
    void
    reduce(Eigen::PlainObjectBase<T0> const &inout, Eigen::MatrixBase<T1> const &received) const {
      FastMatrixMultiply::reconstruct(relocate_receive, received, inout, true);
    }

  protected:
    GraphCommunicator comm;
    std::vector<std::array<t_uint, 3>> relocate_receive;
    std::vector<std::array<t_uint, 3>> relocate_send;
    std::vector<int> receive_counts, send_counts;
    t_uint message_size;
  };

private:
  //! Computed with local input vector
  optimet::FastMatrixMultiply local_fmm_;
  //! Computed with non-local input vector
  optimet::FastMatrixMultiply nonlocal_fmm_;
  //! Computed with local input vector
  optimet::FastMatrixMultiply transpose_local_fmm_;
  //! Computed with non-local input vector
  optimet::FastMatrixMultiply transpose_nonlocal_fmm_;
  //! Distribute input for nonlocal_fmm_
  DistributeInput distribute_input_;
  //! Communicate results from local_fmm_ and reduce all results
  ReduceComputation reduce_computation_;
  //! Reconstruction indices for output of nonlocal_fmm_
  std::vector<std::array<t_uint, 3>> nonlocal_indices_;
  //! Reconstruction indices for input to local_fmm_
  std::vector<std::array<t_uint, 3>> local_indices_;

  //! End-point of the constructor chain
  FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                     std::vector<Scatterer> const &scatterers, Matrix<bool> const &locals,
                     GraphCommunicator const &distribute_comm, GraphCommunicator const &reduce_comm,
                     Vector<t_int> const &vector_distribution,
                     Communicator const &comm = Communicator());
};

template <class T0, class T1>
void FastMatrixMultiply::reconstruct(std::vector<std::array<t_uint, 3>> const &indices,
                                     Eigen::MatrixBase<T0> const &input,
                                     Eigen::PlainObjectBase<T1> const &output, bool sum) {
  auto const out_size = [](t_uint prior, std::array<t_uint, 3> const &item) {
    return std::max<t_uint>(prior, item[0] + item[2]);
  };
  auto const N = std::accumulate(indices.begin(), indices.end(), t_uint(0), out_size);
#ifndef NDEBUG
  auto const in_size = [](t_uint prior, std::array<t_uint, 3> const &item) {
    return std::max<t_uint>(prior, item[0] + item[1]);
  };
  assert(std::accumulate(indices.begin(), indices.end(), t_uint(0), in_size) <= input.size());
#endif
  if(output.size() < N) {
    auto const original = output.derived();
    const_cast<Eigen::PlainObjectBase<T1> &>(output).resize(N);
    const_cast<Eigen::PlainObjectBase<T1> &>(output).head(original.size()) = original;
    const_cast<Eigen::PlainObjectBase<T1> &>(output).tail(output.size() - original.size()).fill(0);
  }
  if(sum)
    for(auto const &location : indices)
      const_cast<Eigen::PlainObjectBase<T1> &>(output).segment(location[2], location[0]) +=
          input.segment(location[1], location[0]).template cast<typename T1::Scalar>();
  else
    for(auto const &location : indices)
      const_cast<Eigen::PlainObjectBase<T1> &>(output).segment(location[2], location[0]) =
          input.segment(location[1], location[0]).template cast<typename T1::Scalar>();
}

template <class T0>
Vector<typename T0::Scalar>
    FastMatrixMultiply::reconstruct(std::vector<std::array<t_uint, 3>> const &indices,
                                    Eigen::MatrixBase<T0> const &input, bool sum) {
  Vector<typename T0::Scalar> result;
  reconstruct(indices, input, result, sum);
  return result;
}

template <class T0, class T1, class FUNCTOR>
std::vector<std::array<t_uint, 3>>
FastMatrixMultiply::reconstruct(std::vector<t_uint> const &neighborhood,
                                Eigen::DenseBase<T0> const &in_allowed,
                                Eigen::DenseBase<T1> const &sizes, FUNCTOR const &func) {
  std::vector<std::array<t_uint, 3>> result;
  for(std::vector<t_uint>::size_type i(0), rloc(0); i < neighborhood.size(); ++i) {
    auto const out_allowed = func(neighborhood[i]);
    assert(out_allowed.size() == in_allowed.size());
    for(t_uint j(0), iloc(0); j < out_allowed.size(); ++j) {
      if(not in_allowed(j))
        continue;
      if(out_allowed(j) == true) {
        result.push_back(std::array<t_uint, 3>{{static_cast<t_uint>(sizes[j]), rloc, iloc}});
        rloc += sizes[j];
      }
      iloc += sizes[j];
    }
  }
  return result;
}

template <class T0, class T1, class T2>
Request
FastMatrixMultiply::ReduceComputation::send(Eigen::MatrixBase<T0> const &input,
                                            Eigen::PlainObjectBase<T1> const &send_buffer,
                                            Eigen::PlainObjectBase<T2> const &receiving) const {
  FastMatrixMultiply::reconstruct(relocate_send, input, send_buffer);
  return comm.ialltoall(send_buffer, receiving, send_counts, receive_counts);
}
}
}

#endif
