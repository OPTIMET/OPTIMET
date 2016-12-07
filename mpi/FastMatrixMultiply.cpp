#include "Types.h"
#include "mpi/FastMatrixMultiply.h"
#include <iostream>

namespace optimet {
namespace mpi {
namespace details {
Matrix<bool> local_interactions(t_int nscatterers) {
  assert(nscatterers >= 0);
  // banded is true where local computations will happen
  Eigen::Array<bool, Eigen::Dynamic, Eigen::Dynamic> band =
      Matrix<bool>::Zero(nscatterers, nscatterers);
  auto const half = std::min((nscatterers + 1) / 2, nscatterers - 1);
  for(t_int i(0); i <= half; ++i) {
    t_int const n = std::max((2 * nscatterers + 1) / 4, 2);
    auto const start = std::min(std::max(nscatterers - n, 0), std::max(i - n / 2, 0));
    auto const size = std::min(start + n, nscatterers) - start;
    band.row(i).segment(start, size).fill(true);
  }
  return (band || band.transpose() || band.reverse() || band.reverse().transpose()).matrix();
}

Vector<t_int> vector_distribution(t_int nscatterers, t_int nprocs) {
  assert(nprocs > 0);
  // horizontal consists of horizontal bands colored with the relevant process
  auto const N = std::min(nscatterers, nprocs);
  Vector<t_int> result = Vector<t_int>::Constant(nscatterers, N);
  for(t_int i(0); i < N; ++i) {
    auto const divided = nscatterers / N;
    auto const remainder = nscatterers % N;
    auto const loc = divided * i + std::min(remainder, i);
    auto const size = divided + (remainder > i ? 1 : 0);
    result.segment(loc, size).fill(i);
  }
  assert(N == 0 or (result.array() < N).all());
  return result;
}

template <class T>
std::vector<std::set<t_uint>>
graph_edges(Eigen::MatrixBase<T> const &considered, Vector<t_int> const &vecdist, bool in_to_out) {
  assert(considered.rows() == considered.cols());
  assert(considered.rows() == vecdist.size());
  auto const nprocs = vecdist.maxCoeff() + 1;
  std::vector<std::set<t_uint>> results(nprocs);
  for(t_uint i(0); i < considered.rows(); ++i)
    for(t_uint j(0); j < considered.cols(); ++j) {
      if(vecdist(i) == vecdist(j) or considered(i, j) == false)
        continue;
      if(in_to_out)
        results[vecdist(j)].insert(vecdist(i));
      else
        results[vecdist(i)].insert(vecdist(j));
    }
  return results;
}

std::vector<std::set<t_uint>>
local_graph_edges(Matrix<bool> const &locals, Vector<t_int> const &vecdist) {
  return graph_edges(locals, vecdist, true);
}
//! Figures out graph connectivity for given distribution
std::vector<std::set<t_uint>>
non_local_graph_edges(Matrix<bool> const &nonlocals, Vector<t_int> const &vecdist) {
  return graph_edges(nonlocals, vecdist, false);
}
}

Vector<t_complex> FastMatrixMultiply::operator()(Vector<t_complex> const &in) const {
  Vector<t_complex> result(rows());
  operator()(in, result);
  return result;
}

void FastMatrixMultiply::operator()(Vector<t_complex> const &input, Vector<t_complex> &out) const {
  // first communicate input data to other processes
  Vector<t_complex> nl_input;
  // auto distribute_request = distribute_comm_.iallgather(input, nl_input, input_counts_);
}

FastMatrixMultiply::FastMatrixMultiply(ElectroMagnetic const &em_background, t_real wavenumber,
                                       std::vector<Scatterer> const &scatterers,
                                       Matrix<bool> const &lnl, Matrix<bool> const &local_dist,
                                       Matrix<bool> const &nonlocal_dist,
                                       Vector<t_int> const &vector_distribution,
                                       mpi::Communicator const &comm)
    : local_fmm_(em_background, wavenumber, scatterers, local_dist),
      nonlocal_fmm_(em_background, wavenumber, scatterers, nonlocal_dist),
      distribute_comm_(comm, GraphCommunicator::symmetrize(
                                 details::local_graph_edges(local_dist, vector_distribution))),
      reduce_comm_(comm, GraphCommunicator::symmetrize(
                             details::non_local_graph_edges(nonlocal_dist, vector_distribution))) {}

namespace {
Vector<int> compute_sizes(std::vector<Scatterer> const &scatterers) {
  Vector<int> result(scatterers.size());
  for(Vector<int>::Index i(0); i < result.size(); ++i)
    result(i) = 2 * scatterers[i].nMax * (scatterers[i].nMax + 2);
  return result;
}
}

FastMatrixMultiply::DistributeInput::DistributeInput(GraphCommunicator const &comm,
                                                     Matrix<bool> const &allowed,
                                                     Vector<int> const &distribution,
                                                     Vector<int> const &sizes)
    : comm(comm) {
  auto const neighborhood = comm.neighborhood();
  // inputs required for local computations
  auto const inputs = [&distribution, &allowed](t_uint rank) {
    return (allowed.array() && (distribution.array() == rank).replicate(1, distribution.size()))
        .colwise()
        .any();
  };
  // inputs required by this process
  auto const local_inputs = inputs(comm.rank());

  auto const owned = (distribution.array() == comm.rank()).eval();
  for(decltype(neighborhood)::size_type i(0), sloc(0); i < neighborhood.size(); ++i) {
    // Helps reconstruct data sent from other procs
    auto const nl_owned = (distribution.array() == neighborhood[i]).eval();
    for(t_uint j(0), iloc(0); j < local_inputs.size(); ++j) {
      if(nl_owned(j) == true and local_inputs(j) == true)
        relocate_receive.emplace_back(sizes[j], sloc, iloc);
      if(nl_owned(j) == true)
        sloc += sizes[j];
      if(local_inputs(j) == true)
        iloc += sizes[j];
    }

    // Amount of data to receive from each proc
    receive_counts.push_back(nl_owned.select(sizes, Vector<int>::Zero(sizes.size())).sum());
  }
  synthesis_size = local_inputs.transpose().select(sizes, Vector<int>::Zero(sizes.size())).sum();
}

FastMatrixMultiply::DistributeInput::DistributeInput(GraphCommunicator const &comm,
                                                     Matrix<bool> const &allowed,
                                                     Vector<int> const &distribution,
                                                     std::vector<Scatterer> const &scatterers)
    : DistributeInput(comm, allowed, distribution, compute_sizes(scatterers)){};

FastMatrixMultiply::ReduceComputation::ReduceComputation(GraphCommunicator const &comm,
                                                         Matrix<bool> const &allowed,
                                                         Vector<int> const &distribution,
                                                         Vector<int> const &sizes)
    : comm(comm) {
  auto const neighborhood = comm.neighborhood();
  auto const comps = [&distribution, &allowed](t_uint rank) {
    return (allowed.array() &&
            (distribution.transpose().array() == rank).replicate(distribution.size(), 1))
        .rowwise()
        .any();
  };
  auto const local_comps = comps(comm.rank());

  auto const owned = (distribution.array() == comm.rank()).eval();
  for(decltype(neighborhood)::size_type i(0), rloc(0), sloc(0); i < neighborhood.size(); ++i) {
    auto const nl_comps = comps(neighborhood[i]);
    // Helps reconstruct data received from other procs
    for(t_uint j(0), iloc(0); j < nl_comps.size(); ++j) {
      if(not owned(j))
        continue;
      if(nl_comps(j) == true) {
        relocate_receive.emplace_back(sizes[j], rloc, iloc);
        rloc += sizes[j];
      }
      iloc += sizes[j];
    }

    // Helps construct data send to other procs
    auto const nl_owned = (distribution.array() == neighborhood[i]).eval();
    for(t_uint j(0), iloc(0); j < local_comps.size(); ++j) {
      if(not local_comps(j))
        continue;
      if(nl_owned(j) == true) {
        relocate_send.emplace_back(sizes[j], sloc, iloc);
        sloc += sizes[j];
      }
      iloc += sizes[j];
    }

    // Amount of data to send to each proc
    send_counts.push_back((local_comps.array() && (distribution.array() == neighborhood[i]))
                              .select(sizes, Vector<int>::Zero(sizes.size()))
                              .sum());
    // Amount of data to receive from each proc
    receive_counts.push_back(
        (nl_comps.array() && owned).select(sizes, Vector<int>::Zero(sizes.size())).sum());
  }
  message.reset(new Vector<int>(std::accumulate(send_counts.begin(), send_counts.end(), 0)));
}

FastMatrixMultiply::ReduceComputation::ReduceComputation(GraphCommunicator const &comm,
                                                         Matrix<bool> const &allowed,
                                                         Vector<int> const &distribution,
                                                         std::vector<Scatterer> const &scatterers)
    : ReduceComputation(comm, allowed, distribution, compute_sizes(scatterers)){};
}
} // optimet namespace
