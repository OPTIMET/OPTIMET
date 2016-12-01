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

std::vector<int> neighborhood_input_counts(Matrix<bool> const &nonlocals,
                                           Vector<t_int> const &vector_distribution,
                                           std::vector<Scatterer> const &scatterers, t_uint rank) {
  // graph should be symmetric
  assert(nonlocals == nonlocals.transpose());
  std::map<int, int> counts;
  for(Matrix<bool>::Index i(0); i < nonlocals.cols(); ++i) {
    auto const other_proc = vector_distribution[i];
    if(other_proc == rank)
      continue;
    bool found = false;
    for(Matrix<bool>::Index j(0); j < nonlocals.rows(); ++j)
      if(vector_distribution[j] == rank and nonlocals(i, j) == true) {
        found = true;
        break;
      }
    if(found) {
      if(counts.find(other_proc) == counts.end())
        counts[other_proc] = 0;
      counts[other_proc] += 2 * scatterers[i].nMax * (scatterers[i].nMax + 2);
    }
  }
  std::vector<int> result;
  std::transform(counts.begin(), counts.end(), std::back_inserter(result),
                 [](std::map<int, int>::const_reference a) { return a.second; });
  return result;
}
}

// optimet::FastMatrixMultiply
// FastMatrixMultiply::create_column(ElectroMagnetic const &em_background, t_real wavenumber,
//                                   std::vector<Scatterer> const &scatterers, t_uint nprocs,
//                                   t_uint rank) {
//   auto const range = details::column_range(scatterers.size(), nprocs, rank);
//   auto const n = range.second - range.first;
//   if(n == 0)
//     return optimet::FastMatrixMultiply(em_background, wavenumber, std::vector<Scatterer>());
//   std::vector<Scatterer> objects;
//   objects.reserve(n);
//   std::copy(scatterers.begin() + range.first, scatterers.begin() + range.second,
//             std::back_inserter(objects));
//   auto const out = rank - range.first;
//   return optimet::FastMatrixMultiply(em_background, wavenumber, objects, {out, out + 1},
//                                      {0, objects.size()});
// }
//
// optimet::FastMatrixMultiply FastMatrixMultiply::create_row(ElectroMagnetic const &em_background,
//                                                            t_real wavenumber,
//                                                            std::vector<Scatterer> const
//                                                            &scatterers,
//                                                            t_uint nprocs, t_uint rank) {
//   auto const range = details::column_range(scatterers.size(), nprocs, rank);
//   auto const n = range.second - range.first;
//   if(n == 0 or n == scatterers.size())
//     return optimet::FastMatrixMultiply(em_background, wavenumber, std::vector<Scatterer>());
//   std::vector<Scatterer> objects;
//   objects.reserve(objects.size() - n + 1);
//   objects.push_back(scatterers[rank]);
//   std::copy(scatterers.begin(), scatterers.begin() + range.first, std::back_inserter(objects));
//   std::copy(scatterers.begin() + range.second, scatterers.end(), std::back_inserter(objects));
//   return optimet::FastMatrixMultiply(em_background, wavenumber, objects, {1, objects.size()},
//                                      {0, 1});
// }
}
} // optimet namespace
