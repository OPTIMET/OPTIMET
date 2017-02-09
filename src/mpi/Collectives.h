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

#ifndef OPTIMET_MPI_COLLECTIVES_H
#define OPTIMET_MPI_COLLECTIVES_H

#include "Types.h"
#ifdef OPTIMET_MPI
#include "mpi/Collectives.hpp"
#include "mpi/Communicator.h"
#include <mpi.h>
#include <numeric>
#include <type_traits>
#include <vector>

namespace optimet {
namespace mpi {
inline void barrier(Communicator const &comm) {
  auto const result = MPI_Barrier(*comm);
  if(result != MPI_SUCCESS)
    throw std::runtime_error("Encoutered mpi error");
}

template <class T>
typename std::enable_if<is_registered_type<T>::value, T>::type
broadcast(T const &value, Communicator const &comm, t_uint root) {
  assert(root < comm.size());
  T result = value;
  MPI_Bcast(&result, 1, registered_type(result), root, *comm);
  return result;
}
template <class T>
typename std::enable_if<is_registered_type<T>::value, T>::type
broadcast(Communicator const &comm, t_uint root) {
  assert(root < comm.size());
  assert(root != comm.rank());
  T result;
  MPI_Bcast(&result, 1, registered_type(result), root, *comm);
  return result;
}

template <class T>
Matrix<T> broadcast(Matrix<T> const &mat, Communicator const &comm, t_uint root) {
  assert(root < comm.size());
  auto const nrows = broadcast(mat.rows(), comm, root);
  auto const ncols = broadcast(mat.cols(), comm, root);
  if(comm.rank() == root) {
    MPI_Bcast(const_cast<T *>(mat.data()), nrows * ncols, Type<T>::value, root, *comm);
    return mat;
  } else {
    Matrix<T> result = Matrix<T>::Zero(nrows, ncols);
    MPI_Bcast(result.data(), nrows * ncols, Type<T>::value, root, *comm);
    return result;
  }
}

template <class T>
Vector<T> broadcast(Vector<T> const &mat, Communicator const &comm, t_uint root) {
  assert(root < comm.size());
  auto const size = broadcast(mat.size(), comm, root);
  if(comm.rank() == root) {
    MPI_Bcast(const_cast<T *>(mat.data()), size, Type<T>::value, root, *comm);
    return mat;
  } else {
    Vector<T> result = Vector<T>::Zero(size);
    MPI_Bcast(result.data(), size, Type<T>::value, root, *comm);
    return result;
  }
}

template <class MATRIX>
typename std::enable_if<std::is_same<Matrix<typename MATRIX::Scalar>, MATRIX>::value, MATRIX>::type
broadcast(Communicator const &comm, t_uint root) {
  assert(root < comm.size());
  assert(root != comm.rank());
  auto const nrows = broadcast<t_uint>(comm, root);
  auto const ncols = broadcast<t_uint>(comm, root);
  MATRIX result = Matrix<typename MATRIX::Scalar>::Zero(nrows, ncols);
  MPI_Bcast(result.data(), nrows * ncols, Type<typename MATRIX::Scalar>::value, root, *comm);
  return result;
}

template <class VECTOR>
typename std::enable_if<std::is_same<Vector<typename VECTOR::Scalar>, VECTOR>::value, VECTOR>::type
broadcast(Communicator const &comm, t_uint root) {
  assert(root < comm.size());
  assert(root != comm.rank());
  auto const size = broadcast<t_uint>(comm, root);
  VECTOR result = Vector<typename VECTOR::Scalar>::Zero(size);
  MPI_Bcast(result.data(), size, Type<typename VECTOR::Scalar>::value, root, *comm);
  return result;
}

template <class T>
typename std::enable_if<is_registered_type<T>::value, std::vector<T>>::type
gather(T const &value, Communicator const &comm, t_uint root) {
  assert(root < comm.size());
  std::vector<T> result(root == comm.rank() ? comm.size() : 1);
  MPI_Gather(&value, 1, registered_type(value), result.data(), 1, registered_type(value), root,
             *comm);
  if(comm.rank() != root)
    result.clear();
  return result;
}

template <class T>
typename std::enable_if<is_registered_type<T>::value, std::vector<T>>::type
all_gather(T const &value, Communicator const &comm) {
  std::vector<T> result(comm.size());
  MPI_Allgather(&value, 1, registered_type(value), result.data(), 1, registered_type(value), *comm);
  return result;
}

template <class T>
typename std::enable_if<is_registered_type<typename T::Scalar>::value,
                        Vector<typename T::Scalar>>::type
all_gather(Eigen::PlainObjectBase<T> const &input, Communicator const &comm) {
  auto const sizes = all_gather<int>(input.size(), comm);
  std::vector<int> displs{0};
  for(t_uint i(1); i < sizes.size(); ++i)
    displs.push_back(displs.back() + sizes[i - 1]);
  Vector<typename T::Scalar> result(std::accumulate(sizes.begin(), sizes.end(), 0));
  MPI_Allgatherv(input.data(), input.size(), Type<typename T::Scalar>::value, result.data(),
                 sizes.data(), displs.data(), Type<typename T::Scalar>::value, *comm);
  return result;
}

} /* optime::mpi */
} /* optimet */
#endif
#endif /* ifndef OPTIMET_MPI_COMMUNICATOR */
