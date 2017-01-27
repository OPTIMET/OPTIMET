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

#ifndef OPTIMET_MPI_COLLECTIVES_HPP
#define OPTIMET_MPI_COLLECTIVES_HPP

#include "Types.h"
#ifdef OPTIMET_MPI
#include "RegisteredTypes.h"
#include <mpi.h>
#include <type_traits>
#include <vector>

namespace optimet {
namespace mpi {
class Communicator;
//! Syncs all procs
void barrier(Communicator const &com);
//! Broadcast from somewhere to somewhere
template <class T>
typename std::enable_if<is_registered_type<T>::value, T>::type
broadcast(T const &, Communicator const &, t_uint root);
//! Broadcast from somewhere to somewhere
template <class T>
typename std::enable_if<is_registered_type<T>::value, T>::type
broadcast(Communicator const &, t_uint root);

//! Gathers data from all procs to root
template <class T>
typename std::enable_if<is_registered_type<T>::value, std::vector<T>>::type
gather(T const &, Communicator const &, t_uint root);
//! Gathers data from all procs to all procs
template <class T>
typename std::enable_if<is_registered_type<T>::value, std::vector<T>>::type
all_gather(T const &, Communicator const &);
template <class T>
typename std::enable_if<is_registered_type<typename T::Scalar>::value,
                        Vector<typename T::Scalar>>::type
all_gather(Eigen::PlainObjectBase<T> const &input, Communicator const &comm);

//! Broadcasts an eigen matrix
template <class T> Matrix<T> broadcast(Matrix<T> const &, Communicator const &, t_uint);
//! Broadcasts an eigen vector
template <class T> Vector<T> broadcast(Vector<T> const &, Communicator const &, t_uint);
//! Broadcasts an eigen matrix
template <class MATRIX>
typename std::enable_if<std::is_same<Matrix<typename MATRIX::Scalar>, MATRIX>::value, MATRIX>::type
broadcast(Communicator const &comm, t_uint root);
//! Broadcasts an eigen vector
template <class VECTOR>
typename std::enable_if<std::is_same<Vector<typename VECTOR::Scalar>, VECTOR>::value, VECTOR>::type
broadcast(Communicator const &comm, t_uint root);
} /* optime::mpi */
} /* optimet */
#endif
#endif /* ifndef OPTIMET_MPI_COMMUNICATOR */
