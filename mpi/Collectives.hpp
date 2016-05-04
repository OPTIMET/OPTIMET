#ifndef OPTIMET_MPI_COLLECTIVES_HPP
#define OPTIMET_MPI_COLLECTIVES_HPP

#include "Types.h"
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
typename std::enable_if<std::is_fundamental<T>::value, T>::type
broadcast(T const &, Communicator const &, t_uint root);
//! Broadcast from somewhere to somewhere
template <class T>
typename std::enable_if<std::is_fundamental<T>::value, T>::type
broadcast(Communicator const &, t_uint root);

//! Gathers data from all procs to root
template <class T>
typename std::enable_if<std::is_fundamental<T>::value, std::vector<T>>::type
gather(T const &, Communicator const &, t_uint root);
//! Gathers data from all procs to all procs
template <class T>
typename std::enable_if<std::is_fundamental<T>::value, std::vector<T>>::type
all_gather(T const &, Communicator const &);

//! Broadcasts an eigen matrix
template <class T> Matrix<T> broadcast(Matrix<T> const &, Communicator const &, t_uint);
//! Broadcasts an eigen matrix
template <class MATRIX>
typename std::enable_if<std::is_same<Matrix<typename MATRIX::Scalar>, MATRIX>::value, MATRIX>::type
broadcast(Communicator const &comm, t_uint root);
} /* optime::mpi */
} /* optimet */
#endif /* ifndef OPTIMET_MPI_COMMUNICATOR */
