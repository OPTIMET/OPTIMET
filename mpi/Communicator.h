#ifndef OPTIMET_MPI_COMMUNICATOR_H
#define OPTIMET_MPI_COMMUNICATOR_H

#include "Types.h"

#ifdef OPTIMET_MPI

#include "mpi/Collectives.hpp"
#include "mpi/RegisteredTypes.h"
#include <memory>
#include <mpi.h>
#include <type_traits>
#include <vector>

namespace optimet {
namespace mpi {

//! \brief A C++ wrapper for an mpi communicator
//! \details All copies made of this communicator are shallow: they reference
//! the same communicator.
class Communicator {
  //! Holds actual data associated with mpi
  struct Impl {
    //! The blacs context
    MPI_Comm comm;
    //! The number of processes
    t_uint size;
    //! The rank of this object
    t_uint rank;
  };

public:
  //! World communicator
  Communicator() : Communicator(MPI_COMM_WORLD){};

  virtual ~Communicator(){};

  //! The number of processes
  decltype(Impl::size) size() const { return impl ? impl->size : 1; }
  //! The rank of this proc
  decltype(Impl::rank) rank() const { return impl ? impl->rank : 0; }
  //! Returns the Blacs context in a way blacs undersands
  decltype(Impl::comm) operator*() const { return impl->comm; }

  //! Split current communicator
  Communicator split(t_int color) const { return split(color, rank()); }
  //! Split current communicator
  Communicator split(t_int color, t_uint rank) const;

  //! True if object is root
  bool is_root() const { return rank() == root_id(); }

  //! \brief Duplicates this communicator
  //! \details Creates a new communicator in all ways equivalent to this one.
  Communicator duplicate() const;
  //! Alias for duplicate
  Communicator clone() const { return duplicate(); }

  //! \brief Helper function for broadcasting
  template <class T>
  decltype(optimet::mpi::broadcast(std::declval<T>(), std::declval<Communicator>(), 0u))
  broadcast(T const &args, t_uint root = Communicator::root_id()) const {
    return optimet::mpi::broadcast(args, *this, static_cast<t_uint>(root));
  }
  //! Helper function for broadcasting
  template <class T>
  decltype(optimet::mpi::broadcast<T>(std::declval<Communicator>(), 0u))
  broadcast(t_uint root = Communicator::root_id()) const {
    return optimet::mpi::broadcast<T>(*this, root);
  }

  //! Helper function for gathering
  template <class T>
  decltype(optimet::mpi::gather(std::declval<T>(), std::declval<Communicator>(), 0u))
  gather(T const &value, t_uint root = Communicator::root_id()) const {
    return optimet::mpi::gather(value, *this, root);
  }
  //! Helper function for gathering
  template <class T>
  decltype(optimet::mpi::all_gather(std::declval<T>(), std::declval<Communicator>()))
  all_gather(T const &value) const {
    return optimet::mpi::all_gather(value, *this);
  }
  //! Helper function for gathering
  template <class T>
  typename std::enable_if<std::is_fundamental<T>::value, T>::type
  all_reduce(T const &value, MPI_Op operation) const {
    T result;
    MPI_Allreduce(&value, &result, 1, registered_type(value), operation, **this);
    return result;
  }

  void barrier() const { return optimet::mpi::barrier(*this); }

  //! Root id for this communicator
  static constexpr t_uint root_id() { return 0; }

private:
  //! Holds data associated with the context
  std::shared_ptr<Impl const> impl;

  //! Deletes an mpi communicator
  static void delete_comm(Impl *impl);

  //! \brief Constructs a communicator
  //! \details Takes ownership of the communicator, unless it is MPI_COMM_WORLD.
  //! This means that once all the shared pointer to the impl are delete, the
  //! communicator will be
  //! released.
  Communicator(MPI_Comm const &comm);
};

} /* optime::mpi */
} /* optimet */
#else
namespace optimet {
namespace mpi {
class Communicator {
public:
  constexpr t_uint size() const { return 1; }
  constexpr t_uint rank() const { return 0; }
  constexpr t_uint root_id() const { return 0; }
};
}
}
#endif /* ifdef OPTIMET_MPI */
#endif /* ifndef OPTIMET_MPI_COMMUNICATOR */
