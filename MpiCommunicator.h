#ifndef OPTIMET_MPI_COMMUNICATOR_H
#define OPTIMET_MPI_COMMUNICATOR_H

#include <mpi.h>
#include "Types.h"
#include <memory>

namespace optimet {
namespace mpi {

//! \brief A C++ wrapper for an mpi communicator
//! \details All copies made of this communicator are shallow: they reference the same communicator.
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
  bool is_root() const { return rank() == 0; }

  //! \brief Duplicates this communicator
  //! \details Creates a new communicator in all ways equivalent to this one.
  Communicator duplicate() const;
  //! Alias for duplicate
  Communicator clone() const { return duplicate(); }

private:
  //! Holds data associated with the context
  std::shared_ptr<Impl const> impl;

  //! Deletes an mpi communicator
  static void delete_comm(Impl *impl);

  //! \brief Constructs a communicator
  //! \details Takes ownership of the communicator, unless it is MPI_COMM_WORLD.
  //! This means that once all the shared pointer to the impl are delete, the communicator will be
  //! released.
  Communicator(MPI_Comm const &comm);
};

} /* optime::mpi */
} /* optimet */
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
