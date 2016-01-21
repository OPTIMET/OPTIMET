#ifndef OPTIMET_MPI_COMMUNICATOR_H
#define OPTIMET_MPI_COMMUNICATOR_H

#include <mpi.h>
#include "Types.h"
#include <memory>

namespace optimet {

//! A C++ wrapper for an mpi communicator
class MpiCommunicator {
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
  MpiCommunicator() : MpiCommunicator(MPI_COMM_WORLD) {};

  virtual ~MpiCommunicator() {};

  //! The number of processes
  decltype(Impl::size) size() const { return impl ? impl->size: 1; }
  //! The rank of this proc
  decltype(Impl::rank) rank() const { return impl ? impl->rank: 0; }
  //! Returns the Blacs context in a way blacs undersands
  decltype(Impl::comm) operator*() const { return impl->comm; }

private:
  //! Holds data associated with the context
  std::shared_ptr<Impl const> impl;

  //! Deletes an mpi communicator
  static void delete_comm(Impl *impl);

  //! Constructs a communicator
  MpiCommunicator(MPI_Comm const &comm);
};

} /* optime */
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
