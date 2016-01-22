#include <exception>
#include <mpi.h>
#include "MpiCommunicator.h"
#include "MpiExit.h"

namespace optimet { namespace mpi {

void Communicator::delete_comm(Communicator::Impl *const impl) {
  if(impl->comm != MPI_COMM_WORLD and initialized() and not finalized())
    MPI_Comm_free(&impl->comm);
  decrement_ref();
  delete impl;
}

Communicator::Communicator(MPI_Comm const& comm) : impl(nullptr) {
  if(not initialized())
    throw std::runtime_error("Mpi was not initialized");

  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  Impl const data{comm, static_cast<t_uint>(size), static_cast<t_uint>(rank)};
  impl = std::shared_ptr<Impl const>(new Impl(data), &delete_comm);
  if(impl)
    increment_ref();
}

Communicator Communicator::split(t_int color, t_uint rank) const {
  MPI_Comm comm;
  MPI_Comm_split(**this, color, static_cast<t_int>(rank), &comm);
  return comm;
}

Communicator Communicator::duplicate() const {
  MPI_Comm comm;
  MPI_Comm_dup(**this, &comm);
  return comm;
}

} /* optimet::mpi */
} /* optimet  */
