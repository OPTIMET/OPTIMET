#include <exception>
#include <mpi.h>
#include "MpiCommunicator.h"
#include "MpiExit.h"

namespace optimet {

void MpiCommunicator::delete_comm(MpiCommunicator::Impl *const impl) {
  if(impl->comm != MPI_COMM_WORLD and mpi_initialized() and not mpi_finalized())
    MPI_Comm_free(&impl->comm);
  decrement_mpi_ref();
  delete impl;
}

MpiCommunicator::MpiCommunicator(MPI_Comm const& comm) : impl(nullptr) {
  if(not mpi_initialized())
    throw std::runtime_error("Mpi was not initialized");

  int size, rank;
  MPI_Comm_size(comm, &size);
  MPI_Comm_rank(comm, &rank);

  Impl const data{comm, static_cast<t_uint>(size), static_cast<t_uint>(rank)};
  impl = std::shared_ptr<Impl const>(new Impl(data), &delete_comm);
  if(impl)
    increment_mpi_ref();
}

} /* optimet  */
