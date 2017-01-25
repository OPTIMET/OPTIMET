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

#include <exception>
#include <mpi.h>
#include "mpi/Communicator.h"
#include "mpi/Session.h"

namespace optimet { namespace mpi {

void Communicator::delete_comm(Communicator::Impl *const impl) {
  if(impl->comm != MPI_COMM_WORLD and initialized() and not finalized())
    MPI_Comm_free(&impl->comm);
  decrement_ref();
  delete impl;
}

Communicator::Communicator(MPI_Comm const& comm) : impl(nullptr) {
  reset(&comm);
}

void Communicator::reset(MPI_Comm const * const comm) {
  if(comm == nullptr) {
    impl.reset();
    return;
  }
  if(not initialized())
    throw std::runtime_error("Mpi was not initialized");

  int size, rank;
  MPI_Comm_size(*comm, &size);
  MPI_Comm_rank(*comm, &rank);

  Impl const data{*comm, static_cast<t_uint>(size), static_cast<t_uint>(rank)};
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
