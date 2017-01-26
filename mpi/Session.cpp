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

#include "mpi/Session.h"
#include <exception>
#include <mpi.h>
#ifdef OPTIMET_BELOS
#include <Tpetra_Core.hpp>
#endif

namespace optimet {
namespace mpi {
namespace {
bool &did_done_do_init() {
  static bool nrefs = false;
  return nrefs;
}

t_uint &global_reference() {
  static t_uint nrefs = 0;
  return nrefs;
}
} // anonymous namespace

void init(int argc, const char **argv) {
  if(did_done_do_init())
    return;
#ifdef OPTIMET_BELOS
  Tpetra::initialize(&argc, const_cast<char ***>(&argv));
#else
  MPI_Init(&argc, const_cast<char ***>(&argv));
#endif
  did_done_do_init() = true;
}

bool initialized() { return did_done_do_init(); }

bool finalized() {
  int finalized;
  MPI_Finalized(&finalized);
  return finalized;
}

void finalize() {
  if(finalized() or not initialized())
    return;
#ifdef OPTIMET_BELOS
  Kokkos::finalize();
# endif
  MPI_Finalize();
}

void increment_ref() { ++global_reference(); }
void decrement_ref() { --global_reference(); }
} /* optimet::mpi  */
} /* optimet  */
