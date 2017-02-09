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

#ifndef OPTIMET_BLACS_EXIT_H
#define OPTIMET_BLACS_EXIT_H

#include "Types.h"

#ifdef OPTIMET_MPI
#include <mpi.h>
#endif

namespace optimet {
namespace mpi {
#ifdef OPTIMET_MPI
//! Calls mpi init
void init(int argc, const char **argv);
//! True if mpi has been initialized
bool initialized();
//! True if mpi has been finalized
bool finalized();
//! Closes mpi stuff
void finalize();
//! Increments number of mpi objects
void increment_ref();
//! Decrements number of mpi objects
void decrement_ref();
#else
inline void init(int argc, const char **argv) {}
inline void finalize() {}
#endif
} /* optimet::mpi */
} /* optimet */
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
