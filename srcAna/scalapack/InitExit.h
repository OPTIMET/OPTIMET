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

#ifndef OPTIMET_SCALAPACK_EXIT_H
#define OPTIMET_SCALAPACK_EXIT_H

#include "Types.h"

#ifdef OPTIMET_SCALAPACK
namespace optimet {
namespace scalapack {
//! Rank of this proc
t_uint global_rank();
//! Size of this proc
t_uint global_size();
//! Checks whether we have exited
bool finalized();
//! Exits all blacs stuff
void finalize(t_int status);
//! Increments number of blacs objects
void increment_ref();
//! Decrements number of blacs objects
void decrement_ref();
} /* scalapack */
} /* optime */
#endif
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
