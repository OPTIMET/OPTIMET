#ifndef OPTIMET_SCALAPACK_EXIT_H
#define OPTIMET_SCALAPACK_EXIT_H

#include "Types.h"

#ifdef OPTIMET_MPI
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
