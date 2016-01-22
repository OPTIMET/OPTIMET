#ifndef OPTIMET_BLACS_EXIT_H
#define OPTIMET_BLACS_EXIT_H

#include "Types.h"

namespace optimet {
namespace mpi {
//! Calls mpi init
void init(int argc, char **argv);
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
} /* optimet::mpi */
} /* optimet */
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
