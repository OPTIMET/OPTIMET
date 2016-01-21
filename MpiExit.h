#ifndef OPTIMET_BLACS_EXIT_H
#define OPTIMET_BLACS_EXIT_H

#include "Types.h"

namespace optimet {
//! Calls mpi init
void mpi_init(int argc, char **argv );
//! True if mpi has been initialized
bool mpi_initialized();
//! True if mpi has been finalized
bool mpi_finalized();
//! Closes mpi stuff
void mpi_finalize();
//! Increments number of mpi objects
void increment_mpi_ref();
//! Decrements number of mpi objects
void decrement_mpi_ref();
} /* optime */
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
