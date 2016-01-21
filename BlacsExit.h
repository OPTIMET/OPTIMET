#ifndef OPTIMET_BLACS_EXIT_H
#define OPTIMET_BLACS_EXIT_H

#include "Types.h"

namespace optimet {
//! Exits all blacs stuff
void blacks_exit(t_int status);
//! Increments number of blacs objects
void increment_blacs_ref();
//! Decrements number of blacs objects
void decrement_blacs_ref();
} /* optime */
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
