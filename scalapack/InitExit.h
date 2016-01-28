#ifndef OPTIMET_SCALAPACK_EXIT_H
#define OPTIMET_SCALAPACK_EXIT_H

#include "Types.h"

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

//! Rows and Colums of the local blocks
struct Sizes {
  t_uint rows, cols;
};
//! Indices of the process starting the distribution
struct Index {
  t_uint row, col;
};

} /* scalapack */
} /* optime */
#endif /* ifndef OPTIMET_BLACS_CONTEXT */
