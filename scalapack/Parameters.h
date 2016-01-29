#ifndef OPTIMET_PARALLEL_PARAMETERS
#define OPTIMET_PARALLEL_PARAMETERS

#include "Types.h"
#include "scalapack/Matrix.h"

namespace optimet {
namespace scalapack {
//! Parameters needed to setup parallel computations
struct Parameters {
  t_uint block_size;
  Sizes grid;

  Parameters(t_uint block_size = 64, Sizes grid = {0, 0})
      : block_size(block_size), grid(grid) {}
};

//! tries and determines the squarest and largest grid size for a given number of procs
Sizes squarest_largest_grid(t_uint n);
}
}
#endif
