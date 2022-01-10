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

#ifndef OPTIMET_PARALLEL_PARAMETERS_H
#define OPTIMET_PARALLEL_PARAMETERS_H

#include "Types.h"

namespace optimet {
namespace scalapack {
//! Rows and Colums of the local blocks
struct Sizes {
  t_uint rows, cols;
};
//! Indices of the process starting the distribution
struct Index {
  t_uint row, col;
};

//! Parameters needed to setup parallel computations
struct Parameters {
  t_uint block_size;
  Sizes grid;

  Parameters(t_uint block_size = 64, Sizes grid = {0, 0})
      : block_size(block_size), grid(grid) {}
};

#ifdef OPTIMET_SCALAPACK
//! \brief tries and determines the squarest and largest grid size for a given number of procs
//! \param[in] nprocs: number of procs to parcel out into a grid
//! \parma[out] skew: determines how much to prefer square (skew large) to rectangular (skew small)
//! grids.
//! \details The cost function is n * m - skew * (n - m) * (n - m). Smaller skews result in grids
//! that will include more processes in the grid at the expense of a less square grid.
Sizes squarest_largest_grid(t_uint nprocs, t_real skew=0.1);
#endif
} // scalapack
} // optimet
#endif
