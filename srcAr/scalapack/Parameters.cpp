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

#include <algorithm>
#include <vector>
#include <numeric>
#include <algorithm>
#include "scalapack/Parameters.h"
#include "scalapack/Matrix.h"

namespace optimet {
namespace scalapack {

Sizes squarest_largest_grid(t_uint nprocs, t_real skew) {
  if(nprocs == 0)
    return {0, 0};
  auto const judgment = [nprocs, &skew](t_uint i) {
    Sizes const s{i, std::max<t_uint>(1, nprocs / i)};
    return s.rows * s.cols - skew * (s.rows - s.cols) * (s.rows - s.cols);
  };
  auto const compare = [judgment](t_uint i, t_uint j) { return judgment(i) < judgment(j); };
  std::vector<t_uint> range(std::max<t_uint>(1, nprocs / 2));
  std::iota(range.begin(), range.end(), 1);
  auto const n0 = std::max_element(range.begin(), range.end(), compare) - range.begin() + 1;
  auto const n1 = std::max<t_uint>(1, nprocs / n0);
  return {std::min<t_uint>(n0, n1), std::max<t_uint>(n0, n1)};
}
} /* scalapac */
} /* optimet */
