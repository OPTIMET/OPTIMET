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
