#include <algorithm>
#include <vector>
#include <numeric>
#include "scalapack/Parameters.h"
#include "scalapack/Matrix.h"

namespace optimet {
namespace scalapack {

Sizes squarest_largest_grid(t_uint nsquared) {
  auto const judgment = [nsquared](t_uint i) {
    Sizes const s{i, nsquared / i};
    return s.rows * s.cols - 0.5 * (s.rows - s.cols) * (s.rows - s.cols);
  };
  auto const compare = [judgment](t_uint i, t_uint j) { return judgment(i) < judgment(j); };
  std::vector<t_uint> range(nsquared / 2);
  std::iota(range.begin(), range.end(), 1);
  auto const result = *std::max_element(range.begin(), range.end(), compare);
  return {result, nsquared / result};
}
} /* scalapac */
} /* optimet */
