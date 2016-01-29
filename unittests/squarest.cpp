#include <iostream>
#include "catch.hpp"

#include "Types.h"
#include "scalapack/Parameters.h"

using namespace optimet;

TEST_CASE("Check squarest") {
  SECTION("n by n") {
    for(t_uint i(0); i < 10; ++i) {
      CHECK(scalapack::squarest_largest_grid(i * i).rows == i);
      CHECK(scalapack::squarest_largest_grid(i * i).cols == i);
    }
  }
  // This is mostly for us to figure out whether the cost function is reasonable
  // Only actual benchmarks can reliably give an answer
  SECTION("n by m") {
    std::tuple<t_uint, t_uint, t_uint> const fixtures[] = {
      // number of procs, number of row, number of cols
      {2, 1, 2}, {3, 1, 3}, {4, 2, 2}, {5, 2, 2}, {6, 2, 3}, {7, 2, 3}, {8, 2, 4}, {9, 3, 3},
      {10, 3, 3}, {11, 3, 3}, {12, 3, 4}, {13, 3, 4}, {14, 3, 4}, {15, 3, 5}, {16, 4, 4},
      {17, 4, 4}, {18, 4, 4}, {19, 4, 4}, {20, 4, 5}, {21, 4, 5}, {22, 4, 5}, {23, 4, 5},
      {24, 4, 6}
    };
    for(auto const fixture: fixtures) {
      CHECK(scalapack::squarest_largest_grid(std::get<0>(fixture)).rows == std::get<1>(fixture));
      CHECK(scalapack::squarest_largest_grid(std::get<0>(fixture)).cols == std::get<2>(fixture));
    }
  }
}
