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
    typedef std::tuple<t_uint, t_uint, t_uint> T;
    T const fixtures[] = {
      // number of procs, number of row, number of cols
      T{2, 1, 2}, T{3, 1, 3}, T{4, 2, 2}, T{5, 2, 2}, T{6, 2, 3}, T{7, 2, 3}, T{8, 2, 4},
      T{9, 3, 3}, T{10, 3, 3}, T{11, 3, 3}, T{12, 3, 4}, T{13, 3, 4}, T{14, 3, 4}, T{15, 3, 5},
      T{16, 4, 4}, T{17, 4, 4}, T{18, 4, 4}, T{19, 4, 4}, T{20, 4, 5}, T{21, 4, 5}, T{22, 4, 5},
      T{23, 4, 5}, T{24, 4, 6}
    };
    for(auto const fixture: fixtures) {
      CHECK(scalapack::squarest_largest_grid(std::get<0>(fixture)).rows == std::get<1>(fixture));
      CHECK(scalapack::squarest_largest_grid(std::get<0>(fixture)).cols == std::get<2>(fixture));
    }
  }
}
