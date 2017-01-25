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
#include "scalapack/Collectives.h"
#include "scalapack/Context.h"

using namespace optimet;

TEST_CASE("broadcasting") {
  auto const world = scalapack::Context::Squarest();

  SECTION("integer") {
    for(t_int i(0); i < world.rows(); ++i)
      for(t_int j(0); j < world.cols(); ++j) {
        auto const value = world.row() == i and world.col() == j ?
                               world.broadcast(2) :
                               world.broadcast<int>(static_cast<t_uint>(i), static_cast<t_uint>(j));
        if(world.is_valid()) {
          CHECK(value == 2);
        }
      }
  }

  SECTION("complex") {
    for(t_int i(0); i < world.rows(); ++i)
      for(t_int j(0); j < world.cols(); ++j) {
        auto const value = world.broadcast<t_complex>(t_complex(world.row(), world.col()), i, j);
        if(world.is_valid()) {
          CHECK(value.real() == Approx(i));
          CHECK(value.imag() == Approx(j));
        }
      }
  }

  SECTION("matrix") {
    auto matrix = [](t_uint i, t_uint j) -> Matrix<t_real> {
      return Matrix<t_real>::Ones(1 + i, 3 + j) * (i * j);
    };
    for(t_int i(0); i < world.rows(); ++i)
      for(t_int j(0); j < world.cols(); ++j) {
        auto const value = world.broadcast(matrix(world.row(), world.col()), i, j);
        if(world.is_valid())
          CHECK(value.isApprox(matrix(i, j)));
      }
  }
}
