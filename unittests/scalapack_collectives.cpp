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
