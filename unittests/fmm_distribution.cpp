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

#include "catch.hpp"
#include "mpi/FastMatrixMultiply.h"
#include <iostream>

TEST_CASE("Column and row distributions") {
  using namespace optimet;
  SECTION("Column distributions") {
    using optimet::mpi::details::vector_distribution;
    CHECK(vector_distribution(0, 4).size() == 0);

    CHECK(vector_distribution(1, 4) == Vector<t_int>::Zero(1));

    CHECK(vector_distribution(2, 4)(0) == 0);
    CHECK(vector_distribution(2, 4)(1) == 1);

    CHECK(vector_distribution(4, 4)(0) == 0);
    CHECK(vector_distribution(4, 4)(1) == 1);
    CHECK(vector_distribution(4, 4)(2) == 2);
    CHECK(vector_distribution(4, 4)(3) == 3);

    CHECK(vector_distribution(5, 4).head(2) == Vector<t_int>::Zero(2));
    CHECK(vector_distribution(5, 4)(2) == 1);
    CHECK(vector_distribution(5, 4)(3) == 2);
    CHECK(vector_distribution(5, 4)(4) == 3);

    CHECK(vector_distribution(5, 3).head(2) == Vector<t_int>::Zero(2));
    CHECK(vector_distribution(5, 3).segment(2, 2) == Vector<t_int>::Ones(2));
    CHECK(vector_distribution(5, 3)(4) == 2);

    CHECK(vector_distribution(5, 2).head(3) == Vector<t_int>::Zero(3));
    CHECK(vector_distribution(5, 2).segment(3, 2) == Vector<t_int>::Ones(2));
  }

  SECTION("Local vs non-local") {
    using optimet::mpi::details::local_interactions;
    CHECK(local_interactions(0).rows() == 0);
    CHECK(local_interactions(0).cols() == 0);

    CHECK(local_interactions(1).all());
    CHECK(local_interactions(2).all());

    CHECK(local_interactions(3).cast<int>().sum() == 7);
    CHECK(local_interactions(3)(2, 0) == false);
    CHECK(local_interactions(3)(0, 2) == false);

    CHECK(local_interactions(3, 0).cast<int>().sum() == 3);
    CHECK(local_interactions(3, 0).diagonal().cast<int>().sum() == 3);

    CHECK(local_interactions(13, 2) == local_interactions(13, 2).transpose());
  }
}

TEST_CASE("Graph communicators creation") {
  using namespace optimet;
  using mpi::details::vector_distribution;
  using mpi::details::graph_edges;

  auto const nprocs = 3;
  auto const nscatt = 6;
  auto const vecdist = vector_distribution(nscatt, nprocs);
  SECTION("Nothing allowed") {
    CHECK(graph_edges(Matrix<bool>::Zero(nscatt, nscatt), vecdist).size() == nprocs);
    for(auto const &connections : graph_edges(Matrix<bool>::Zero(nscatt, nscatt), vecdist))
      CHECK(connections.size() == 0);
  }

  SECTION("All allowed") {
    auto const locals = graph_edges(Matrix<bool>::Ones(nscatt, nscatt), vecdist);
    CHECK(locals.size() == nprocs);
    for(t_int i(0); i < locals.size(); ++i) {
      CHECK(locals[i].size() == nprocs);
      for(t_uint n(0); n < nprocs; ++n)
        CHECK((locals[i].count(n) == 1));
    }
  }

  SECTION("Quarter") {
    Matrix<bool> allowed = Matrix<bool>::Zero(nscatt, nscatt);
    allowed.topLeftCorner(nscatt / 2, nscatt / 2).fill(true);
    auto const locals = graph_edges(allowed, vecdist);
    CHECK(locals.size() == nprocs);
    CHECK(locals[0].size() == 2);
    CHECK(locals[0].count(0) == 1);
    CHECK(locals[0].count(1) == 1);
    CHECK(locals[1].size() == 2);
    CHECK(locals[1].count(0) == 1);
    CHECK(locals[1].count(1) == 1);
    CHECK(locals[2].size() == 0);
  }
}
