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
#include "scalapack/Context.h"
#include "scalapack/InitExit.h"
#include "mpi/Communicator.h"
#include "mpi/Collectives.h"

using namespace optimet;

TEST_CASE("Check pinfo") {
  mpi::Communicator world;
  CHECK(world.size() == scalapack::global_size());

  auto ranks = world.gather(scalapack::global_rank());
  if(not world.is_root())
    return;

  REQUIRE(ranks.size() == world.size());
  std::sort(ranks.begin(), ranks.end(), std::less<t_uint>());
  for(std::size_t i(0); i < ranks.size(); ++i)
    CHECK(i == ranks[i]);
}

TEST_CASE("Creates a blacs context 1x1") {
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(1, 1);
  REQUIRE(context.is_valid() == (rank == 0));
  if(rank == 0) {
    CHECK(context.rows() == 1);
    CHECK(context.cols() == 1);
    CHECK(context.row() == 0);
    CHECK(context.col() == 0);
  }
}

void check_nxm(t_uint n, t_uint m) {
  if(mpi::Communicator().size() < m * n)
    return;
  auto const rank = scalapack::global_rank();
  scalapack::Context const context(n, m);
  REQUIRE(context.is_valid() == (rank < n * m));
  if(rank < n * m) {
    CHECK(context.rows() == n);
    CHECK(context.cols() == m);
    // Fortran is column major
    CAPTURE(rank);
    CAPTURE(m);
    CHECK(context.col() == rank % m);
    CHECK(context.row() == rank / m);
  }
}

TEST_CASE("Create different blacs context") {
  check_nxm(1, 1);
  check_nxm(2, 1);
  check_nxm(1, 2);
  check_nxm(2, 2);
  check_nxm(3, 2);
}

void check_matrix(scalapack::Context const &base, Matrix<t_uint> const &mat) {
  scalapack::Context context(base, mat);
  if(not base.is_valid())
    CHECK(not context.is_valid());
  auto const in_mat = std::find(mat.data(), mat.data() + mat.size(), scalapack::global_rank());
  CHECK(context.is_valid() == (in_mat != (mat.data() + mat.size())));
  if(not context.is_valid())
    return;
  REQUIRE(context.is_valid());
  REQUIRE(context.rows() == mat.rows());
  REQUIRE(context.cols() == mat.cols());
  CHECK(mat(context.row(), context.col()) == scalapack::global_rank());
}
void check_matrix(Matrix<t_uint> const &mat) {
  check_matrix(scalapack::Context(), mat);
}

TEST_CASE("Explicit organisation") {
  Matrix<t_uint> mat;
  SECTION("1x2") {
    if(scalapack::global_size() < 2)
      return;
    mat.resize(1, 2);
    mat << 1, 0;
    check_matrix(mat);
    if(scalapack::global_size() < 3)
      return;
    mat << 1, 2;
    check_matrix(mat);
  }
  SECTION("At least 4") {
    if(scalapack::global_size() < 4)
      return;
    mat.resize(2, 2);
    mat << 1, 2, 0, 3;
    check_matrix(mat);
  }
  SECTION("Sub-context") {
    // checks we can use one context to create another
    // this is more interesting if there are more than 4 procs
    if(scalapack::global_size() < 4)
      return;
    mat.resize(2, 2);
    mat << 1, 2, 0, 3;
    scalapack::Context base(mat);

    if(not base.is_valid())
      return;

    Matrix<t_uint> other;
    other.resize(2, 1);
    other << 1, 2;
    scalapack::Context context(base, other);
    bool is_valid = (scalapack::global_rank() == 0) or (scalapack::global_rank() == 2);
    CHECK(context.is_valid() == is_valid);
    if(context.is_valid()) {
      mat.resize(1, mat.size());
      CHECK(context.rows() == other.rows());
      CHECK(context.cols() == other.cols());
    }
  }
}

TEST_CASE("Splitting grid") {
  if(scalapack::global_size() < 4)
    return;

  SECTION("Split 1 by n") {
    scalapack::Context linear(1, scalapack::global_size());
    auto const split = linear.split(1, 2);
    CHECK(split.is_valid());
    CHECK(split.rows() == 1);
    if(linear.cols() % 2 == 0 or linear.col() >= linear.cols() / 2 + 1)
      CHECK(split.cols() == linear.cols() / 2);
    else
      CHECK(split.cols() == linear.cols() / 2 + 1);
  }

  SECTION("Split n by 1") {
    scalapack::Context linear(scalapack::global_size(), 1);
    auto const split = linear.split(2, 1);
    CHECK(split.is_valid());
    CHECK(split.cols() == 1);
    if(linear.rows() % 2 == 0 or linear.row() >= linear.rows() / 2 + 1)
      CHECK(split.rows() == linear.rows() / 2);
    else
      CHECK(split.rows() == linear.rows() / 2 + 1);
  }
}

TEST_CASE("Process number and coords") {
  scalapack::Context context(scalapack::squarest_largest_grid(scalapack::global_size()));
  mpi::Communicator split(mpi::Communicator().split(context.is_valid()));

  if(not context.is_valid())
    return;

  auto const rows = split.all_gather(context.row());
  auto const cols = split.all_gather(context.col());
  auto const pnum = split.all_gather(context.process_number());
  auto mpi_rank = [&rows, &cols](t_int row, t_int col) {
    for(std::size_t i(0); i < rows.size(); ++i)
      if(rows[i] == row and cols[i] == col)
        return i;
    return rows.size() + 1;
  };

  for(t_int i(0); i < context.rows(); ++i) {
    CHECK(std::find(rows.begin(), rows.end(), i) != rows.end());
    for(t_int j(0); j < context.cols(); ++j) {
      CHECK(std::find(cols.begin(), cols.end(), j) != cols.end());
      auto const r = mpi_rank(i, j);
      CHECK(pnum[r] == context.process_number(i, j));
      CHECK(i == context.process_coordinates(pnum[r]).row);
      CHECK(j == context.process_coordinates(pnum[r]).col);
    }
  }
}
