#include "catch.hpp"
#include <iostream>
#include <numeric>

#include "Types.h"
#include "mpi/Collectives.h"
#include "mpi/Communicator.h"
#include "scalapack/Belos.h"
#include "scalapack/InitExit.h"
#include "scalapack/LinearSystemSolver.h"

using namespace optimet;

template <class SCALAR>
void check_map(scalapack::Sizes const &grid, scalapack::Sizes const &size,
               scalapack::Sizes const &blocks) {
  if(mpi::Communicator().size() < grid.rows * grid.cols) {
    WARN("Not enough processes to run test: " << grid.rows << "x" << grid.cols << " < "
                                              << mpi::Communicator().size());
    return;
  }
  scalapack::Context const context(grid.rows, grid.cols);
  scalapack::Matrix<SCALAR> A(context, size, blocks);
  auto const split = mpi::Communicator().split(context.is_valid());
  if(not context.is_valid())
    return;
  auto const map = matrix_map(A, split);
  CHECK(map->isContiguous());
  CHECK(map->getGlobalNumElements() == A.size());
  CHECK(map->getMinLocalIndex() == 0);
  CHECK(map->getMaxLocalIndex() + 1 == A.local().size());
  // Check global starting indices are different
  auto const indices = split.all_gather(map->getGlobalElement(0));
  std::set<int> const index_set(indices.begin(), indices.end());
  CHECK(index_set.size() == context.size());
}

TEST_CASE("Map a matrix") {
  SECTION("1x1") { check_map<double>({1, 1}, {1024, 2048}, {32, 64}); }
  SECTION("1x2") { check_map<double>({1, 2}, {1024, 2048}, {32, 64}); }
  SECTION("2x1") { check_map<double>({2, 1}, {1024, 2048}, {32, 64}); }
  SECTION("2x3") { check_map<double>({2, 3}, {1024, 2048}, {32, 64}); }
}

template <class SCALAR> void check_gmres(scalapack::Sizes const &grid, t_uint n, t_uint nb) {
  scalapack::Sizes const size = {n, n};
  scalapack::Sizes const blocks = {nb, nb};
  if(mpi::Communicator().size() < grid.rows * grid.cols) {
    WARN("Not enough processes to run test: " << grid.rows << "x" << grid.cols << " < "
                                              << mpi::Communicator().size());
    return;
  }
  scalapack::Context const parallel_context(grid.rows, grid.cols);
  scalapack::Context const serial_context(1, 1);
  scalapack::Matrix<SCALAR> Aserial(serial_context, size, blocks);
  scalapack::Matrix<SCALAR> bserial(serial_context, {size.cols, 1}, blocks);
  if(serial_context.is_valid()) {
    Aserial.local() = Matrix<SCALAR>::Random(size.rows, size.cols) -
                      0.5 * Matrix<SCALAR>::Ones(size.rows, size.cols) +
                      10e0 * Matrix<SCALAR>::Identity(size.rows, size.cols);
    bserial.local() = Matrix<SCALAR>::Random(size.rows, 1).array() - 0.5;
  }

  auto const Aparallel = Aserial.transfer_to(parallel_context);
  auto const bparallel = bserial.transfer_to(parallel_context);

  auto const result = gmres_linear_system(Aparallel, bparallel);
  REQUIRE(std::get<1>(result) == 0);
  auto const x = std::get<0>(result).transfer_to(serial_context);
  if(serial_context.is_valid())
    CHECK((Aserial.local() * x.local()).isApprox(bserial.local(), 1e-8));
}

TEST_CASE("Compute solver in parallel, check result in serial") {
  SECTION("double") {
    SECTION("1x1") { check_gmres<double>({1, 1}, 1024, 64); }
    SECTION("1x3") { check_gmres<double>({1, 3}, 1024, 64); }
    SECTION("3x1") { check_gmres<double>({3, 1}, 1024, 64); }
    SECTION("3x2") { check_gmres<double>({3, 2}, 1024, 64); }
  }
  SECTION("complex double") {
    SECTION("1x1") { check_gmres<std::complex<double>>({1, 1}, 1024, 64); }
    SECTION("1x3") { check_gmres<std::complex<double>>({1, 3}, 1024, 64); }
    SECTION("3x1") { check_gmres<std::complex<double>>({3, 1}, 1024, 64); }
    SECTION("3x2") { check_gmres<std::complex<double>>({3, 2}, 1024, 64); }
  }
}
