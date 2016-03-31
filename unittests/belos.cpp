#include <iostream>
#include <numeric>
#include "catch.hpp"

#include "Types.h"
#include "scalapack/LinearSystemSolver.h"
#include "scalapack/InitExit.h"
#include "mpi/Communicator.h"
#include "mpi/Collectives.h"

using namespace optimet;

template<class SCALAR> void check(scalapack::Sizes const &grid, t_uint n, t_uint nb) {
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
    SECTION("1x1") { check<double>({1, 1}, 1024, 64); }
    SECTION("1x3") { check<double>({1, 3}, 1024, 64); }
    SECTION("3x1") { check<double>({3, 1}, 1024, 64); }
    SECTION("3x2") { check<double>({3, 2}, 1024, 64); }
  }
  SECTION("complex double") {
    SECTION("1x1") { check<std::complex<double>>({1, 1}, 1024, 64); }
    SECTION("1x3") { check<std::complex<double>>({1, 3}, 1024, 64); }
    SECTION("3x1") { check<std::complex<double>>({3, 1}, 1024, 64); }
    SECTION("3x2") { check<std::complex<double>>({3, 2}, 1024, 64); }
  }
}
