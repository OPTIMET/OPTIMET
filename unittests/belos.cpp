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
#include <iostream>
#include <numeric>

#include "Reader.h"
#include "Run.h"
#include "Types.h"
#include "mpi/Collectives.h"
#include "mpi/Communicator.h"
#include "scalapack/Belos.h"
#include "scalapack/InitExit.h"
#include "scalapack/LinearSystemSolver.h"
#include <string>

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

template <class SCALAR>
void check_vector(scalapack::Sizes const &grid, scalapack::Sizes const &size,
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
  A.local() = Matrix<SCALAR>::Ones(size.rows, size.cols);
  auto const vector = tpetra_vector(A, split);
  CHECK(vector->getNumVectors() == 1);
  CHECK(vector->getGlobalLength() == A.size());
  CHECK(vector->getLocalLength() == A.local().size());

  CHECK(vector->getVector(0)->norm1() == A.size());
  A.local() *= 2;
  CHECK(vector->getVector(0)->norm1() == A.size());
}

TEST_CASE("Tpetra vector from a matrix") {
  SECTION("1x1") { check_vector<double>({1, 1}, {1024, 2048}, {32, 64}); }
  SECTION("1x2") { check_vector<double>({1, 2}, {1024, 2048}, {32, 64}); }
  SECTION("2x1") { check_vector<double>({2, 1}, {1024, 2048}, {32, 64}); }
  SECTION("2x3") { check_vector<double>({2, 3}, {1024, 2048}, {32, 64}); }
}

template <class SCALAR>
void check_linear_problem(scalapack::Sizes const &grid, t_uint const &size, t_uint const &blocks) {
  if(mpi::Communicator().size() < grid.rows * grid.cols) {
    WARN("Not enough processes to run test: " << grid.rows << "x" << grid.cols << " < "
                                              << mpi::Communicator().size());
    return;
  }
  scalapack::Context const context(grid.rows, grid.cols);
  scalapack::Matrix<SCALAR> A(context, {size, size}, {blocks, blocks});
  scalapack::Matrix<SCALAR> b(context, {size, 1}, A.blocks());
  scalapack::Matrix<SCALAR> expected(context, b.sizes(), b.blocks());
  auto const split = mpi::Communicator().split(context.is_valid());
  if(not context.is_valid())
    return;

  if(A.local().size() > 0)
    A.local() = Matrix<SCALAR>::Random(A.local().rows(), A.local().cols());
  if(b.local().size() > 0)
    b.local() = Matrix<SCALAR>::Random(b.local().rows(), b.local().cols());

  Teuchos::RCP<scalapack::BelosOperator<SCALAR>> Aptr =
      Teuchos::rcp(new scalapack::BelosOperator<SCALAR>(A));
  auto rhs = tpetra_vector(b, split);
  auto lhs = tpetra_vector(b, split);
  // auto belos_residuals = tpetra_vector(residuals, split);
  Teuchos::RCP<scalapack::BelosLinearProblem<SCALAR>> problem =
      rcp(new scalapack::BelosLinearProblem<SCALAR>(Aptr, lhs, rhs));
  CHECK(problem->setProblem());
  problem->apply(*lhs, *rhs);

  pdgemm(1e0, A, b, 0e0, expected);
  auto const actual = scalapack::as_matrix(*rhs, A);
  CHECK(actual.local().isApprox(expected.local()));
}

TEST_CASE("Create and apply belos linear problem") {
  SECTION("double 1x1") { check_linear_problem<double>({1, 1}, 1024, 64); }
  SECTION("double 1x2") { check_linear_problem<double>({1, 2}, 1024, 64); }
  SECTION("double 2x1") { check_linear_problem<double>({2, 1}, 1024, 64); }
  SECTION("double 2x3") { check_linear_problem<double>({2, 3}, 1024, 64); }
  SECTION("complex double 1x1") { check_linear_problem<std::complex<double>>({1, 1}, 1024, 64); }
  SECTION("complex double 1x2") { check_linear_problem<std::complex<double>>({1, 2}, 1024, 64); }
  SECTION("complex double 2x1") { check_linear_problem<std::complex<double>>({2, 1}, 1024, 64); }
  SECTION("complex double 2x3") { check_linear_problem<std::complex<double>>({2, 3}, 1024, 64); }
}

template <class SCALAR> void check_gmres(scalapack::Sizes const &grid, t_uint n, t_uint nb) {
  scalapack::Sizes const size = {n, n};
  scalapack::Sizes const blocks = {nb, nb};
  if(mpi::Communicator().size() < grid.rows * grid.cols) {
    WARN("Not enough processes to run test: " << grid.rows << "x" << grid.cols << " < "
                                              << mpi::Communicator().size());
    return;
  }
  scalapack::Context const context(grid.rows, grid.cols);
  scalapack::Matrix<SCALAR> b(context, {size.cols, 1}, blocks);

  // create a simple matrix to diagonalize
  scalapack::Context const serial_context(1, 1);
  scalapack::Matrix<SCALAR> Aserial(serial_context, size, blocks);
  if(serial_context.is_valid())
    Aserial.local() =
        (Matrix<SCALAR>::Random(Aserial.local().rows(), Aserial.local().cols()).array() - 0.5)
            .matrix() +
        Matrix<SCALAR>::Identity(Aserial.local().rows(), Aserial.local().cols()) * 5;
  auto const A = Aserial.transfer_to(context);
  b.local() = Matrix<SCALAR>::Random(b.local().rows(), b.local().cols()).array() - 0.5;

  auto params = Teuchos::rcp(new Teuchos::ParameterList);
  params->set("Solver", "GMRES");
  params->set<int>("Num Blocks", Aserial.rows());
  params->set("Maximum Iterations", 4000);
  params->set("Convergence Tolerance", 1.0e-14);
  auto const result = gmres_linear_system(A, b, params);
  REQUIRE(std::get<1>(result) == 0);
  if(context.is_valid()) {
    scalapack::Matrix<SCALAR> x = b;
    scalapack::pdgemm(1e0, A, std::get<0>(result), -1e0, x);
    CHECK(x.local().isZero(1e-6));
  }
}

TEST_CASE("Compute solver in parallel, check result in serial") {
  SECTION("double") {
    SECTION("1x1") { check_gmres<double>({1, 1}, 500, 64); }
    SECTION("2x1") { check_gmres<double>({2, 1}, 500, 64); }
    SECTION("3x2") { check_gmres<double>({3, 2}, 1024, 64); }
  }
  SECTION("complex double") {
    SECTION("1x1") { check_gmres<std::complex<double>>({1, 1}, 560, 64); }
    SECTION("3x1") { check_gmres<std::complex<double>>({3, 1}, 560, 64); }
    SECTION("3x2") { check_gmres<std::complex<double>>({3, 2}, 560, 64); }
  }
}

TEST_CASE("Read XML") {
  std::istringstream buffer(
      "<simulation>\n"
      "  <harmonics nmax=\"6\" />\n"
      "</simulation>\n"
      "<source type=\"planewave\">\n"
      "  <wavelength value=\"1460\" />\n"
      "  <propagation theta=\"90\" phi=\"90\" />\n"
      "  <polarization Etheta.real=\"1.0\" Etheta.imag=\"0.0\" Ephi.real=\"0.0\" "
      "Ephi.imag=\"0.0\" />\n"
      "</source>\n"
      "<geometry>\n"
      "  <object type=\"sphere\">\n"
      "    <cartesian x=\"0.0\" y=\"0.0\" z=\"0.0\" />\n"
      "    <properties radius=\"500.0\" />\n"
      "    <epsilon type=\"relative\" value.real=\"13.0\" value.imag=\"0.0\" />\n"
      "    <mu type=\"relative\" value.real=\"1.0\" value.imag=\"0.0\" />\n"
      "  </object>\n"
      "</geometry>\n"
      "<output type=\"field\">\n"
      "  <grid type=\"cartesian\">\n"
      "    <x min=\"-1000\" max=\"1000\" steps=\"21\" />\n"
      "    <y min=\"-0.1\" max=\"0.1\" steps=\"2\" />\n"
      "    <z min=\"-1000\" max=\"1000\" steps=\"21\" />\n"
      "  </grid>\n"
      "</output>\n"
      "<ParameterList name=\"Belos\">\n"
      "  <Parameter name=\"Solver\" type=\"string\" value=\"GMRES\"/>\n"
      "  <Parameter name=\"Maximum Iterations\" type=\"int\" value=\"4000\"/>\n"
      "  <Parameter name=\"Output Frequency\" type=\"int\" value=\"20\"/>\n"
      "</ParameterList>\n");

  auto const run = optimet::simulation_input(buffer);
  CHECK(run.belos_params->get<std::string>("Solver") == "GMRES");
  CHECK(run.belos_params->get<int>("Maximum Iterations") == 4000);
  CHECK(run.belos_params->get<int>("Output Frequency") == 20);
}
