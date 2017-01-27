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

#include "HarmonicsIterator.h"
#include "PreconditionedMatrix.h"
#include "Tools.h"
#include "Types.h"
#include "constants.h"
#include "mpi/Communicator.h"
#include "scalapack/Matrix.h"

using namespace optimet;
constexpr t_real default_length() { return 2000e-9; }

Matrix<t_real> fcc_cell() {
  Matrix<t_real> result = Matrix<t_real>::Ones(3, 3) * 0.5;
  result.diagonal().fill(0);
  return result;
}

Scatterer default_scatterer(t_int nHarmonics) {
  auto const radius = 500e-9;
  ElectroMagnetic const elmag{13.1, 1.0};
  return {{0, 0, 0}, elmag, radius, nHarmonics};
}

std::tuple<Geometry, std::shared_ptr<Excitation>>
fcc_system(std::tuple<int, int, int> const &range, t_real length, Scatterer const &scatterer) {
  // setup geometry
  auto const cell = fcc_cell();
  Geometry geometry;
  for(int i(0); i < std::get<0>(range); ++i)
    for(int j(0); j < std::get<1>(range); ++j)
      for(int k(0); k < std::get<2>(range); ++k) {
        Eigen::Matrix<t_real, 3, 1> pos = cell * Eigen::Matrix<t_real, 3, 1>(i, j, k) * length;
        geometry.pushObject({Tools::toSpherical({pos(0), pos(1), pos(2)}), scatterer.elmag,
                             scatterer.radius, scatterer.nMax});
      }

  // Create excitation
  auto const wavelength = 750e-9;
  Spherical<t_real> const vKinc{2 * consPi / wavelength, 90 * consPi / 180.0, 90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto const excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, scatterer.nMax);
  excitation->populate();
  geometry.update(excitation);
  return std::tuple<Geometry, std::shared_ptr<Excitation>>(geometry, excitation);
}

void check(Geometry const &geometry, std::shared_ptr<Excitation const> excitation) {
  scalapack::Context context = scalapack::Context::Squarest();
  scalapack::Sizes const blocks = {16, 16};
  auto const nHarmonics = geometry.objects.size() ? geometry.objects.front().nMax : 1;
  auto const N = geometry.objects.size() * 2 * (HarmonicsIterator::max_flat(nHarmonics) - 1);
  // Compute matrix in parallel
  auto const parallel_local =
      preconditioned_scattering_matrix(geometry, excitation, context, blocks);
  // Compute serial matrix
  auto const serial_context = context.serial();
  auto const N1 = std::max(N, t_uint(1));
  scalapack::Matrix<t_complex> serial(serial_context, {N, N}, {N1, N1});
  if(serial_context.is_valid()) {
    serial.local() = preconditioned_scattering_matrix(geometry, excitation);
    REQUIRE(serial.local().rows() == N);
    REQUIRE(serial.local().cols() == N);
  }
  // transfer serial matrix to parallel context
  scalapack::Matrix<t_complex> parallel(context, {N, N}, blocks);
  serial.transfer_to(parallel.context(), parallel);
  // Check the local matrices are identical
  REQUIRE(parallel.local().rows() == parallel_local.rows());
  REQUIRE(parallel.local().cols() == parallel_local.cols());
  CHECK(parallel.local().isApprox(parallel_local));
  if(context.size() == 1)
    CHECK(parallel.local().isApprox(serial.local()));
}

TEST_CASE("Scattering matrix without remainder") {
  mpi::Communicator world;
  auto const nHarmonics = 5;
  auto const scatterer = default_scatterer(nHarmonics);
  auto const length = default_length();
  auto const input = fcc_system(std::make_tuple(world.size(), 2, 1), length, scatterer);
  check(std::get<0>(input), std::get<1>(input));
}

TEST_CASE("Scattering matrix with remainder") {
  mpi::Communicator world;
  auto const nHarmonics = 3;
  auto const scatterer = default_scatterer(nHarmonics);
  auto const length = default_length();
  auto const cell = fcc_cell();
  SECTION("Remainder > number of local objects") {
    auto input = fcc_system(std::make_tuple(world.size(), 1, 1), length, scatterer);
    auto &geometry = std::get<0>(input);
    for(int i(1); i < std::min(3, static_cast<int>(world.size())); ++i) {
      Eigen::Matrix<t_real, 3, 1> pos =
          cell * Eigen::Matrix<t_real, 3, 1>(-1 - i, -1 - i, -1 - i) * length;
      geometry.pushObject({Tools::toSpherical({pos(0), pos(1), pos(2)}), scatterer.elmag,
                           scatterer.radius, scatterer.nMax});
      geometry.update(std::get<1>(input));
      check(std::get<0>(input), std::get<1>(input));
    }
  }

  SECTION("Remainder < number of local objects") {
    auto input = fcc_system(std::make_tuple(world.size() + 1, world.size(), 1), length, scatterer);
    auto &geometry = std::get<0>(input);
    for(int i(1); i < std::min(3, static_cast<int>(world.size())); ++i) {
      Eigen::Matrix<t_real, 3, 1> pos =
          cell * Eigen::Matrix<t_real, 3, 1>(-1 - i, -1 - i, -1 - i) * length;
      geometry.pushObject({Tools::toSpherical({pos(0), pos(1), pos(2)}), scatterer.elmag,
                           scatterer.radius, scatterer.nMax});
      geometry.update(std::get<1>(input));
      check(std::get<0>(input), std::get<1>(input));
    }
  }
}

TEST_CASE("Zero objects") {
  auto const nHarmonics = 3;
  auto const scatterer = default_scatterer(nHarmonics);
  auto const length = default_length();
  auto const cell = fcc_cell();
  auto input = fcc_system(std::make_tuple(0, 1, 1), length, scatterer);
  check(std::get<0>(input), std::get<1>(input));
}

TEST_CASE("Very small number of objects") {
  mpi::Communicator world;
  auto const nHarmonics = 3;
  auto const scatterer = default_scatterer(nHarmonics);
  auto const length = default_length();
  auto const cell = fcc_cell();
  auto input = fcc_system(std::make_tuple(world.size() - 1, 1, 1), length, scatterer);
  check(std::get<0>(input), std::get<1>(input));
}

TEST_CASE("Distributed source vector") {
  mpi::Communicator world;
  auto const nHarmonics = 5;
  auto const scatterer = default_scatterer(nHarmonics);
  auto const length = default_length();
  auto const input = fcc_system(std::make_tuple(world.size(), 2, 1), length, scatterer);

  scalapack::Context context = scalapack::Context::Squarest();
  scalapack::Sizes const blocks = {16, 16};

  // Compute matrix in parallel
  auto const parallel_matrix =
      preconditioned_scattering_matrix(std::get<0>(input), std::get<1>(input), context, blocks);
  auto const serial_matrix =
      preconditioned_scattering_matrix(std::get<0>(input), std::get<1>(input));

  Vector<t_complex> const serial_vector = Vector<t_complex>::Random(serial_matrix.rows());
  auto const parallel_vector = distributed_source_vector(serial_vector, context, blocks);

  t_uint const N(serial_matrix.rows());
  scalapack::Matrix<t_complex> sca_matrix(context, {N, N}, blocks);
  CHECK(sca_matrix.local().rows() == parallel_matrix.rows());
  CHECK(sca_matrix.local().cols() == parallel_matrix.cols());
  sca_matrix.local() = parallel_matrix;
  scalapack::Matrix<t_complex> sca_vector(context, {N, 1}, blocks);
  CHECK(sca_vector.local().size() == parallel_vector.size());
  sca_vector.local() = parallel_vector;

  scalapack::Matrix<t_complex> sca_result(context, {N, 1}, blocks);
  pdgemm(1e0, sca_matrix, sca_vector, 0, sca_result);
  if(context.is_valid()) {
    CHECK(sca_result.local().rows() ==
          scalapack::Matrix<t_complex>::local_rows(context, {N, 1}, blocks, {0, 0}));
    CHECK(sca_result.local().cols() ==
          scalapack::Matrix<t_complex>::local_cols(context, {N, 1}, blocks, {0, 0}));
  }
  auto const result = gather_all_source_vector(sca_result);
  if(context.is_valid()) {
    REQUIRE(result.rows() == N);
    REQUIRE(result.cols() == 1);
    REQUIRE(result.isApprox(serial_matrix * serial_vector));
  } else {
    REQUIRE(result.rows() == 0);
  }
}
