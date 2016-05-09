#include "catch.hpp"
#include <iostream>

#include "HarmonicsIterator.h"
#include "Solver.h"
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
  geometry.update(excitation.get());
  return {geometry, excitation};
}

void check(Geometry const &geometry, Excitation const &excitation) {
  scalapack::Context context;
  scalapack::Sizes const blocks = {16, 16};
  auto const nHarmonics = geometry.objects.front().nMax;
  auto const N = geometry.objects.size() * 2 * (HarmonicsIterator::max_flat(nHarmonics) - 1);
  // Compute matrix in parallel
  auto const parallel_local =
      preconditioned_scattering_matrix(geometry, excitation, context, blocks);
  // Compute serial matrix
  auto const serial_context = context.serial();
  scalapack::Matrix<t_complex> serial(serial_context, {N, N}, {N, N});
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
  CAPTURE(parallel.local());
  CAPTURE(parallel_local);
  CHECK(parallel.local().isApprox(parallel_local));
}

TEST_CASE("Scattering matrix without remainder") {
  mpi::Communicator world;
  auto const nHarmonics = 5;
  auto const scatterer = default_scatterer(nHarmonics);
  auto const length = default_length();
  auto const input = fcc_system({world.size(), 2, 1}, length, scatterer);
  check(std::get<0>(input), *std::get<1>(input));
}

TEST_CASE("Scattering matrix with remainder") {
  mpi::Communicator world;
  auto const nHarmonics = 3;
  auto const scatterer = default_scatterer(nHarmonics);
  auto const length = default_length();
  auto const cell = fcc_cell();
  SECTION("Remainder > number of local objects") {
    auto input = fcc_system({world.size(), 1, 1}, length, scatterer);
    auto &geometry = std::get<0>(input);
    for(int i(1); i < std::min(3, static_cast<int>(world.size())); ++i) {
      Eigen::Matrix<t_real, 3, 1> pos =
          cell * Eigen::Matrix<t_real, 3, 1>(-1 - i, -1 - i, -1 - i) * length;
      geometry.pushObject({Tools::toSpherical({pos(0), pos(1), pos(2)}), scatterer.elmag,
                           scatterer.radius, scatterer.nMax});
      geometry.update(std::get<1>(input).get());
      check(std::get<0>(input), *std::get<1>(input));
    }
  }

  SECTION("Remainder < number of local objects") {
    auto input = fcc_system({world.size() + 1, world.size(), 1}, length, scatterer);
    auto &geometry = std::get<0>(input);
    for(int i(1); i < std::min(3, static_cast<int>(world.size())); ++i) {
      Eigen::Matrix<t_real, 3, 1> pos =
          cell * Eigen::Matrix<t_real, 3, 1>(-1 - i, -1 - i, -1 - i) * length;
      geometry.pushObject({Tools::toSpherical({pos(0), pos(1), pos(2)}), scatterer.elmag,
                           scatterer.radius, scatterer.nMax});
      geometry.update(std::get<1>(input).get());
      check(std::get<0>(input), *std::get<1>(input));
    }
  }
}
