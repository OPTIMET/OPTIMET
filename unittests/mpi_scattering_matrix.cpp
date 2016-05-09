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

TEST_CASE("Scattering matrix without remainder") {
  mpi::Communicator world;
  if(world.size() == 1)
    return;
  scalapack::Sizes const blocks = {16, 16};
  auto const nHarmonics = 5;
  auto const input =
      fcc_system({world.size(), 2, 1}, default_length(), default_scatterer(nHarmonics));

  scalapack::Context context;
  auto const parallel_local =
      preconditioned_scattering_matrix(std::get<0>(input), *std::get<1>(input), context, blocks);
  auto const serial_context = context.serial();
  if(serial_context.is_valid())
    REQUIRE(serial_context.size() == 1);
  else
    REQUIRE(serial_context.size() == 0);
  scalapack::Matrix<t_complex> serial(
      serial_context,
      {std::get<0>(input).objects.size() * 2 * (HarmonicsIterator::max_flat(nHarmonics) - 1),
       std::get<0>(input).objects.size() * 2 * (HarmonicsIterator::max_flat(nHarmonics) - 1)},
      {std::get<0>(input).objects.size() * 2 * (HarmonicsIterator::max_flat(nHarmonics) - 1),
       std::get<0>(input).objects.size() * 2 * (HarmonicsIterator::max_flat(nHarmonics) - 1)});
  if(serial_context.is_valid()) {
    serial.local() = preconditioned_scattering_matrix(std::get<0>(input), *std::get<1>(input));
    REQUIRE(serial.local().rows() ==
            std::get<0>(input).objects.size() * 2 * (HarmonicsIterator::max_flat(nHarmonics) - 1));
    REQUIRE(serial.local().cols() ==
            std::get<0>(input).objects.size() * 2 * (HarmonicsIterator::max_flat(nHarmonics) - 1));
  }
  scalapack::Matrix<t_complex> parallel(
      context,
      {std::get<0>(input).objects.size() * 2 * (HarmonicsIterator::max_flat(nHarmonics) - 1),
       std::get<0>(input).objects.size() * 2 * (HarmonicsIterator::max_flat(nHarmonics) - 1)},
      blocks);
  serial.transfer_to(parallel.context(), parallel);
  CHECK(parallel.local().rows() == parallel_local.rows());
  CHECK(parallel.local().cols() == parallel_local.cols());
  CHECK(parallel.local().isApprox(parallel_local));
}
