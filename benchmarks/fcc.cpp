#include "Aliases.h"
#include "Geometry.h"
#include "Scatterer.h"
#include "Solver.h"
#include "Tools.h"
#include "Types.h"
#include "constants.h"
#include <benchmark/benchmark.h>

namespace {
using namespace optimet;
//! Gets an fcc cell
Matrix<t_real> fcc_cell() {
  Matrix<t_real> result = Matrix<t_real>::Ones(3, 3) * 0.5;
  result.diagonal().fill(0);
  return result;
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
  return {geometry, excitation};
}

Scatterer default_scatterer(t_int nHarmonics) {
  auto const radius = 500e-9;
  ElectroMagnetic const elmag{13.1, 1.0};
  return {{0, 0, 0}, elmag, radius, nHarmonics};
}

constexpr t_real default_length() { return 2000e-9; }

void problem_setup(benchmark::State &state) {
  auto const range = std::make_tuple(state.range_x(), state.range_x(), state.range_x());
  auto const nHarmonics = state.range_y();
  ElectroMagnetic const elmag{13.1, 1.0};
  auto input = fcc_system(range, default_length(), default_scatterer(nHarmonics));

  while(state.KeepRunning())
    Solver(&std::get<0>(input), std::get<1>(input), O3DSolverIndirect, nHarmonics);
}

void serial(benchmark::State &state) {
  auto const range = std::make_tuple(state.range_x(), state.range_x(), state.range_x());
  auto const nHarmonics = state.range_y();
  ElectroMagnetic const elmag{13.1, 1.0};
  auto input = fcc_system(range, default_length(), default_scatterer(nHarmonics));
  Solver solver(&std::get<0>(input), std::get<1>(input), O3DSolverIndirect, nHarmonics);
  Result result(&std::get<0>(input), std::get<1>(input), nHarmonics);

  while(state.KeepRunning()) {
    result.internal_coef.fill(0);
    solver.solve(result.scatter_coef, result.internal_coef);
  }
  state.SetBytesProcessed(int64_t(state.iterations()) * int64_t(result.internal_coef.size()));
}
}

std::tuple<int, int> const nGeom = {1, 4};
std::tuple<int, int> const nHarm = {1, 20};
BENCHMARK(problem_setup)
    ->RangePair(std::get<0>(nGeom), std::get<1>(nGeom), std::get<0>(nHarm), std::get<1>(nHarm));
BENCHMARK(serial)
    ->RangePair(std::get<0>(nGeom), std::get<1>(nGeom), std::get<0>(nHarm), std::get<1>(nHarm));
BENCHMARK_MAIN();
