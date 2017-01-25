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

#include "Geometry.h"
#include "Run.h"
#include "Scatterer.h"
#include "Tools.h"
#include "Types.h"
#include "constants.h"
#include "fcc.h"
#include <benchmark/benchmark.h>

namespace optimet {
namespace {

//! Gets an fcc cell
Matrix<t_real> fcc_cell() {
  Matrix<t_real> result = Matrix<t_real>::Ones(3, 3) * 0.5;
  result.diagonal().fill(0);
  return result;
}

constexpr t_real default_wavelength() { return 750e-9; }
constexpr t_real default_length() { return 2000e-9; }

Scatterer default_scatterer(Teuchos::RCP<Teuchos::ParameterList> const &parameters) {
  auto const nMax = parameters->get<int>("nMax");
  auto const radius = parameters->get<optimet::t_real>("radius", 1.0) * default_length();
  ElectroMagnetic const elmag{parameters->get<optimet::t_complex>("epsilon_r", 13.1), 1.0};
  return {{0, 0, 0}, elmag, radius, nMax};
}
}

Run fcc_input(Teuchos::RCP<Teuchos::ParameterList> const &parameters) {
  auto const N = parameters->get<int>("nObjects");
  auto const length = (parameters->get<optimet::t_real>("radius", 0.5) + 0.5) * default_length();
  auto const scatterer = default_scatterer(parameters);
  // setup geometry
  auto const cell = fcc_cell();
  t_int n = std::pow(N, 1e0 / 3e0);
  auto range = std::make_tuple(n, n, n);
  if(n * n * n < N)
    ++std::get<0>(range);
  if((n + 1) * n * n < N)
    ++std::get<1>(range);
  auto geometry = std::make_shared<Geometry>();
  n = 0;
  for(int i(0); i < std::get<0>(range); ++i)
    for(int j(0); j < std::get<1>(range); ++j)
      for(int k(0); k < std::get<2>(range); ++k, ++n) {
        if(n == N)
          break;
        Eigen::Matrix<t_real, 3, 1> pos = cell * Eigen::Matrix<t_real, 3, 1>(i, j, k) * length;
        geometry->pushObject({Tools::toSpherical({pos(0), pos(1), pos(2)}), scatterer.elmag,
                              scatterer.radius, scatterer.nMax});
      }

  // Create excitation
  Spherical<t_real> const vKinc{2 * consPi / default_wavelength(), 90 * consPi / 180.0,
                                90 * consPi / 180.0};
  SphericalP<t_complex> const Eaux{0e0, 1e0, 0e0};
  auto const excitation =
      std::make_shared<Excitation>(0, Tools::toProjection(vKinc, Eaux), vKinc, scatterer.nMax);
  excitation->populate();
  geometry->update(excitation);
  Run result;
  result.geometry = geometry;
  result.excitation = excitation;

  result.belos_params = parameters;
  result.do_fmm = parameters->get<bool>("do_fmm");
  result.fmm_subdiagonals = parameters->get<t_int>("fmm_subdiagonals");

  return result;
}

void CustomArguments(std::vector<t_uint> const &nharmonics, std::vector<t_uint> const &nparticles,
                     benchmark::internal::Benchmark *b) {
  for(auto const i : nparticles)
    for(auto const j : nharmonics)
      b->ArgPair(i, j);
}

std::unique_ptr<::benchmark::BenchmarkReporter> parse_benchmark_cmdl(int argc, char **argv) {
  for(int i(0); i < argc; ++i)
    if(std::string(argv[i]) == "--benchmark_format=console")
      return std::unique_ptr<::benchmark::BenchmarkReporter>{new ::benchmark::ConsoleReporter()};
    else if(std::string(argv[i]) == "--benchmark_format=json")
      return std::unique_ptr<::benchmark::BenchmarkReporter>{new ::benchmark::JSONReporter()};
    else if(std::string(argv[i]) == "--benchmark_format=csv")
      return std::unique_ptr<::benchmark::BenchmarkReporter>{new ::benchmark::CSVReporter()};
  return std::unique_ptr<::benchmark::BenchmarkReporter>{new ::benchmark::ConsoleReporter()};
}
}
