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

#include "Excitation.h"
#include "Geometry.h"
#include "Result.h"
#include "Run.h"
#include "Solver.h"
#include "Tools.h"
#include "Types.h"
#include "cmdl.h"
#include "mpi/FastMatrixMultiply.h"
#if defined(OPTIMET_MPI) && !defined(OPTIMET_JUST_DO_SERIAL)
#include "mpi/Session.h"
#endif
#include <chrono>
#include <iostream>
#include <sstream>

namespace {
constexpr optimet::t_real default_wavelength() { return 750e-9; }
constexpr optimet::t_real default_length() { return 2000e-9; }
optimet::Matrix<optimet::t_real> fcc_cell() {
  optimet::Matrix<optimet::t_real> result = optimet::Matrix<optimet::t_real>::Ones(3, 3) * 0.5;
  result.diagonal().fill(0);
  return result;
}
}

template <class T>
T find_arg(int argc, char *const argv[], std::string const &arg, T const &default_) {
  for(int i(0); i < argc - 1; ++i)
    if(std::string(argv[i]) == ("--" + arg)) {
      std::istringstream sstr(argv[i + 1]);
      T result;
      sstr >> result;
      return result;
    }
  return default_;
}

int main(int argc, char *argv[]) {
  using namespace optimet;
  mpi::init(argc, const_cast<const char **>(argv));
  mpi::Communicator const world;

  auto const parameters = parse_cmdl(argc, argv);
  auto const iterations = parameters->get<t_int>("iterations", 10);
  auto const warmup = parameters->get<t_int>("warmup", 2);
  auto const nparticles = std::stoi(parameters->get("nparticles", "100"));
  auto const radius = parameters->get<t_real>("radius", 0.25);
  auto const nMax = std::stoi(parameters->get("nharmonics", "10"));
  ElectroMagnetic const elmag{13.1, 1.0};
  auto const length = (radius + 0.5) * default_length();
  Scatterer const scatterer = {{0, 0, 0}, elmag, radius * default_length(), nMax};

  // setup geometry
  auto const cell = fcc_cell();
  t_int n = std::pow(nparticles, 1e0 / 3e0);
  auto range = std::make_tuple(n, n, n);
  if(n * n * n < nparticles)
    ++std::get<0>(range);
  if((n + 1) * n * n < nparticles)
    ++std::get<1>(range);
  auto geometry = std::make_shared<Geometry>();
  n = 0;
  for(int i(0); i < std::get<0>(range); ++i)
    for(int j(0); j < std::get<1>(range); ++j)
      for(int k(0); k < std::get<2>(range); ++k, ++n) {
        if(n == nparticles)
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

  Run input;
  input.geometry = geometry;
  input.excitation = excitation;
  input.fmm_subdiagonals = parameters->get<t_int>(
      "fmm_subdiagonals", std::max<int>(1, geometry->objects.size() / 2 - 2));
  input.belos_params = parameters;
	input.do_fmm = input.belos_params->get<bool>("do fmm", true);
  if(input.do_fmm and input.belos_params->get("Solver", "scalapack") == "scalapack")
    input.belos_params->set("Solver", "GMRES");

  Result result(input.geometry, input.excitation);
  auto const solver = solver::factory(input);

  for(int i(0); i < warmup; ++i) {
    result.scatter_coef.fill(0);
    result.internal_coef.fill(0);
    solver->solve(result.scatter_coef, result.internal_coef);
  }
  t_real elapsed(0);
  for(int i(0); i < iterations; ++i) {
    result.scatter_coef.fill(0);
    result.internal_coef.fill(0);
    auto start = std::chrono::high_resolution_clock::now();
    solver->solve(result.scatter_coef, result.internal_coef);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

    elapsed += world.all_reduce(elapsed_seconds.count(), MPI_MAX);
  }

  if(world.is_root()) {
    std::cout << (input.do_fmm ? "fmm": "scalapack") << " solver:\n";
#ifdef __APPLE__
    std::cout << "    os: Apple\n";
#else
    std::cout << "    os: Unix\n";
#endif
#ifdef __INTEL_COMPILER
    std::cout << "    compiler: intel " << __VERSION__ << "\n";
#elif defined(__APPLE_CC__)
    std::cout << "    compiler: clang " << __VERSION__ << "\n";
#elif defined(__GNUC__)
    std::cout << "    compiler: gnu " << __VERSION__ << "\n";
#else
    std::cout << "    compiler: unknown " << __VERSION__ << "\n";
#endif
    std::cout << "    program: " << argv[0] << "\n";
    std::cout << "    nprocs: " << world.size() << "\n";
    std::cout << "    nharmonics: " << nMax << "\n";
    std::cout << "    nobjects: " << nparticles << "\n";
    std::cout << "    iterations: " << iterations << "\n";
    std::cout << "    tolerance: " << input.belos_params->get<t_real>("Convergence Tolerance")
              << "\n";
    std::cout << "    solver: " << input.belos_params->get<std::string>("Solver") << "\n";
    std::cout << "    Total time: " << elapsed << " seconds\n";
    std::cout << "    Timing: " << elapsed / iterations << " seconds\n";
    std::cout << "---\n";
  }

  scalapack::finalize(1);
  mpi::finalize();
  return 0;
}
