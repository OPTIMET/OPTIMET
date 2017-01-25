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
#include "PreconditionedMatrix.h"
#include "Tools.h"
#include "Types.h"
#include "mpi/Communicator.h"
#include "mpi/Session.h"
#include "scalapack/Context.h"
#include "scalapack/InitExit.h"
#include "scalapack/Matrix.h"
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

int main(int argc, char *const argv[]) {
  using namespace optimet;
  mpi::init(argc, const_cast<const char **>(argv));
  mpi::Communicator const world;

  auto const iterations = find_arg<t_int>(argc, argv, "iterations", 50);
  auto const warmup = find_arg<t_int>(argc, argv, "warmup", 3);
  auto const nobjects = find_arg<t_int>(argc, argv, "nobjects", 100);
  auto const radius = find_arg<t_real>(argc, argv, "radius", 0.25);
  auto const nMax = find_arg<t_int>(argc, argv, "nharmonics", 10);
  auto const size = 2 * nobjects * nMax * (nMax + 2);
  ElectroMagnetic const elmag{13.1, 1.0};
  auto const length = (radius + 0.5) * default_length();
  Scatterer const scatterer = {{0, 0, 0}, elmag, radius * default_length(), nMax};

  // setup geometry
  auto const cell = fcc_cell();
  t_int n = std::pow(nobjects, 1e0 / 3e0);
  auto range = std::make_tuple(n, n, n);
  if(n * n * n < nobjects)
    ++std::get<0>(range);
  if((n + 1) * n * n < nobjects)
    ++std::get<1>(range);
  auto geometry = std::make_shared<Geometry>();
  n = 0;
  for(int i(0); i < std::get<0>(range); ++i)
    for(int j(0); j < std::get<1>(range); ++j)
      for(int k(0); k < std::get<2>(range); ++k, ++n) {
        if(n == nobjects)
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

  if(static_cast<unsigned long long>(size * size) * 16ull >= 2ull * 1024ull * 1024ull * 1024ull) {
    std::cerr << "Matrix too large for scalapack\n";
    scalapack::finalize(1);
    mpi::finalize();
    return 1;
  }
  auto const N = geometry->scatterer_size();
  auto const context = scalapack::Context::Squarest();
  scalapack::Sizes const block_size = {64, 64};
  auto const Q =
      distributed_source_vector(source_vector(*geometry, excitation), context, block_size);
  auto const S = preconditioned_scattering_matrix(*geometry, excitation, context, block_size);
  scalapack::Matrix<t_complex> Aparallel(context, {N, N}, block_size);
  if(Aparallel.size() > 0)
    Aparallel.local() = S;
  scalapack::Matrix<t_complex> bparallel(context, {N, 1}, block_size);
  if(bparallel.local().size() > 0)
    bparallel.local() = Q;
  scalapack::Matrix<t_complex> result(context, {N, 1}, block_size);
  if(result.local().size() > 0)
    result.local() = Q;

  for(int i(0); i < warmup; ++i)
    pdgemm(1e0, Aparallel, bparallel, 0.0, result);
  t_real elapsed(0);
  for(int i(0); i < iterations; ++i) {
    auto start = std::chrono::high_resolution_clock::now();
    pdgemm(1e0, Aparallel, bparallel, 0.0, result);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);

    elapsed += world.all_reduce(elapsed_seconds.count(), MPI_MAX);
  }

  if(world.is_root()) {
    std::cout << "scalapack multiplication:\n";
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
    std::cout << "    nobjects: " << nobjects << "\n";
    std::cout << "    iterations: " << iterations << "\n";
    std::cout << "    Total time: " << elapsed << " seconds\n";
    std::cout << "    Timing: " << elapsed / iterations << " seconds\n";
    std::cout << "---\n";
  }

  scalapack::finalize(1);
  mpi::finalize();
  return 0;
}
