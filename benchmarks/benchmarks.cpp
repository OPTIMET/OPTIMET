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
#include "PreconditionedMatrix.h"
#include "Scatterer.h"
#include "Solver.h"
#include "Tools.h"
#include "Types.h"
#include "cmdl.h"
#include "constants.h"
#include "fcc.h"
#include "mpi/Collectives.hpp"
#include "mpi/Communicator.h"
#include "mpi/FastMatrixMultiply.h"
#include "mpi/Session.h"
#include <benchmark/benchmark.h>
#include <chrono>
#include <iostream>

namespace optimet {
namespace {

std::shared_ptr<std::vector<t_uint> const> nharmonics;
std::shared_ptr<std::vector<t_uint> const> nparticles;

void NeedExtraDataArgument(benchmark::internal::Benchmark *b) {
  return CustomArguments(*nharmonics, *nparticles, b);
}

#ifndef OPTIMET_MPI
OPTIMET_BENCHMARK(serial_problem_setup) {

  OPTIMET_BENCHMARK_TIME_START;
  source_vector(*input.geometry, input.excitation);
  preconditioned_scattering_matrix(*input.geometry, input.excitation);
  OPTIMET_BENCHMARK_TIME_END;
}
#endif

#ifdef OPTIMET_SCALAPACK
OPTIMET_BENCHMARK(scalapack_problem_setup) {
  scalapack::Sizes const block_size = {64, 64};
  auto const context = scalapack::Context::Squarest();
  if(static_cast<unsigned long long>(size * size) * 16ull >= 2ull * 1024ull * 1024ull * 1024ull)
    state.SkipWithError("Matrix too large for scalapack");

  OPTIMET_BENCHMARK_TIME_START;
  distributed_source_vector(source_vector(*input.geometry, input.excitation), context, block_size);
  preconditioned_scattering_matrix(*input.geometry, input.excitation, context, block_size);
  OPTIMET_BENCHMARK_TIME_END;
}
#endif

OPTIMET_BENCHMARK(fmm_problem_setup) {
  auto const wavenumber = input.excitation->wavenumber();

  OPTIMET_BENCHMARK_TIME_START;
#ifdef OPTIMET_MPI
  mpi::FastMatrixMultiply(input.geometry->bground, wavenumber, input.geometry->objects,
                          input.fmm_subdiagonals, input.communicator);
#else
  optimet::FastMatrixMultiply(input.geometry->bground, wavenumber, input.geometry->objects);
#endif
  OPTIMET_BENCHMARK_TIME_END;
}

#ifndef OPTIMET_MPI
OPTIMET_BENCHMARK(serial_solver) {
  input.belos_params->set("Solver", "eigen");
  Result result(input.geometry, input.excitation);
  auto const solver = solver::factory(input);

  OPTIMET_BENCHMARK_TIME_START;
  result.internal_coef.fill(0);
  solver->solve(result.scatter_coef, result.internal_coef);
  OPTIMET_BENCHMARK_TIME_END;
}
#endif

#ifdef OPTIMET_SCALAPACK
OPTIMET_BENCHMARK(scalapack_solver) {
  if(static_cast<unsigned long long>(size * size) * 16ull >= 2ull * 1024ull * 1024ull * 1024ull) {
    state.SkipWithError("Matrix too large for scalapack");
    return;
  }

  input.belos_params->set("Solver", "scalapack");
  input.belos_params->set("fmm", false);
  Result result(input.geometry, input.excitation);
  auto const solver = solver::factory(input);

  OPTIMET_BENCHMARK_TIME_START;
  result.internal_coef.fill(0);
  solver->solve(result.scatter_coef, result.internal_coef);
  OPTIMET_BENCHMARK_TIME_END;
}

OPTIMET_BENCHMARK(gmres_solver) {
  if(static_cast<unsigned long long>(size * size) * 16ull >= 2ull * 1024ull * 1024ull * 1024ull) {
    state.SkipWithError("Matrix too large for scalapack");
    return;
  }

  input.belos_params->set("Solver", "GMRES");
  input.belos_params->set("fmm", false);
  Result result(input.geometry, input.excitation);
  auto const solver = solver::factory(input);

  OPTIMET_BENCHMARK_TIME_START;
  result.internal_coef.fill(0);
  solver->solve(result.scatter_coef, result.internal_coef);
  OPTIMET_BENCHMARK_TIME_END;
}
#endif

OPTIMET_BENCHMARK(fmm_solver) {
  input.belos_params->set("Solver", "GMRES");
  input.belos_params->set("fmm", true);
  Result result(input.geometry, input.excitation);
  auto const solver = solver::factory(input);

  OPTIMET_BENCHMARK_TIME_START;
  result.internal_coef.fill(0);
  solver->solve(result.scatter_coef, result.internal_coef);
  OPTIMET_BENCHMARK_TIME_END;
}

#ifndef OPTIMET_MPI
OPTIMET_BENCHMARK(eigen_multiplication) {
  auto const Q = source_vector(*input.geometry, input.excitation);
  auto const S = preconditioned_scattering_matrix(*input.geometry, input.excitation);
  Vector<t_complex> result(Q.size());

  OPTIMET_BENCHMARK_TIME_START;
  result = Q * S;
  OPTIMET_BENCHMARK_TIME_END;
}
#endif

#ifdef OPTIMET_SCALAPACK
OPTIMET_BENCHMARK(scalapack_multiplication) {
  if(static_cast<unsigned long long>(size * size) * 16ull >= 2ull * 1024ull * 1024ull * 1024ull) {
    state.SkipWithError("Matrix too large for scalapack");
    return;
  }
  auto const N = input.geometry->scatterer_size();
  auto const &context = input.context;
  scalapack::Sizes const block_size = {64, 64};
  auto const Q = distributed_source_vector(source_vector(*input.geometry, input.excitation),
                                           context, block_size);
  auto const S =
      preconditioned_scattering_matrix(*input.geometry, input.excitation, context, block_size);
  scalapack::Matrix<t_complex> Aparallel(context, {N, N}, block_size);
  if(Aparallel.size() > 0)
    Aparallel.local() = S;
  scalapack::Matrix<t_complex> bparallel(context, {N, 1}, block_size);
  if(bparallel.local().size() > 0)
    bparallel.local() = Q;
  scalapack::Matrix<t_complex> result(context, {N, 1}, block_size);
  if(result.local().size() > 0)
    result.local() = Q;

  OPTIMET_BENCHMARK_TIME_START;
  pdgemm(1e0, Aparallel, bparallel, 0.0, result);
  OPTIMET_BENCHMARK_TIME_END;
}
#endif

OPTIMET_BENCHMARK(fmm_multiplication) {
#ifdef OPTIMET_MPI
  mpi::FastMatrixMultiply const fmm(input.geometry->bground, input.excitation->wavenumber(),
                                    input.geometry->objects, input.fmm_subdiagonals,
                                    input.communicator);
  auto const distribution =
      mpi::details::vector_distribution(input.geometry->objects.size(), input.communicator.size());
  auto const first = std::find(distribution.data(), distribution.data() + distribution.size(),
                               input.communicator.rank()) -
                     distribution.data();
  auto const last =
      std::find_if(distribution.data() + first, distribution.data() + distribution.size(),
                   [&input](t_int value) { return value != input.communicator.rank(); }) -
      distribution.data();
  auto const Q = source_vector(input.geometry->objects.begin() + first,
                               input.geometry->objects.begin() + last, input.excitation);
#else
  optimet::FastMatrixMultiply const fmm(input.geometry->bground, input.excitation->wavenumber(),
                                        input.geometry->objects);
  auto const Q = source_vector(*input.geometry, input.excitation);
#endif
  Vector<t_complex> result(Q.size());

  OPTIMET_BENCHMARK_TIME_START;
  fmm(Q, result);
  OPTIMET_BENCHMARK_TIME_END;
}
}
}

int main(int argc, char **argv) {
  using namespace optimet;
  mpi::init(argc, const_cast<const char **>(argv));

  MPIReporter reporter(parse_benchmark_cmdl(argc, argv));

  auto const parameters = parse_cmdl(argc, argv);
  nharmonics = std::make_shared<const std::vector<t_uint>>(
      convert_string<t_uint>(parameters->get<std::string>("nharmonics")));
  nparticles = std::make_shared<const std::vector<t_uint>>(
      convert_string<t_uint>(parameters->get<std::string>("nparticles")));

#ifndef OPTIMET_MPI
  OPTIMET_REGISTER_BENCHMARK(serial_problem_setup)->Unit(benchmark::kMicrosecond);
#elif defined(OPTIMET_SCALAPACK)
  OPTIMET_REGISTER_BENCHMARK(scalapack_problem_setup)->Unit(benchmark::kMicrosecond);
#endif
  OPTIMET_REGISTER_BENCHMARK(fmm_problem_setup)->Unit(benchmark::kMicrosecond);

#ifndef OPTIMET_MPI
    OPTIMET_REGISTER_BENCHMARK(serial_solver)->Unit(benchmark::kMicrosecond);
#else
#ifdef OPTIMET_SCALAPACK
    OPTIMET_REGISTER_BENCHMARK(scalapack_solver)->Unit(benchmark::kMicrosecond);
    OPTIMET_REGISTER_BENCHMARK(gmres_solver)->Unit(benchmark::kMicrosecond);
#endif
    OPTIMET_REGISTER_BENCHMARK(fmm_solver)->Unit(benchmark::kMicrosecond);
#endif

#ifndef OPTIMET_MPI
  OPTIMET_REGISTER_BENCHMARK(eigen_multiplication)->Unit(benchmark::kMicrosecond);
#elif defined(OPTIMET_SCALAPACK)
  OPTIMET_REGISTER_BENCHMARK(scalapack_multiplication)->Unit(benchmark::kMicrosecond);
#endif
  OPTIMET_REGISTER_BENCHMARK(fmm_multiplication)->Unit(benchmark::kMicrosecond);

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks(&reporter);
  mpi::finalize();
  return 0;
}
