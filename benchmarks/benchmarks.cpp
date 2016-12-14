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

OPTIMET_BENCHMARK(serial_problem_setup) {

  OPTIMET_BENCHMARK_TIME_START;
  source_vector(*input.geometry, input.excitation);
  preconditioned_scattering_matrix(*input.geometry, input.excitation);
  OPTIMET_BENCHMARK_TIME_END;
}

#ifdef OPTIMET_SCALAPACK
OPTIMET_BENCHMARK(scalapack_problem_setup) {
  scalapack::Sizes const block_size = {64, 64};
  auto const context = scalapack::Context::Squarest();

  OPTIMET_BENCHMARK_TIME_START;
  distributed_source_vector(source_vector(*input.geometry, input.excitation), context, block_size);
  preconditioned_scattering_matrix(*input.geometry, input.excitation, context, block_size);
  OPTIMET_BENCHMARK_TIME_END;
}
#endif

#ifdef OPTIMET_MPI
OPTIMET_BENCHMARK(fmm_problem_setup) {
  auto const wavenumber = input.excitation->wavenumber();

  OPTIMET_BENCHMARK_TIME_START;
  mpi::FastMatrixMultiply(input.geometry->bground, wavenumber, input.geometry->objects,
                          input.fmm_subdiagonals, world);
  OPTIMET_BENCHMARK_TIME_END;
}
#endif

OPTIMET_BENCHMARK(serial_solver) {
  input.belos_params->set("Solver", "eigen");
  Result result(input.geometry, input.excitation);
  auto const solver = solver::factory(input);

  OPTIMET_BENCHMARK_TIME_START;
  result.internal_coef.fill(0);
  solver->solve(result.scatter_coef, result.internal_coef);
  OPTIMET_BENCHMARK_TIME_END;
}

OPTIMET_BENCHMARK(eigen_multiplication) {
  auto const Q = source_vector(*input.geometry, input.excitation);
  auto const S = preconditioned_scattering_matrix(*input.geometry, input.excitation);
  Vector<t_complex> result(Q.size());

  OPTIMET_BENCHMARK_TIME_START;
  result = Q * S;
  OPTIMET_BENCHMARK_TIME_END;
}

#ifdef OPTIMET_SCALAPACK
OPTIMET_BENCHMARK(scalapack_solver) {
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

#ifdef OPTIMET_MPI
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
#endif

//
// void benchmark_solver(benchmark::State &state) {
//   mpi::Communicator world;
//   auto const nHarmonics = state.range_y();
//   auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
//   auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
//   int64_t const size = input.geometry->scatterer_size();
//   auto const solver = solver::factory(input);
//
//   Result result(input.geometry, input.excitation);
//   while(state.KeepRunning()) {
//     result.internal_coef.fill(0);
//     auto start = std::chrono::high_resolution_clock::now();
//     solver->solve(result.scatter_coef, result.internal_coef, world);
//     auto end = std::chrono::high_resolution_clock::now();
//     auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end -
//     start);
//     auto proc_max = world.all_reduce(elapsed_seconds.count(), MPI_MAX);
//     state.SetIterationTime(proc_max);
//   }
//   state.SetItemsProcessed(int64_t(state.iterations()) * size);
// }
//
// void serial_matrix_multiplication(benchmark::State &state) {
//   mpi::Communicator world;
//   auto const nHarmonics = state.range_y();
//   auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
//   auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
//   int64_t const size = input.geometry->scatterer_size();
//   auto const Q = source_vector(*input.geometry, input.excitation);
//   auto const S = preconditioned_scattering_matrix(*input.geometry, input.excitation);
//   Vector<t_complex> result(Q.size());
//
//   while(state.KeepRunning()) {
//     auto start = std::chrono::high_resolution_clock::now();
//     result = Q * S;
//     auto end = std::chrono::high_resolution_clock::now();
//     auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end -
//     start);
//     auto proc_max = world.all_reduce(elapsed_seconds.count(), MPI_MAX);
//     state.SetIterationTime(proc_max);
//   }
//   state.SetItemsProcessed(int64_t(state.iterations()) * size);
// }
//
// #ifdef OPTIMET_SCALAPACK
// void scalapack_matrix_multiplication(benchmark::State &state) {
//   mpi::Communicator world;
//   scalapack::Context context = scalapack::Context::Squarest();
//   auto const nHarmonics = state.range_y();
//   auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
//   auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
//   int64_t const size = input.geometry->scatterer_size();
//   scalapack::Sizes const block_size = {64, 64};
//   auto const Q = distributed_source_vector(source_vector(*input.geometry, input.excitation),
//                                            context, block_size);
//   auto const S = preconditioned_scattering_matrix(*input.geometry, input.excitation);
//   auto const N = input.geometry->scatterer_size();
//   scalapack::Matrix<t_complex> Aparallel(context, {N, N}, block_size);
//   if(Aparallel.size() > 0)
//     Aparallel.local() = S;
//   scalapack::Matrix<t_complex> bparallel(context, {N, 1}, block_size);
//   if(bparallel.local().size() > 0)
//     bparallel.local() = Q;
//   scalapack::Matrix<t_complex> result(context, {N, 1}, block_size);
//   if(result.local().size() > 0)
//     result.local() = Q;
//
//   while(state.KeepRunning()) {
//     auto start = std::chrono::high_resolution_clock::now();
//     pdgemm(1.5e0, Aparallel, bparallel, 0.0, result);
//     auto end = std::chrono::high_resolution_clock::now();
//     auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end -
//     start);
//     auto proc_max = world.all_reduce(elapsed_seconds.count(), MPI_MAX);
//     state.SetIterationTime(proc_max);
//   }
//   state.SetItemsProcessed(int64_t(state.iterations()) * size);
// }
// #else
// void scalapack_matrix_multiplication(benchmark::State &state) {
//   throw std::runtime_error("Optimet not compiled with scalapack");
// }
// #endif
//
// void fmm_multiplication(benchmark::State &state) {
//   mpi::Communicator world;
//   auto const nHarmonics = state.range_y();
//   auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
//   auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
//   int64_t const size = input.geometry->scatterer_size();
//   mpi::FastMatrixMultiply const fmm(input.geometry->bground, input.excitation->wavenumber(),
//                                     input.geometry->objects, input.fmm_subdiagonals, world);
//   auto const distribution =
//       mpi::details::vector_distribution(input.geometry->objects.size(), world.size());
//   auto const first =
//       std::find(distribution.data(), distribution.data() + distribution.size(), world.rank()) -
//       distribution.data();
//   auto const last =
//       std::find_if(distribution.data() + first, distribution.data() + distribution.size(),
//                    [&world](t_int value) { return value != world.rank(); }) -
//       distribution.data();
//   auto const Q = source_vector(input.geometry->objects.begin() + first,
//                                input.geometry->objects.begin() + last, input.excitation);
//   Vector<t_complex> result(Q.size());
//
//   while(state.KeepRunning()) {
//     auto start = std::chrono::high_resolution_clock::now();
//     fmm(Q, result);
//     auto end = std::chrono::high_resolution_clock::now();
//     auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end -
//     start);
//     auto proc_max = world.all_reduce(elapsed_seconds.count(), MPI_MAX);
//     state.SetIterationTime(proc_max);
//   }
//   state.SetItemsProcessed(int64_t(state.iterations()) * size);
// }
//
// void matrix_multiplication(benchmark::State &state) {
//   if(parameters->get<std::string>("Solver") == "eigen")
//     return serial_matrix_multiplication(state);
//   if(parameters->get<std::string>("Solver") != "eigen" and not parameters->get<bool>("do_fmm"))
//     return scalapack_matrix_multiplication(state);
//   if(parameters->get<std::string>("Solver") != "eigen" and
//      parameters->get<std::string>("Solver") != "scalapack" and parameters->get<bool>("do_fmm"))
//     return fmm_multiplication(state);
//   throw std::runtime_error("Invalid benchmark request");
// }
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

  OPTIMET_REGISTER_BENCHMARK(serial_problem_setup)->Unit(benchmark::kMicrosecond);
#ifdef OPTIMET_SCALAPACK
  OPTIMET_REGISTER_BENCHMARK(fmm_problem_setup)->Unit(benchmark::kMicrosecond);
#endif
#ifdef OPTIMET_MPI
  OPTIMET_REGISTER_BENCHMARK(scalapack_problem_setup)->Unit(benchmark::kMicrosecond);
#endif

  OPTIMET_REGISTER_BENCHMARK(serial_solver)->Unit(benchmark::kMicrosecond);
#ifdef OPTIMET_SCALAPACK
  OPTIMET_REGISTER_BENCHMARK(scalapack_solver)->Unit(benchmark::kMicrosecond);
  OPTIMET_REGISTER_BENCHMARK(gmres_solver)->Unit(benchmark::kMicrosecond);
#endif
#ifdef OPTIMET_MPI
  OPTIMET_REGISTER_BENCHMARK(fmm_solver)->Unit(benchmark::kMicrosecond);
#endif

  OPTIMET_REGISTER_BENCHMARK(eigen_multiplication)->Unit(benchmark::kMicrosecond);

  // #ifdef OPTIMET_MPI
  //   BENCHMARK(matrix_multiplication)
  //       ->Apply(CustomArguments)
  //       ->Unit(benchmark::kMillisecond)
  //       ->UseManualTime();
  //   BENCHMARK(problem_setup)->Apply(CustomArguments)->Unit(benchmark::kMillisecond)->UseManualTime();
  //   BENCHMARK(benchmark_solver)
  //       ->Apply(CustomArguments)
  //       ->Unit(benchmark::kMillisecond)
  //       ->UseManualTime();
  // #endif

  ::benchmark::Initialize(&argc, argv);
  ::benchmark::RunSpecifiedBenchmarks(&reporter);
  mpi::finalize();
  return 0;
}
