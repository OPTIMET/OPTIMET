#include "Aliases.h"
#include "Geometry.h"
#include "Scatterer.h"
#include "Solver.h"
#include "Tools.h"
#include "Types.h"
#include "cmdl.h"
#include "constants.h"
#include "mpi/Collectives.hpp"
#include "mpi/Communicator.h"
#include "mpi/Session.h"
#include <benchmark/benchmark.h>
#include <chrono>
#include <iostream>

namespace {
#ifdef OPTIMET_BELOS
//! Holds solver parameters
Teuchos::RCP<Teuchos::ParameterList> parameters;
template <class T> T get_param(std::string const &name, T const &default_) {
  return parameters->get<T>(name, default_);
}
#else
template <class T> T get_param(std::string const &, T const &default_) { return default_; }
#endif

using namespace optimet;
//! Gets an fcc cell
Matrix<t_real> fcc_cell() {
  Matrix<t_real> result = Matrix<t_real>::Ones(3, 3) * 0.5;
  result.diagonal().fill(0);
  return result;
}

t_real default_wavelength() { return 750e-9; }

Run fcc_system(t_int const &N, t_real length, Scatterer const &scatterer) {
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
#ifdef OPTIMET_BELOS
  result.belos_params = parameters;
#endif
  return result;
}

Scatterer default_scatterer(t_int nHarmonics) {
  auto const radius = get_param<t_real>("radius", 1e0) * default_wavelength();
  ElectroMagnetic const elmag{get_param<t_complex>("epsilon_r", 13.1), 1.0};
  return {{0, 0, 0}, elmag, radius, nHarmonics};
}

constexpr t_real default_length() { return 2000e-9; }

#ifndef OPTIMET_MPI
void problem_setup(benchmark::State &state) {
  auto const nHarmonics = state.range_y();
  auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
  auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));

  while(state.KeepRunning())
    optimet::solver::factory(input);
  state.SetItemsProcessed(int64_t(state.iterations()) *
                          int64_t(std::get<0>(input).scatterer_size()));
}

void benchmark_solver(benchmark::State &state) {
  auto const nHarmonics = state.range_y();
  auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
  auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
  auto const solver(input);
  Result result(std::get<0>(input), std::get<1>(input));

  while(state.KeepRunning()) {
    result.internal_coef.fill(0);
    solver.solve(result.scatter_coef, result.internal_coef);
  }
  state.SetItemsProcessed(int64_t(state.iterations()) *
                          int64_t(std::get<0>(input).scatterer_size()));
}
#endif

#ifdef OPTIMET_MPI
void benchmark_solver(benchmark::State &state) {
  mpi::Communicator world;
  auto const nHarmonics = state.range_y();
  auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
  auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
  auto const solver = optimet::solver::factory(input);

  Result result(input.geometry, input.excitation);
  while(state.KeepRunning()) {
    result.internal_coef.fill(0);
    auto start = std::chrono::high_resolution_clock::now();
    solver->solve(result.scatter_coef, result.internal_coef, world);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    auto proc_max = world.all_reduce(elapsed_seconds.count(), MPI_MAX);
    state.SetIterationTime(proc_max);
  }
  state.SetItemsProcessed(int64_t(state.iterations()) * int64_t(solver->scattering_size()));
}
#endif
}

static void CustomArguments(benchmark::internal::Benchmark *b) {
  for(auto const i : {1, 2, 3, 4, 5, 10, 20, 30, 40, 50})
    for(auto const j : {1, 2, 3, 4, 6, 8, 10})
      b->ArgPair(i, j);
}

extern std::string FLAG_benchmark_format;

#ifndef OPTIMET_MPI
BENCHMARK(problem_setup)->Apply(CustomArguments)->Unit(benchmark::kMicrosecond);
BENCHMARK(benchmark_solver)->Apply(CustomArguments)->Unit(benchmark::kMicrosecond);
#else
BENCHMARK(benchmark_solver)->Apply(CustomArguments)->Unit(benchmark::kMicrosecond)->UseManualTime();
#endif

namespace benchmark {
class MPIReporter : public BenchmarkReporter {
public:
  MPIReporter(std::unique_ptr<BenchmarkReporter> reporter, bool doreport)
      : reporter(std::move(reporter)), doreport(doreport) {}
  MPIReporter(std::unique_ptr<BenchmarkReporter> reporter)
      : MPIReporter(std::move(reporter), optimet::mpi::Communicator().rank() == 0) {}
  virtual bool ReportContext(const Context &context) {
    if(doreport)
      return reporter->ReportContext(context);
    else
      return true;
  }
  virtual void ReportRuns(const std::vector<Run> &report) {
    if(doreport)
      reporter->ReportRuns(report);
  }
  virtual void Finalize() {
    if(doreport)
      reporter->Finalize();
  }

private:
  std::unique_ptr<BenchmarkReporter> reporter;
  bool doreport;
};

std::unique_ptr<BenchmarkReporter> parse_cmdl(int argc, char **argv) {
  for(int i(0); i < argc; ++i)
    if(std::string(argv[i]) == "--benchmark_format=console")
      return std::unique_ptr<BenchmarkReporter>{new ConsoleReporter()};
    else if(std::string(argv[i]) == "--benchmark_format=json")
      return std::unique_ptr<BenchmarkReporter>{new JSONReporter()};
    else if(std::string(argv[i]) == "--benchmark_format=csv")
      return std::unique_ptr<BenchmarkReporter>{new CSVReporter()};
  return std::unique_ptr<BenchmarkReporter>{new ConsoleReporter()};
}
}

int main(int argc, char **argv) {
  optimet::mpi::init(argc, const_cast<const char **>(argv));

  ::benchmark::MPIReporter reporter(::benchmark::parse_cmdl(argc, argv));
  ::benchmark::Initialize(&argc, argv);

#ifdef OPTIMET_BELOS
  parameters = parse_cmdl(argc, argv);
#endif

  ::benchmark::RunSpecifiedBenchmarks(&reporter);
  optimet::mpi::finalize();
  return 0;
}
