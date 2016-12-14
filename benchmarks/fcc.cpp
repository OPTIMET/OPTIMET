#include "Aliases.h"
#include "Geometry.h"
#include "PreconditionedMatrix.h"
#include "Scatterer.h"
#include "Solver.h"
#include "Tools.h"
#include "Types.h"
#include "cmdl.h"
#include "constants.h"
#include "mpi/Collectives.hpp"
#include "mpi/Communicator.h"
#include "mpi/FastMatrixMultiply.h"
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
  result.do_fmm = parameters->get<bool>("do_fmm");
  result.fmm_subdiagonals = parameters->get<t_int>("fmm_subdiagonals");
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
  int64_t const size = input.geometry->scatterer_size();

  while(state.KeepRunning())
    solver::factory(input);
  state.SetItemsProcessed(int64_t(state.iterations()) * size);
}

void benchmark_solver(benchmark::State &state) {
  auto const nHarmonics = state.range_y();
  auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
  auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
  int64_t const size = input.geometry->scatterer_size();
  auto const solver = solver::factory(input);
  Result result(std::get<0>(input), std::get<1>(input));

  while(state.KeepRunning()) {
    result.internal_coef.fill(0);
    solver.solve(result.scatter_coef, result.internal_coef);
  }
  state.SetItemsProcessed(int64_t(state.iterations()) * size);
}
#else
void problem_setup(benchmark::State &state) {
  mpi::Communicator world;
  auto const nHarmonics = state.range_y();
  auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
  auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
  int64_t const size = input.geometry->scatterer_size();

  while(state.KeepRunning()) {
    auto start = std::chrono::high_resolution_clock::now();
    solver::factory(input);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    auto proc_max = world.all_reduce(elapsed_seconds.count(), MPI_MAX);
    state.SetIterationTime(proc_max);
  }
  state.SetItemsProcessed(int64_t(state.iterations()) * size);
}

void benchmark_solver(benchmark::State &state) {
  mpi::Communicator world;
  auto const nHarmonics = state.range_y();
  auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
  auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
  int64_t const size = input.geometry->scatterer_size();
  auto const solver = solver::factory(input);

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
  state.SetItemsProcessed(int64_t(state.iterations()) * size);
}

void serial_matrix_multiplication(benchmark::State &state) {
  mpi::Communicator world;
  auto const nHarmonics = state.range_y();
  auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
  auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
  int64_t const size = input.geometry->scatterer_size();
  auto const Q = source_vector(*input.geometry, input.excitation);
  auto const S = preconditioned_scattering_matrix(*input.geometry, input.excitation);
  Vector<t_complex> result(Q.size());

  while(state.KeepRunning()) {
    auto start = std::chrono::high_resolution_clock::now();
    result = Q * S;
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    auto proc_max = world.all_reduce(elapsed_seconds.count(), MPI_MAX);
    state.SetIterationTime(proc_max);
  }
  state.SetItemsProcessed(int64_t(state.iterations()) * size);
}

#ifdef OPTIMET_SCALAPACK
void scalapack_matrix_multiplication(benchmark::State &state) {
  mpi::Communicator world;
  scalapack::Context context = scalapack::Context::Squarest();
  auto const nHarmonics = state.range_y();
  auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
  auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
  int64_t const size = input.geometry->scatterer_size();
  scalapack::Sizes const block_size = {64, 64};
  auto const Q = distributed_source_vector(source_vector(*input.geometry, input.excitation),
                                           context, block_size);
  auto const S = preconditioned_scattering_matrix(*input.geometry, input.excitation);
  auto const N = input.geometry->scatterer_size();
  scalapack::Matrix<t_complex> Aparallel(context, {N, N}, block_size);
  if(Aparallel.size() > 0)
    Aparallel.local() = S;
  scalapack::Matrix<t_complex> bparallel(context, {N, 1}, block_size);
  if(bparallel.local().size() > 0)
    bparallel.local() = Q;
  scalapack::Matrix<t_complex> result(context, {N, 1}, block_size);
  if(result.local().size() > 0)
    result.local() = Q;

  while(state.KeepRunning()) {
    auto start = std::chrono::high_resolution_clock::now();
    pdgemm(1.5e0, Aparallel, bparallel, 0.0, result);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    auto proc_max = world.all_reduce(elapsed_seconds.count(), MPI_MAX);
    state.SetIterationTime(proc_max);
  }
  state.SetItemsProcessed(int64_t(state.iterations()) * size);
}
#else
void scalapack_matrix_multiplication(benchmark::State &state) {
  throw std::runtime_error("Optimet not compiled with scalapack");
}
#endif

void fmm_multiplication(benchmark::State &state) {
  mpi::Communicator world;
  auto const nHarmonics = state.range_y();
  auto const length = (get_param<t_real>("radius", 0.5) + 0.5) * default_length();
  auto const input = fcc_system(state.range_x(), length, default_scatterer(nHarmonics));
  int64_t const size = input.geometry->scatterer_size();
  mpi::FastMatrixMultiply const fmm(input.geometry->bground, input.excitation->wavenumber(),
                                    input.geometry->objects, input.fmm_subdiagonals, world);
  auto const distribution =
      mpi::details::vector_distribution(input.geometry->objects.size(), world.size());
  auto const first =
      std::find(distribution.data(), distribution.data() + distribution.size(), world.rank()) -
      distribution.data();
  auto const last =
      std::find_if(distribution.data() + first, distribution.data() + distribution.size(),
                   [&world](t_int value) { return value != world.rank(); }) -
      distribution.data();
  auto const Q = source_vector(input.geometry->objects.begin() + first,
                               input.geometry->objects.begin() + last, input.excitation);
  Vector<t_complex> result(Q.size());

  while(state.KeepRunning()) {
    auto start = std::chrono::high_resolution_clock::now();
    fmm(Q, result);
    auto end = std::chrono::high_resolution_clock::now();
    auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    auto proc_max = world.all_reduce(elapsed_seconds.count(), MPI_MAX);
    state.SetIterationTime(proc_max);
  }
  state.SetItemsProcessed(int64_t(state.iterations()) * size);
}

void matrix_multiplication(benchmark::State &state) {
  if(parameters->get<std::string>("Solver") == "eigen")
    return serial_matrix_multiplication(state);
  if(parameters->get<std::string>("Solver") != "eigen" and not parameters->get<bool>("do_fmm"))
    return scalapack_matrix_multiplication(state);
  if(parameters->get<std::string>("Solver") != "eigen" and
     parameters->get<std::string>("Solver") != "scalapack" and parameters->get<bool>("do_fmm"))
    return fmm_multiplication(state);
  throw std::runtime_error("Invalid benchmark request");
}
#endif
}

template <class T> std::vector<T> convert_string(std::string const &input) {
  std::istringstream sstr(input);
  std::vector<T> result;
  T value;
  while(sstr.good()) {
    sstr >> value;
    result.push_back(value);
  }
  return result;
}

static void CustomArguments(benchmark::internal::Benchmark *b) {
  auto const harmonics = convert_string<t_uint>(parameters->get<std::string>("nharmonics"));
  auto const nparticles = convert_string<t_uint>(parameters->get<std::string>("nparticles"));
  for(auto const i : nparticles)
    for(auto const j : harmonics)
      b->ArgPair(i, j);
}

extern std::string FLAG_benchmark_format;

namespace benchmark {
class MPIReporter : public BenchmarkReporter {
public:
  MPIReporter(std::unique_ptr<BenchmarkReporter> reporter, bool doreport)
      : reporter(std::move(reporter)), doreport(doreport) {}
  MPIReporter(std::unique_ptr<BenchmarkReporter> reporter)
      : MPIReporter(std::move(reporter), mpi::Communicator().rank() == 0) {}
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
  mpi::init(argc, const_cast<const char **>(argv));

  ::benchmark::MPIReporter reporter(::benchmark::parse_cmdl(argc, argv));
  ::benchmark::Initialize(&argc, argv);

#ifdef OPTIMET_BELOS
  parameters = parse_cmdl(argc, argv);
#endif

#ifndef OPTIMET_MPI
  BENCHMARK(problem_setup)->Apply(CustomArguments)->Unit(benchmark::kMicrosecond);
  BENCHMARK(benchmark_solver)->Apply(CustomArguments)->Unit(benchmark::kMicrosecond);
#else
  BENCHMARK(matrix_multiplication)
      ->Apply(CustomArguments)
      ->Unit(benchmark::kMillisecond)
      ->UseManualTime();
  BENCHMARK(problem_setup)->Apply(CustomArguments)->Unit(benchmark::kMillisecond)->UseManualTime();
  BENCHMARK(benchmark_solver)
      ->Apply(CustomArguments)
      ->Unit(benchmark::kMillisecond)
      ->UseManualTime();
#endif

  ::benchmark::RunSpecifiedBenchmarks(&reporter);
  mpi::finalize();
  return 0;
}
