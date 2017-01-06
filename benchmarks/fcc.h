#ifndef OPTIMET_BENCHMARK_FCC_H
#define OPTIMET_BENCHMARK_FCC_H
#include "Run.h"
#include "Types.h"
#include <Teuchos_ParameterList.hpp>
#ifndef BENCHMARK_HAS_CXX11
#define BENCHMARK_HAS_CXX11
#endif
#include <benchmark/benchmark.h>
#include <chrono>

#ifndef OPTIMET_MPI

#define OPTIMET_BENCHMARK(NAME)                                                                    \
  void NAME(benchmark::State &state, Teuchos::RCP<Teuchos::ParameterList> const &parameters) {     \
    parameters->set("nObjects", state.range(0));                                                   \
    parameters->set("nMax", state.range(1));                                                       \
    auto input = optimet::fcc_input(parameters);                                                   \
    auto const size = input.geometry->scatterer_size();
#define OPTIMET_BENCHMARK_TIME_START while(state.KeepRunning()) {
#define OPTIMET_BENCHMARK_TIME_END                                                                 \
  }                                                                                                \
  state.SetItemsProcessed(int64_t(state.iterations()) * size);                                     \
  state.SetLabel("nprocs(0)");                                                                     \
  }
#define OPTIMET_REGISTER_BENCHMARK(NAME)                                                           \
  auto const NAME##_impl = [parameters](::benchmark::State &state) {                               \
    return NAME(state, parameters);                                                                \
  };                                                                                               \
  ::benchmark::RegisterBenchmark(#NAME, NAME##_impl)                                               \
      ->Apply(NeedExtraDataArgument)                                                               \
      ->Unit(benchmark::kMicrosecond)

#else

#define OPTIMET_BENCHMARK(NAME)                                                                    \
  void NAME(benchmark::State &state, Teuchos::RCP<Teuchos::ParameterList> const &parameters) {     \
    parameters->set("nObjects", state.range(0));                                                   \
    parameters->set("nMax", state.range(1));                                                       \
    auto input = optimet::fcc_input(parameters);                                                   \
    auto const size = input.geometry->scatterer_size();
#define OPTIMET_BENCHMARK_TIME_START                                                               \
  while(state.KeepRunning()) {                                                                     \
    auto start = std::chrono::high_resolution_clock::now();
#define OPTIMET_BENCHMARK_TIME_END                                                                 \
  auto end = std::chrono::high_resolution_clock::now();                                            \
  auto elapsed_seconds = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);   \
  auto proc_max = input.communicator.all_reduce(elapsed_seconds.count(), MPI_MAX);                 \
  state.SetIterationTime(proc_max);                                                                \
  }                                                                                                \
  state.SetItemsProcessed(int64_t(state.iterations()) * size);                                     \
  state.SetLabel("nprocs(" + std::to_string(input.communicator.size()) + ")");                     \
  }

#define OPTIMET_REGISTER_BENCHMARK(NAME)                                                           \
  auto const NAME##_impl = [parameters](::benchmark::State &state) {                               \
    return NAME(state, parameters);                                                                \
  };                                                                                               \
  ::benchmark::RegisterBenchmark(#NAME, NAME##_impl)                                               \
      ->Apply(NeedExtraDataArgument)                                                               \
      ->Unit(benchmark::kMicrosecond)                                                              \
      ->UseManualTime()
#endif

namespace optimet {
Run fcc_input(Teuchos::RCP<Teuchos::ParameterList> const &parameters);
void CustomArguments(std::vector<t_uint> const &nharmonics, std::vector<t_uint> const &nparticles,
                     benchmark::internal::Benchmark *b);
std::unique_ptr<::benchmark::BenchmarkReporter> parse_benchmark_cmdl(int argc, char **argv);

template <class T> std::vector<T> convert_string(std::string const &input) {
  std::string replaced(input);
  for(t_uint i(0); i < replaced.size(); ++i)
    if(replaced[i] == '_')
      replaced[i] = ' ';
  std::istringstream sstr(replaced);
  std::vector<T> result;
  T value;
  while(sstr.good()) {
    sstr >> value;
    result.push_back(value);
  }
  return result;
}

class MPIReporter : public ::benchmark::BenchmarkReporter {
public:
  MPIReporter(std::unique_ptr<::benchmark::BenchmarkReporter> reporter, bool doreport)
      : reporter(std::move(reporter)), doreport(doreport) {}
  MPIReporter(std::unique_ptr<::benchmark::BenchmarkReporter> reporter)
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
  std::unique_ptr<::benchmark::BenchmarkReporter> reporter;
  bool doreport;
};
}

#endif
