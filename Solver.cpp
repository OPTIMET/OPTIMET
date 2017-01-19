#include "Solver.h"

#include "ElectroMagnetic.h"
#include "FMMBelosSolver.h"
#include "MatrixBelosSolver.h"
#include "PreconditionedMatrixSolver.h"
#include "ScalapackSolver.h"
#include "Scatterer.h"
#include "Solver.h"
#include "Types.h"

namespace optimet {
namespace solver {
std::shared_ptr<AbstractSolver> factory(Run const &run) {
#ifndef OPTIMET_MPI
  return std::make_shared<PreconditionedMatrix>(run);
#elif defined(OPTIMET_SCALAPACK) && !defined(OPTIMET_BELOS)
  if(run.do_fmm)
    throw std::runtime_error("Scalapack and Fast Matrix Multiplication are not compatible");
  return std::make_shared<Scalapack>(run);
#elif defined(OPTIMET_BELOS) && defined(OPTIMET_SCALAPACK)
  if((run.belos_params()->get<std::string>("Solver") == "scalapack" or
      run.belos_params()->get<std::string>("Solver") == "eigen") and
     run.do_fmm)
    throw std::runtime_error("Cannot run FMM with scalapack or eigen solver");
  if(run.belos_params()->get<std::string>("Solver") == "eigen")
    return std::make_shared<PreconditionedMatrix>(run);
  if(run.belos_params()->get<std::string>("Solver") == "scalapack")
    return std::make_shared<Scalapack>(run);
  if(run.do_fmm)
    return std::make_shared<FMMBelos>(run);
  return std::make_shared<MatrixBelos>(run);
#elif defined(OPTIMET_BELOS)
  if(run.belos_params()->get<std::string>("Solver") == "scalapack")
    throw std::runtime_error("Optimet was not compiled with scalapack");
  if(not run.do_fmm)
    throw std::runtime_error("Optimet was not compiled with scalapack, please choose FMM matrix");
  return std::make_shared<FMMBelos>(run);
#else
#error Need at least Belos to run MPI solvers
#endif
}
} // namespace solver

Vector<t_complex> convertInternal(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {
  Vector<t_complex> result(scattered.size());
  size_t i = 0;
  for(auto const &object : objects) {
    auto const N = 2 * object.nMax * (object.nMax + 2);
    result.segment(i, N).array() =
        scattered.segment(i, N).array() * object.getIaux(omega, bground).array();
    i += N;
  }
  return result;
}

Vector<t_complex> convertIndirect(Vector<t_complex> const &scattered, t_real const &omega,
                                  ElectroMagnetic const &bground,
                                  std::vector<Scatterer> const &objects) {
  Vector<t_complex> result(scattered.size());
  size_t i(0);
  for(auto const &object : objects) {
    auto const N = 2 * object.nMax * (object.nMax + 2);
    result.segment(i, N) =
        object.getTLocal(omega, bground).array() * scattered.segment(i, N).array();
    i += N;
  }
  return result;
}
} // optimet namespace
