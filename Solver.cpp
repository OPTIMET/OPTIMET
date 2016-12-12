#include "Solver.h"

#include "ElectroMagnetic.h"
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
  return std::make_shared<Scalapack>(run);
#elif defined(OPTIMET_BELOS) && defined(OPTIMET_SCALAPACK)
  if(run.belos_params()->template get<std::string>("Solver") == "scalapack")
    return std::make_shared<Scalapack>(run);
  return std::make_shared<MatrixBelos>(run);
#elif defined(OPTIMET_BELOS)
#error FMM will go here
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
