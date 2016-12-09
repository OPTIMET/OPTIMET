#include "MatrixBelosSolver.h"
#include "scalapack/LinearSystemSolver.h"

namespace optimet {
namespace solver {

void MatrixBelos::solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_,
                        mpi::Communicator const &comm) const {
  if(not context().is_valid())
    return;
  auto const solver = belos_parameters()->get<std::string>("Solver");
  if(solver == "scalapack") {
    Scalapack::solve(X_sca_, X_int_);
    return;
  }
  auto input = parallel_input();
  // Now the actual work
  auto const gls_result = scalapack::gmres_linear_system(std::get<0>(input), std::get<1>(input),
                                                         belos_parameters(), comm);
  if(std::get<1>(gls_result) != 0)
    throw std::runtime_error("Error encountered while solving the linear system");
  // Transfer back to root
  X_sca_ = gather_all_source_vector(std::get<0>(gls_result));
  PreconditionedMatrix::unprecondition(X_sca_, X_int_);
}

void MatrixBelos::solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const {
  auto const solver = belos_parameters()->get<std::string>("Solver");
  return solver == "scalapack" ? Scalapack::solve(X_sca_, X_int_) :
                                 solve(X_sca_, X_int_, mpi::Communicator());
}
}
}
