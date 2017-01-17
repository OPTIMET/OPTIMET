#include "MatrixBelosSolver.h"
#include "scalapack/LinearSystemSolver.h"
#include "scalapack/BroadcastToOutOfContext.h"

namespace optimet {
namespace solver {

void MatrixBelos::solve(Vector<t_complex> &X_sca_, Vector<t_complex> &X_int_) const {
  if(belos_parameters()->get<std::string>("Solver", "GMRES") == "scalapack")
    return Scalapack::solve(X_sca_, X_int_);
  auto const splitcomm = communicator().split(context().is_valid());
  if(context().is_valid()) {
    auto const solver = belos_parameters()->get<std::string>("Solver");
    if(solver == "scalapack") {
      Scalapack::solve(X_sca_, X_int_);
      return;
    }
    auto input = parallel_input();
    // Now the actual work
    auto const gls_result = scalapack::gmres_linear_system(std::get<0>(input), std::get<1>(input),
                                                           belos_parameters(), splitcomm);
    if(std::get<1>(gls_result) != 0)
      throw std::runtime_error("Error encountered while solving the linear system");
    // Transfer back to root
    X_sca_ = gather_all_source_vector(std::get<0>(gls_result));
    PreconditionedMatrix::unprecondition(X_sca_, X_int_);
  }
  if(context().size() != communicator().size()) {
    broadcast_to_out_of_context(X_sca_, context(), communicator());
    broadcast_to_out_of_context(X_int_, context(), communicator());
  }
}
}
}
