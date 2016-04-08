#include <Types.h>
#ifndef OPTIMET_BELOS
#error This executable needs Belos
#endif
#include <scalapack/Belos.h>
#include "scalapack/LinearSystemSolver.h"

int main(int nargs, char const *args[]) {
  optimet::scalapack::BelosSolverFactory<double> factory;
  auto const avail = factory.supportedSolverNames();
  for(auto const solver_name: avail) {
    std::cout << "Parameters for solver " << solver_name << ":\n";
    Teuchos::RCP<Teuchos::ParameterList> parameters = Teuchos::rcp(new Teuchos::ParameterList());
    auto const solver_manager = factory.create(solver_name, parameters);
    auto const valid_parameters = solver_manager->getValidParameters();
    valid_parameters->print();
    std::cout << "\n";
  }
  return 0;
}
