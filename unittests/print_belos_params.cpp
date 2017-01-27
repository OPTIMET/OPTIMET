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
