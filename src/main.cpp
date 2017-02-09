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

// OPTIMET 3D - main.cpp
//
// Restricted to testing for the moment.
//

#include "Simulation.h"
#include "mpi/Session.h"
#include <iostream>

int main(int argc, const char *argv[]) {
  optimet::mpi::init(argc, argv);
  if(argc <= 1) {
    std::cerr << "Usage: " << argv[0] << " <path/to/xml/file/without/extension>" << std::endl;
    return 1;
  }
  optimet::Simulation simulation(argv[1]);
  simulation.run();
  simulation.done();

  optimet::mpi::finalize();

  return 0;
}
