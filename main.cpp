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
