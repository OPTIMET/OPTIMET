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

#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "mpi/Communicator.h"
#include <memory>
#include <string>
#include <vector>

namespace optimet {
class Run;
namespace solver {
class AbstractSolver;
}
/**
 * The Simulation class implements a full simulation.
 * A Simulation object will create a set of Cases and Requests based
 * on the input file.
 */
class Simulation {
public:
  /**
   * Initialization constructor for the Simulation class.
   * @param caseFile the name of the case file (NO extension).
   */
  Simulation(std::string const &filename) : Simulation(filename, mpi::Communicator()) {}
  Simulation(std::string const &filename, mpi::Communicator const &comm)
      : caseFile(filename), communicator_(comm) {}

  /**
   * Default destructor for the Simulation class.
   */
  virtual ~Simulation() {}

  /**
   * Starts a simulation.
   * @return 0 if successful, 1 otherwise.
   */
  int run();

  /**
   * Finishes a simulation.
   * Placeholder method. Not needed at the moment.
   * @return 0 if succesful, 1 otherwise.
   */
  int done();

  mpi::Communicator const &communicator() const { return communicator_; }
  Simulation &communicator(mpi::Communicator const &c) {
    communicator_ = c;
    return *this;
  }

protected:
  #ifdef OPTIMET_MPI
  void scan_wavelengths_parallel(Run &run, std::shared_ptr<solver::AbstractSolver> solver);
  void field_simulation_parallel(Run &run, std::shared_ptr<solver::AbstractSolver> solver);
  void All2all(std::vector<double *> CLGcoeff, std::vector<double *> CLGcoeff_par, int sizeVec);
  #endif
  void scan_wavelengths(Run &run, std::shared_ptr<solver::AbstractSolver> solver);
  void field_simulation(Run &run, std::shared_ptr<solver::AbstractSolver> solver);


private:
  std::string caseFile; /**< Name of the case without extensions. */
  mpi::Communicator communicator_;
};
}
#endif /* SIMULATION_H_ */
