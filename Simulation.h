#ifndef SIMULATION_H_
#define SIMULATION_H_

#include "mpi/Communicator.h"
#include <memory>
#include <string>

class Run;
namespace optimet {
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
  void scan_wavelengths(Run &run, std::shared_ptr<solver::AbstractSolver> solver);
  void field_simulation(Run &run, std::shared_ptr<solver::AbstractSolver> solver);
  void radius_scan(Run &run, std::shared_ptr<solver::AbstractSolver> solver);
  void radius_and_wavelength_scan(Run &run, std::shared_ptr<solver::AbstractSolver> solver);
  void coefficients(Run &run, std::shared_ptr<solver::AbstractSolver> solver);

private:
  std::string caseFile; /**< Name of the case without extensions. */
  //! \details Fake if not compiled with MPI
  mpi::Communicator communicator_;
};
}
#endif /* SIMULATION_H_ */
