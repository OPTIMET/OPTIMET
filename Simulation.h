#ifndef SIMULATION_H_
#define SIMULATION_H_

/**
 * The Simulation class implements a full simulation.
 * A Simulation object will create a set of Cases and Requests based
 * on the input file.
 */
class Simulation
{
private:
	char *caseFile; /**< Name of the case without extensions. */
	bool initDone;	/**< Specifies initialization status. */
public:
	/**
	 * Default constructor for the Simulation class.
	 * Does NOT initialize the object.
	 */
	Simulation();

	/**
	 * Initialization constructor for the Simulation class.
	 * @param caseFile the name of the case file (NO extension).
	 */
	Simulation(char *caseFile_);

	/**
	 * Default destructor for the Simulation class.
	 */
	virtual ~Simulation();

	/**
	 * Initialization method for the Simulation class.
	 * @param caseFile_ the name of the case file (NO extension).
	 * @return 0 if successful, 1 otherwise.
	 */
	void init(char *caseFile_);

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
};

#endif /* SIMULATION_H_ */
