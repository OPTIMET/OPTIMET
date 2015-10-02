#ifndef READER_H_
#define READER_H_

#include "Geometry.h"
#include "Tools.h"
#include "Scatterer.h"
#include "Spherical.h"
#include "Cartesian.h"
#include "ElectroMagnetic.h"
#include "constants.h"
#include "pugi/pugixml.hpp"
#include "Run.h"
#include "aliases.h"
#include <iostream>
#include <cstring>

using namespace pugi;

/**
 * The Reader class is used to read a simulation case file.
 * Provides an interface for reading a Run object.
 */

class Reader
{
private:
	Run *run;				/**< Pointer to the Run object. */
	xml_document inputFile;	/**< The input file. */
	bool initDone;			/**< Specifies if object was initialized. */
public:
	/**
	 * Default constructor for the Reader class.
	 * Does NOT initialize the object.
	 */
	Reader();

	/**
	 * Initializing constructor for the Reader class.
	 * @param geometry_ the pointer to a geometry.
	 */
	Reader(Run *run_);

	/**
	 * Default destructor for the Reader class.
	 */
	virtual ~Reader();

	/**
	 * Initializes the Reader class.
	 * @param geometry_ the pointer to a geometry.
	 */
	void init(Run *run_);

	/**
	 * Read and validate a geometry into the geometry variable.
	 * @return 0 geometry was read and validate, 1 otherwise
	 */
	int readGeometry();

    /**
     * Read an incoming excitation.
     * @return 0 if valid, 1 otherwise.
     */
    int readExcitation();

    /**
     * Read a pre-defined geometric structure.
     * @param geo_node_ the geometry node.
     * @return number of objects pushed, 0 if failure.
     */
    int readStructure(xml_node geo_node_);

	int readSimulation(char* fileName_);

	/**
	 * Read the output requests.
	 * @return 0 if valid, 1 otherwise.
	 */
	int readOutput();
};

#endif /* READER_H_ */
