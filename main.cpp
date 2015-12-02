//OPTIMET 3D - main.cpp
//
//Restricted to testing for the moment.
//

#include "Simulation.h"
#include <iostream>

int main(int argc, char *argv[])
{
	if (argc < 1)
	{
		std::cerr << "Usage: " << argv[0] << " <path/to/xml/file/without/extension>" << std::endl;
		return 1;
	}
	Simulation simulation(argv[1]);
	simulation.run();
	simulation.done();
	return 0;
}
