//OPTIMET 3D - main.cpp
//
//Restricted to testing for the moment.
//

#include "PeriodicCoupling.h"
#include "Simulation.h"

int main(int argc, char *argv[])
{
	Simulation simulation;
	simulation.init(argv[1]);
	simulation.run();
	simulation.done();
	return 0;
}
