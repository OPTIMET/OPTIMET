#include "Scatterer.h"
#include "Tools.h"

Scatterer::Scatterer()
{
	initDone = false;
}

Scatterer::Scatterer(Spherical<double> vR_, ElectroMagnetic elmag_, double radius_, int nMax_)
{
	init(vR_, elmag_, radius_, nMax_);
}

void Scatterer::init(int nMax_)
{
	nMax = nMax_;
	sourceCoef = new complex<double>[2*Tools::iteratorMax(nMax)];
	initDone = true;
}

Scatterer::~Scatterer()
{
	if(initDone)
		delete[] sourceCoef;
}

void Scatterer::init(Spherical<double> vR_, ElectroMagnetic elmag_, double radius_, int nMax_)
{
	vR = vR_;
	elmag = elmag_;
	radius = radius_;
	nMax = nMax_;
	sourceCoef = new complex<double>[2*Tools::iteratorMax(nMax)];
	initDone = true;
}
