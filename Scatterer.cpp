#include "Scatterer.h"
#include "Tools.h"


Scatterer::Scatterer(Spherical<double> vR_, ElectroMagnetic elmag_, double radius_, int nMax_)
  : vR(vR_), elmag(elmag_), radius(radius_), nMax(nMax_), sourceCoef(2*Tools::iteratorMax(nMax))
{
}

Scatterer::Scatterer(int nMax)
  : Scatterer({0, 0, 0}, {0, 0}, 0e0, nMax)
{
}

Scatterer::~Scatterer()
{
}
