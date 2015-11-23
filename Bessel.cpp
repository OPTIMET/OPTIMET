#include "Bessel.h"
#include "constants.h"

#include <iostream>
#include <cmath>

using std::cerr;
using std::abs;
using std::sqrt;

extern "C"
{
  int zbesj_(double *, double *, double *, long int *, long int *, double *, double *, long int *, long int *);
  int zbesh_(double *, double *, double *, long int *, long int *, long int *, double *, double *, long int *, long int *);
}

Bessel::Bessel() {
	initDone = false;
	ierr = 19;
}

Bessel::Bessel(complex<double> argument_, int besselType_, int scale_, int maxOrder_)
{
	init(argument_, besselType_, scale_, maxOrder_);
}

int Bessel::init(complex<double> argument_, int besselType_, int scale_, int maxOrder_)
{
	argument = argument_;
	besselType = besselType_;
	if((besselType < 0) || (besselType > 2))
	{
		cerr << "Wrong besselType:" << besselType <<"!";
		ierr = 11;
		return ierr;
	}
	scale = scale_ + 1;	//Scale must be increased by 1 to correspond to AMOS notations
	maxOrder = maxOrder_;
	if(maxOrder < 1)
	{
		cerr << "Maximum order required for Bessel smaller than 1!";
		ierr = 11;
		return ierr;
	}

	data = new complex<double>[maxOrder+1];
	ddata = new complex<double>[maxOrder+1];

	initDone = true;
	ierr = 0;
	return ierr;
}

int Bessel::populate(void)
{
	if (!initDone) //Parameters were not initialized
	{
		ierr = 10;
		cerr << "Bessel engine was not initialized!";
		return ierr;
	}

	double zr = argument.real();
	double zi = argument.imag();

	double order = 0.5;
	long int noOrders_local = maxOrder+2; //Local no of orders (+1 for zero, +1 for derivative).

	double * cyr = new double[noOrders_local];	//Real part return vector (+1 for derivative)
	double * cyi = new double[noOrders_local];	//Imaginary part return vector (+1 for derivative)

	if (besselType == 0) //Calculate the Bessel function of the first kind
	{
		zbesj_(&zr, &zi, &order, &scale, &noOrders_local, cyr, cyi, &zeroUnderflow, &ierr);
	}
	else //Calculate the Hankel function of the first or second kind
	{
		zbesh_(&zr, &zi, &order, &scale, &besselType, &noOrders_local, cyr, cyi, &zeroUnderflow, &ierr);
	}

	//Assemble the direct functions
	for(int i=0; i<=maxOrder; i++)
	{
		if(abs(argument) <= errEpsilon)
			data[i] = complex<double>(0.0, 0.0);
		else
			data[i] = sqrt(consPi / (2.0*argument)) * complex<double>(cyr[i], cyi[i]);
	}

	//The last derivative
	ddata[maxOrder] = consCm1 * sqrt(consPi / (2.0*argument)) * complex<double>(cyr[maxOrder+1], cyi[maxOrder+1]) + ((double)maxOrder/argument) * data[maxOrder];

	for(int i=0; i<maxOrder; i++)
	{
		ddata[i] = consCm1 * data[i+1] + (((double)i)/argument) * data[i];
	}

	//Clear local vectors
	delete [] cyr;
	delete [] cyi;

	return ierr;
}

Bessel::~Bessel() {
	if(initDone)
	{
		delete [] data;
		delete [] ddata;
	}
}
