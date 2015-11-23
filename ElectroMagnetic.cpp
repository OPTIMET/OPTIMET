#include "ElectroMagnetic.h"
#include "constants.h"

//#include "Reader.h"
//using std::cout;
//using std::endl;

ElectroMagnetic::ElectroMagnetic()
{
	init_r(complex<double>(1.0, 0.0), complex<double>(1.0, 0.0));
}

ElectroMagnetic::ElectroMagnetic(complex<double> epsilon_r_,
		complex<double> mu_r_)
{
	init_r(epsilon_r_, mu_r_);
}

ElectroMagnetic::~ElectroMagnetic()
{
	//
}

void ElectroMagnetic::init(complex<double> epsilon_, complex<double> mu_)
{
	epsilon = epsilon_;
	mu = mu_;

	epsilon_r = epsilon / consEpsilon0;
	mu_r = mu / consMu0;

	epsilon_SH = epsilon_;
	mu_SH = mu_;

	epsilon_r_SH = epsilon_SH / consEpsilon0;
	mu_r_SH = mu_SH / consMu0;

	a_SH = complex<double>(0.5, -0.25);
	b_SH = complex<double>(0.1, 0.0);
	d_SH = complex<double>(1.0, 0.0);

	modelType = 0;
}

void ElectroMagnetic::init_r(complex<double> epsilon_r_, complex<double> mu_r_)
{
	epsilon_r = epsilon_r_;
	mu_r = mu_r_;

	epsilon = epsilon_r * consEpsilon0;
	mu = mu_r * consMu0;

	epsilon_r_SH = epsilon_r_;
	mu_r_SH = mu_r_;

	epsilon_SH = epsilon_r_SH * consEpsilon0;
	mu_SH = mu_r_SH * consMu0;

	a_SH = complex<double>(0.5, -0.25);
	b_SH = complex<double>(0.1, 0.0);
	d_SH = complex<double>(1.0, 0.0);

	modelType = 0;
}

// Sellmeier
void ElectroMagnetic::initSellmeier_r(double B1_, double C1_, double B2_,
		double C2_, double B3_, double C3_, double B4_, double C4_, double B5_, double C5_, complex<double> mu_r_)
{
	B1 = B1_;
	C1 = C1_;
	B2 = B2_;
	C2 = C2_;
	B3 = B3_;
	C3 = C3_;
	B4 = B4_;
	C4 = C4_;
	B5 = B5_;
	C5 = C5_;

	mu_r = mu_r_;

	modelType = 1;
}

void ElectroMagnetic::populateSellmeier()
{
	double lambda_aux = lambda * 1e6;

	epsilon_r = ((( (B1 * lambda_aux * lambda_aux) / (lambda_aux * lambda_aux - C1*C1) )
				 +( (B2 * lambda_aux * lambda_aux) / (lambda_aux * lambda_aux - C2*C2) )
				 +( (B3 * lambda_aux * lambda_aux) / (lambda_aux * lambda_aux - C3*C3) )
				 +( (B4 * lambda_aux * lambda_aux) / (lambda_aux * lambda_aux - C4*C4) )
				 +( (B5 * lambda_aux * lambda_aux) / (lambda_aux * lambda_aux - C5*C5) )) + 1);

	epsilon = epsilon_r * consEpsilon0;
}

// Drude Model
void ElectroMagnetic::initDrudeModel_r(complex<double> plasma_freq_, complex<double> damping_freq_, complex<double> mu_r_)
{
	plasma_freq = plasma_freq_;
	damping_freq = damping_freq_;

	mu_r = mu_r_;

	modelType = 2;
}
void ElectroMagnetic::populateDrudeModel()
{
	double input_freq = consC / lambda;
	// epsilon_r = ...
	epsilon_r = 1. - ( (plasma_freq*plasma_freq) / (input_freq*(input_freq+damping_freq)) );

	epsilon = epsilon_r * consEpsilon0;
}

void ElectroMagnetic::update(double lambda_)
{
	lambda = lambda_;

	if(modelType == 0) //Static model
	{
		// Do nothing as we are in the static case.
	}

	if(modelType == 1) //Sellmeier model
	{
		populateSellmeier();
	}

	if(modelType == 2) //Drude model
	{
		populateDrudeModel();
	}
}
