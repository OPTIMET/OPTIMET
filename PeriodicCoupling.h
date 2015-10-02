/*
 * Periodic_Couplin.h
 *
 *  Created on: June 27, 2013
 *      Author: uceealj
 */

#ifndef Periodic_Coupling_H_
#define Periodic_Coupling_H_

#include "AuxCoefficients.h"
#include "Spherical.h"
#include "Tools.h"
#include "Legendre.h"
#include "Bessel.h"
#include "Symbol.h"
#include "CompoundIterator.h"
#include "Excitation.h"
#include "constants.h"
#include "gsl/gsl_sf_gamma.h"
#include <complex>
#include <cmath>
#include <iostream>
#include "AJ_AuxFuns.h"

class PeriodicCoupling
{
private:

public:
	static void compute_AlBe_nmlk(Spherical<double> R, complex<double> waveK, int BHreg, int n_max, complex<double> ****AlBe_nmlk);

	static complex<double> compute_Bnm(int n1, int n2, int n3, int m1, int m2, int m3);

	static int compute_Znmlk(double a1, double a2, complex<double> waveK, int Rmax, int n_max, complex<double> ****Z_nmlk);

	static int compute_OMEGAnmlk(complex<double> waveK, int n_max, complex<double> **dataOMEGA_1pq, complex<double> **dataOMEGA_2pq);

	static int compute_Ap_vkg(int n_max, Cartesian<double> vecKg, complex<double> waveK, double Ao, Cartesian<complex<double> > *dataAp_vkg);

	static double a_nm_p(double n, double m);
	static double a_nm_m(double n, double m);
	static double b_nm_p(double n, double m);
	static double b_nm_m(double n, double m);

	static double compute_alpha_nm(int n, int m);
	static double compute_beta_nm(int n, int m);
};

#endif /*Periodic_Coupling_H_ */
