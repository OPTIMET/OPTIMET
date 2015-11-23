
#include "AJ_AuxFuns.h"

#include "AuxCoefficients.h"
#include "Bessel.h"
#include "CompoundIterator.h"
// #include "Scatterer.h"
#include "Spherical.h"
#include "constants.h"

#include <cmath>
#include <iostream>
#include <cstdlib>

//#define PI 3.14159265358979323846

// -----------------------------------------------------------------------
// computes all corresponding Spherical Harmonics values given (nMax, m)
int compute_YJn_m(Spherical<double> R, std::complex<double> waveK, int BHreg, int nMax, int m, std::complex<double> *YJnm){


	int i(0), n(0);
	double d_n(0.);
	double d_temp(0.);

	double *dn;
	dn = new double[nMax+1];

	double *Wigner, *dWigner;
	Wigner 	= new double[nMax+1];
	dWigner = new double[nMax+1];

	AuxCoefficients auxCoefficients;
	auxCoefficients.compute_dn(nMax, dn);
	auxCoefficients.VIGdVIG(nMax, m, R, Wigner, dWigner);

	// Initialize and populate Bessel object
	Bessel besselH_R;
	besselH_R.init(R.rrr*waveK, BHreg, 0, nMax);
	besselH_R.populate();

	double dm = pow(-1., double(m));					// Legendre to Wigner function
	std::complex<double> exp_imphi(cos(double(m) * R.phi), sin(double(m) * R.phi));

	for (i = 0; i <= nMax; i++) {

		// obtain Spherical hormonics - Ynm
		n=i;
		d_n=double(i);

		d_temp  = dm*dn[i]*sqrt(d_n*(d_n+1.));

		if(m==0 && n==0){
			YJnm[i]=sqrt(1./(4*consPi));
		}
		else{
			YJnm[i] = d_temp * exp_imphi * Wigner[i];
		}
		// obtain Ynm*Jn
		YJnm[i]*=besselH_R.data[i];

	}

	delete [] dn;
	delete [] Wigner;
	delete [] dWigner;

	return 0;
}
// -----------------------------------------------------------------------



// -----------------------------------------------------------------------
// computes all corresponding Spherical Harmonics values given (nMax) - m=[-nMax, nMax]
int compute_YJp(Spherical<double> R, std::complex<double> waveK, int BHreg, int nMax, std::complex<double> *dataYJp){

	CompoundIterator p;
	CompoundIterator q;
	std::complex<double> *YJnm;
	YJnm = new std::complex<double> [nMax+1];

	// (n,m)=(0,0)
	compute_YJn_m(R, waveK, BHreg, nMax, 0, YJnm);
	dataYJp[p.max(nMax)]=YJnm[0];

	// all other elements compute and map into CompoundIterator format
	for (q = CompoundIterator(nMax, nMax); q < q.max(nMax); q++) {
		for (int n = abs(q.second); n <= nMax; n++) {

			if (n != 0) {
				compute_YJn_m(R, waveK, BHreg, nMax, q.second, YJnm);
				p.init(n, q.second);
				dataYJp[p] = YJnm[n];
			}

		}
	}

	delete [] YJnm;
	return 0;
}
// -----------------------------------------------------------------------


// -----------------------------------------------------------------------
// computes all corresponding Spherical Harmonics values given (nMax, m)
int compute_Yn_m(Spherical<double> R, std::complex<double> waveK, int nMax, int m, std::complex<double> *Ynm){


	int i(0), n(0);
	double d_n(0.);
	double d_temp(0.);

	double *dn;
	dn = new double[nMax+1];

	double *Wigner, *dWigner;
	Wigner 	= new double[nMax+1];
	dWigner = new double[nMax+1];

	AuxCoefficients auxCoefficients;
	auxCoefficients.compute_dn(nMax, dn);
	auxCoefficients.VIGdVIG(nMax, m, R, Wigner, dWigner);

	double dm = pow(-1., double(m));					// Legendre to Wigner function
	std::complex<double> exp_imphi(cos(double(m) * R.phi), sin(double(m) * R.phi));

	for (i = 0; i <= nMax; i++) {

		// obtain Spherical hormonics - Ynm
		n=i;
		d_n=double(i);

		d_temp  = dm*dn[i]*sqrt(d_n*(d_n+1.));

		if(m==0 && n==0){
			Ynm[i]=sqrt(1./(4*consPi));
		}
		else{
			Ynm[i] = d_temp * exp_imphi * Wigner[i];
		}

	}

	delete [] dn;
	delete [] Wigner;
	delete [] dWigner;

	return 0;
}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// computes all corresponding Spherical Harmonics values given (nMax) - m=[-nMax, nMax]
int compute_Yp(Spherical<double> R, std::complex<double> waveK, int nMax, std::complex<double> *dataYp){

	CompoundIterator p;
	CompoundIterator q;
	std::complex<double> *Yn_m;
	Yn_m = new std::complex<double> [nMax+1];

	// (n,m)=(0,0)
	compute_Yn_m(R, waveK, nMax, 0, Yn_m);
	dataYp[p.max(nMax)]=Yn_m[0];

	// all other elements compute and map into CompoundIterator format
	for (q = CompoundIterator(nMax, nMax); q < q.max(nMax); q++) {
		for (int n = abs(q.second); n <= nMax; n++) {

			if (n != 0) {
				compute_Yn_m(R, waveK, nMax, q.second, Yn_m);
				p.init(n, q.second);
				dataYp[p] = Yn_m[n];
			}

		}
	}

	delete [] Yn_m;
	return 0;
}
/*
// T local
int getTLocal(double omega_, Scatterer particle, int nMax_, std::complex<double> ** T_local_)
{

	std::cout<<"So far so good";

	std::complex<double> k_s = omega_ * sqrt(particle.elmag.epsilon * particle.elmag.mu);
	std::complex<double> k_b = omega_ * sqrt(consEpsilon0 * consMu0);

	std::complex<double> rho = k_s / k_b;
	std::complex<double> r_0 = k_b * particle.radius;
	std::complex<double> mu_sob = particle.elmag.mu /consMu0;

	std::complex<double> psi(0., 0.), ksi(0., 0.);
	std::complex<double> dpsi(0., 0.), dksi(0., 0.);

	std::complex<double> psirho(0., 0.), ksirho(0., 0.);
	std::complex<double> dpsirho(0., 0.), dksirho(0., 0.);
	// AJ -----------------------------------------------------------------------------------------------------

	Bessel J_n;
	Bessel Jrho_n;
	Bessel H_n;
	Bessel Hrho_n;

	J_n.init(r_0, 0, 0, nMax_);
	if(J_n.populate())
	{
		std::cerr << "Error computing Bessel functions. Amos said: " << J_n.ierr << "!";
		return 1;
	}

	Jrho_n.init(rho*r_0, 0, 0, nMax_);
	if(Jrho_n.populate())
	{
		std::cerr << "Error computing Bessel functions. Amos said: " << Jrho_n.ierr << "!";
		return 1;
	}

	H_n.init(r_0, 1, 0, nMax_);
	if(H_n.populate())
	{
		std::cerr << "Error computing Hankel functions. Amos said: " << H_n.ierr << "!";
	}

	Hrho_n.init(rho*r_0, 1, 0, nMax_);
	if(Hrho_n.populate())
	{
		std::cerr << "Error computing Hankel functions. Amos said: " << Hrho_n.ierr << "!";
	}

	CompoundIterator p;
	CompoundIterator q;

	int pMax = p.max(nMax_); //Recalculating these would be pointless.
	int qMax = q.max(nMax_);

	for(p=0; (int)p<pMax; p++)
		for(q=0; (int)q<qMax; q++)
		{
			if((p.first == q.first) && (p.second == q.second)) //Kronicker symbols
			{

				// AJ ------------------------------------------------------------
				// obtain aux functions
				psi = r_0*J_n.data[p.first];
				dpsi= r_0*J_n.ddata[p.first] + J_n.data[p.first];

				ksi = r_0*H_n.data[p.first];
				dksi= r_0*H_n.ddata[p.first] + H_n.data[p.first];

				psirho = r_0*rho*Jrho_n.data[p.first];
				dpsirho= r_0*rho*Jrho_n.ddata[p.first] + Jrho_n.data[p.first];

				ksirho = r_0*rho*Hrho_n.data[p.first];
				dksirho= r_0*rho*Hrho_n.ddata[p.first] + Hrho_n.data[p.first];

				//TE Part
				T_local_[p][q] = (psi/ksi) 		* (mu_sob*dpsi/psi - rho*dpsirho/psirho)
												/ (rho*dpsirho/psirho - mu_sob*dksi/ksi);

				//TM part
				T_local_[(int)p+pMax][(int)q+qMax] = (psi/ksi)	* (mu_sob*dpsirho/psirho - rho*dpsi/psi)
																/ (rho*dksi/ksi - mu_sob*dpsirho/psirho);
				// AJ ------------------------------------------------------------
			}
			else
			{
				T_local_[p][q] = std::complex<double>(0.0, 0.0);
				T_local_[(int)p+pMax][(int)q+qMax] == std::complex<double>(0.0, 0.0);
			}

			//Set the rest of the matrix to (0.0, 0.0)
			T_local_[(int)p+pMax][q] = std::complex<double>(0.0, 0.0);
			T_local_[p][(int)q+qMax] = std::complex<double>(0.0, 0.0);
		}

	return 0;
}
*/
