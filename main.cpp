
#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <time.h>
#include <complex>
#include "Spherical.h"
#include "AuxCoefficients.h"
#include "Bessel.h"
#include "Reader.h"
#include "CompoundIterator.h"
#include "Legendre.h"
#include "Excitation.h"
#include "Coupling.h"
#include "SphericalP.h"
#include "Solver.h"
#include "Result.h"
#include "Symbol.h"
#include "Tools.h"

#include "AJ_AuxFuns.h"
#include "PeriodicCoupling.h"
#include "Simulation.h"

using namespace std;

int getTLocal(double omega_, Scatterer particle, int nMax_, complex<double> ** T_local_)
{

	complex<double> k_s = omega_ * sqrt(particle.elmag.epsilon * particle.elmag.mu);
	complex<double> k_b = omega_ * sqrt(consEpsilon0 * consMu0);

	complex<double> rho = k_s / k_b;
	complex<double> r_0 = k_b * particle.radius;
	complex<double> mu_sob = particle.elmag.mu /consMu0;

	complex<double> psi(0., 0.), ksi(0., 0.);
	complex<double> dpsi(0., 0.), dksi(0., 0.);

	complex<double> psirho(0., 0.), ksirho(0., 0.);
	complex<double> dpsirho(0., 0.), dksirho(0., 0.);
	// AJ -----------------------------------------------------------------------------------------------------

	Bessel J_n;
	Bessel Jrho_n;
	Bessel H_n;
	Bessel Hrho_n;

	J_n.init(r_0, 0, 0, nMax_);
	if(J_n.populate())
	{
		cerr << "Error computing Bessel functions. Amos said: " << J_n.ierr << "!";
		return 1;
	}

	Jrho_n.init(rho*r_0, 0, 0, nMax_);
	if(Jrho_n.populate())
	{
		cerr << "Error computing Bessel functions. Amos said: " << Jrho_n.ierr << "!";
		return 1;
	}

	H_n.init(r_0, 1, 0, nMax_);
	if(H_n.populate())
	{
		cerr << "Error computing Hankel functions. Amos said: " << H_n.ierr << "!";
	}

	Hrho_n.init(rho*r_0, 1, 0, nMax_);
	if(Hrho_n.populate())
	{
		cerr << "Error computing Hankel functions. Amos said: " << Hrho_n.ierr << "!";
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
				T_local_[p][q] = complex<double>(0.0, 0.0);
				T_local_[(int)p+pMax][(int)q+qMax] == complex<double>(0.0, 0.0);
			}

			//Set the rest of the matrix to (0.0, 0.0)
			T_local_[(int)p+pMax][q] = complex<double>(0.0, 0.0);
			T_local_[p][(int)q+qMax] = complex<double>(0.0, 0.0);
		}
	return 0;
}

int getIaux (double omega_, Scatterer particle, int nMax_, complex<double> *I_aux_){

	complex<double> k_s = omega_ * sqrt(particle.elmag.epsilon * particle.elmag.mu);
	complex<double> k_b = omega_ * sqrt(consEpsilon0 * consMu0);

	complex<double> rho = k_s / k_b;
	complex<double> r_0 = k_b * particle.radius;
	complex<double> mu_j  = particle.elmag.mu;
	complex<double> mu_0  = consMu0;
	complex<double> mu_sob = particle.elmag.mu /consMu0;

	complex<double> psi(0., 0.), ksi(0., 0.);
	complex<double> dpsi(0., 0.), dksi(0., 0.);

	complex<double> psirho(0., 0.);
	complex<double> dpsirho(0., 0.);

	Bessel J_n;
	Bessel Jrho_n;

	J_n.init(r_0, 0, 0, nMax_);
	if(J_n.populate())
	{
		cerr << "Error computing Bessel functions. Amos said: " << J_n.ierr << "!";
		return 1;
	}

	Jrho_n.init(rho*r_0, 0, 0, nMax_);
	if(Jrho_n.populate())
	{
		cerr << "Error computing Bessel functions. Amos said: " << Jrho_n.ierr << "!";
		return 1;
	}

	// AJ - direct
	Bessel H_n;
	H_n.init(r_0, 1, 0, nMax_);
	if(H_n.populate())
	{
		cerr << "Error computing Hankel functions. Amos said: " << H_n.ierr << "!";
	}

	CompoundIterator p;

	int pMax = p.max(nMax_);

	/**
	 * @todo Modify this for non-spherical objects (need matrix, also see Solver::solveInternal).
	 */

	for(p=0; (int)p<pMax; p++)
	{
			// obtain Riccati-Bessel functions
			psi = r_0*J_n.data[p.first];
			dpsi= r_0*J_n.ddata[p.first] + J_n.data[p.first];

			psirho = r_0*rho*Jrho_n.data[p.first];
			dpsirho= r_0*rho*Jrho_n.ddata[p.first] + Jrho_n.data[p.first];


			// Stout 2002 - from scattered
			//TE Part
			I_aux_[p] 			= (mu_j*rho) / (mu_0*rho*dpsirho*psi - mu_j*psirho*dpsi);
			I_aux_[p] 		   *= complex<double>(0., 1.);

			//TM part
			I_aux_[(int)p+pMax]	= (mu_j*rho) / (mu_j*psi*dpsirho - mu_0*rho*psirho*dpsi);
			I_aux_[(int)p+pMax]*= complex<double>(0., 1.);
	}
	return 0;
}


int getCabsAux (double omega_, Scatterer particle, int nMax_, complex<double> *Cabs_aux_){

    complex<double> temp1(0., 0.), temp2(0., 0.);
    complex<double> k_s = omega_ * sqrt(particle.elmag.epsilon * particle.elmag.mu);
    complex<double> k_b = omega_ * sqrt(consEpsilon0 * consMu0);

    complex<double> rho = k_s / k_b;
    complex<double> r_0 = k_b * particle.radius;
    complex<double> mu_j  = particle.elmag.mu;
    complex<double> mu_0  = consMu0;
    complex<double> mu_sob = particle.elmag.mu /consMu0;

    complex<double> psi(0., 0.), ksi(0., 0.);
    complex<double> dpsi(0., 0.), dksi(0., 0.);

    complex<double> psirho(0., 0.);
    complex<double> dpsirho(0., 0.);

    Bessel J_n;
    Bessel Jrho_n;

    J_n.init(r_0, 0, 0, nMax_);
    if(J_n.populate())
    {
        cerr << "Error computing Bessel functions. Amos said: " << J_n.ierr << "!";
        return 1;
    }

    Jrho_n.init(rho*r_0, 0, 0, nMax_);
    if(Jrho_n.populate())
    {
        cerr << "Error computing Bessel functions. Amos said: " << Jrho_n.ierr << "!";
        return 1;
    }

    // AJ - direct
    Bessel H_n;
    H_n.init(r_0, 1, 0, nMax_);
    if(H_n.populate())
    {
        cerr << "Error computing Hankel functions. Amos said: " << H_n.ierr << "!";
    }

    CompoundIterator p;

    int pMax = p.max(nMax_);

    /**
     * @todo Modify this for non-spherical objects (need matrix, also see Solver::solveInternal).
     */

    for(p=0; (int)p<pMax; p++)
    {
            // obtain Riccati-Bessel functions
            psi = r_0*J_n.data[p.first];
            dpsi= r_0*J_n.ddata[p.first] + J_n.data[p.first];

            psirho = r_0*rho*Jrho_n.data[p.first];
            dpsirho= r_0*rho*Jrho_n.ddata[p.first] + Jrho_n.data[p.first];

            // Stout 2002 - from scattered
            //TE Part
            temp1 = complex<double>(0., 1.)*rho*mu_0*conj(mu_j)*conj(psirho)*dpsirho;
            temp2 = abs((mu_j*psirho*dpsi - mu_0*rho*dpsirho*psi));
            temp2*=temp2;
            Cabs_aux_[p]             = real(temp1)/temp2;

            //TM part
            temp1 = complex<double>(0., 1.)*conj(rho)*mu_0*mu_j*conj(psirho)*dpsirho;
            temp2 = abs((mu_0*rho*psirho*dpsi - mu_j*dpsirho*psi));
            temp2*=temp2;
            Cabs_aux_[(int)p+pMax]    = real(temp1)/temp2;
    }
    return 0;
}

/*
int checkInner(Spherical<double> R_, Scatterer *particles, int np)
{
	Spherical<double> Rrel;

	for(int j=0; j<np; j++)
	{
		//Translate to object j
		Rrel = Tools::toPoint(R_, particles[j].vR);

		//Check if R_.rrr is inside radius of object j
		if(Rrel.rrr <= particles[j].radius)
			return j;
		else
			return -1;
	}
}
*/

int checkInner(Spherical<double> R_, Scatterer *particles, int Np) {
	int j, Np_in = -1;						// assume it is always outside
	Spherical<double> Rrel;

	for(j=0; j<Np; j++){
		//Translate to object j
		Rrel = Tools::toPoint(R_, particles[j].vR);
		//Check if R_.rrr is inside radius of object j
		if(Rrel.rrr <= (particles[j].radius))
			Np_in=j;
	}
	return Np_in;

}

int main(int argc, char *argv[])
{


	// Ahmed -------------------------------------------------------------------------------
	int nMax=20;							// max no of harmonics
	int Obs_part=0;							// cont check - choose observation particle
	int steps = 200;						// op steps for 2D plot
	int the_range = 180;					// cont check steps
	int phi_range = 180;					// cont check steps
	int Cext_nMax=15;						// individual TE/TM components to calculate

	// scan parameters ---------------------------------------------------------------------
	double inc =1;
	// lam based loop
	double lam_start=390.;							// Cext scan start at
	double lam_final=5390.;							// Cext scan finish at
	int Cext_steps=lam_final-lam_start+1;						// Cext steps - user defined
	double lam_steps=((lam_final-lam_start)/double(Cext_steps-1));	// Cext scan steps
	// radius based loop
	double rad_start=10.;							// Cext scan start at
	double rad_final=1000.;							// Cext scan finish at
	double rad_steps=((rad_final-rad_start)/double(Cext_steps-1));	// Cext scan steps
	// -------------------------------------------------------------------------------------

	// Translation Coefficients test - from Periodic_Coupling.cpp + Coupling.cpp
	int l, i, j, ii, jj;
	complex<double> c_temp(0., 0.);
	CompoundIterator p, q;
	int pMax = p.max(nMax);
	int BHreg_sca=0;					// compute regular values of BH - non-regular
	int BHreg_inc=1;					// compute regular values of BH - regular

	// input variables
	double lam_0=500*1e-9;				// input wavelength
	double lam_ip=lam_0;				// default scan
	complex<double> eps_i=complex<double>(13, 0.);
	double freq = consC/lam_0;
	double omega=2*consPi*freq;
	complex<double> waveK_0 = omega*sqrt(consEpsilon0 * consMu0);	// background wave-number
	complex<double> waveK_j;
	double x(0.), y(0.), z(0.);

	cout<<"AJ intake on OPTIMET-3D -------------------- "<<endl;
	cout<<"Scan Cext_setps = "<<Cext_steps<<endl;

	// AJ - Chiral Auto build --------------------------------------------------
	//	Spiral structure building block
	double radius = 500.0;
	double rad_ip=radius;				// default scan
	double sep =100.0;
	int arms =3;							// Number of arms in spiral structure
	int No = 4; 							// Number of particles on each arm including center
	//Distance based simulation
	double Theta = consPi/(No-1); 			// Calculate the separation angle
	double d = 2.0*radius + sep;			// Separation between two particles
	double Phi = (consPi-Theta)/2.0;		// remaining angles on equal sided triangle
	double R = d * sin(Phi)/sin(Theta);		// Determine R - given by triangle rule : d/sin(Theta) = R/sin(Phi)
	if(No==2)								// Special case where Theta==consPi
		R=d;

	int Np =(No-1)*arms + 1;				// total number of particles in system
	Np=2;
	cout<<"Total number of particles = "<<Np<<endl;
	Scatterer *particles;
	particles = new Scatterer[Np];
	Spherical<double> Origin=Spherical<double>(0., 0., 0.);				// cluster origin
	Origin = Spherical<double>(Tools::toSpherical(Cartesian<double>(0., 0., 0.)));

	double Theta_rot = 2*consPi/arms;
	// 1 - determine fundamental arm particles' locations
	double x_arm[No-1]; // to store x-axis location of particles for each arm
	double y_arm[No-1]; // to store y-axis location of particles for each arm
	j=0;
	for(i=0; i<No-1; i++){
		double x_;
		double y_;
		Tools::Pol2Cart(R, i*Theta, x_, y_);
		x_arm[j] = x_ + R;
		y_arm[j] = y_;
		j++;
	}

	//Create vectors for r, theta, x and y, X and Y
	double X_= 0.;							// to store x-axis location of particles for all arm
	double Y_= 0.;							// to store y-axis location of particles for all arm
	eps_i=complex<double>(1.61, 0.);			// fixed here - Sellmier is used below
	// particle 1
	particles[0].vR = Spherical<double>(Tools::toSpherical(Cartesian<double>(0., 0., -radius*1e-9)));	// The coordinates of the center of the scatterer
	particles[0].elmag.epsilon=consEpsilon0*eps_i;						// The electromagnetic properties of the scatterer
	particles[0].elmag.mu=consMu0;										// The electromagnetic properties of the scatterer.
	particles[0].nMax=nMax;												// Maximum value of the n iterator.
	particles[0].radius=radius*1e-9;									// The radius of a sphere encompassing the scatterer.
	// particle 2
	particles[1].vR = Spherical<double>(Tools::toSpherical(Cartesian<double>(0., 0., radius*1e-9)));	// The coordinates of the center of the scatterer
	particles[1].elmag.epsilon=consEpsilon0*eps_i;						// The electromagnetic properties of the scatterer
	particles[1].elmag.mu=consMu0;										// The electromagnetic properties of the scatterer.
	particles[1].nMax=nMax;												// Maximum value of the n iterator.
	particles[1].radius=radius*1e-9;									// The radius of a sphere encompassing the scatterer.
	// AJ ------------------------------------------------------------------------------------------------------------------

	// 2 - Global incident field -------------------------------------------------------------------------------------------
	complex<double> *Inc;
	Inc = new complex<double> [2*p.max(nMax)];
	SphericalP<complex<double> > Einc(complex<double>(1.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
	Spherical<double> vKInc;
	vKInc.rrr = waveK_0.real();					// stores incoming wavenumber
	// z-axis propagation
	vKInc.the = consPi*(0.0001/180.);
	vKInc.phi = consPi*(0.0001/180.);
	// y-axis propagation - does not exist for LG beams
	// x-axis propagation

	// Set Eaux Polarisation -----------------------------------------------------------------------------------------------------
	Spherical<double> vAux = Spherical<double>(0.0, vKInc.the, vKInc.phi);
	SphericalP<complex<double> > Eaux;
	// linear
	Eaux = SphericalP<complex<double> >(complex<double>(0.0, 0.0),
										complex<double>(0.0, 0.0),
										complex<double>(1.0, 0.0));
	// circular
	int n_op_start=1;
	int n_op_final=nMax;
	Einc = Tools::toProjection(vAux, Eaux);
	// ---------------------------------------------------------------------------------------------------------------------------

	Excitation excitation;
	int Type=0;														// compute regular values of BH - regular
	excitation.init(Type, Einc, vKInc, nMax);
	excitation.populate();
	for(i=0; i<p.max(nMax); i++){
		Inc[i]=excitation.dataIncAp[i];
		Inc[i+1*p.max(nMax)]=excitation.dataIncBp[i];
	}

	// latest changes -------------------------------------------------------------------------
	// System of equations vectors ------------------------------------------------------------
	complex<double> *F;
	F = new complex<double> [Np*2*pMax];							// Scattered coeffs solution

	ofstream op_CextCabs("op_CextCabs.txt");
	ofstream op_CextCabs_EH("op_CextCabs_EH.txt");
	ofstream op_CextCabs_EH_1("op_CextCabs_EH_1.txt");
	ofstream op_CextCabs_EH_2("op_CextCabs_EH_2.txt");
	ofstream op_CextCabs_EH_3("op_CextCabs_EH_3.txt");
	ofstream op_CextCabs_EH_4("op_CextCabs_EH_4.txt");
	ofstream op_CextCabs_EH_5("op_CextCabs_EH_5.txt");
	time_t itime, ftime;
	double drift(0.);
	int hours(0), mins(0), secs(0);
	// latest changes -------------------------------------------------------------------------

	// -----------------------------------------------------------------------------------------
	for(l=0, lam_ip=lam_start; l<Cext_steps; l++, lam_ip+=lam_steps){

		// 3 - obtain Translation coeffs matrices - alpha & beta -------------------------------------------------------
		Spherical<double> R_rel;
		Coupling TC_alpha;										// alpha terms
		Coupling TC_beta;										// beta terms - g==origin

		// 4 - construct system of equations ------------------------------------------------------
		// 4.1 individual terms -------------------------------------------------------------------
		complex<double> *B0;
		B0 = new complex<double> [2*pMax];								// individual Incident coeffs
		// -----------------------------------------------------------------------------------------
		complex<double> **T_matrix;										// individual T-matrix
		T_matrix = Tools::Get_2D_c_double(2*p.max(nMax), 2*p.max(nMax));
		// ----------------------------------------------------------------------------------------

		// 4.2 overall system of equations --------------------------------------------------------
		// Auxiliary TC arrays --------------------------------------------------------------------
		complex<double> **TC_betaa_t, **TC_alpha_t;
		complex<double> **mT_TC_alpha;
		TC_betaa_t = Tools::Get_2D_c_double(2*pMax, 2*pMax);
		TC_alpha_t = Tools::Get_2D_c_double(2*pMax, 2*pMax);
		mT_TC_alpha = Tools::Get_2D_c_double(2*pMax, 2*pMax);
		// -----------------------------------------------------------------------------------------
		complex<double> *B, *E;// *F;
		B = new complex<double> [Np*2*pMax];							// translated Incident coeffs
		E = new complex<double> [Np*2*pMax];							// Excitation coeffs solution
		// -----------------------------------------------------------------------------------------
		// System of equations global matrices -----------------------------------------------------
		complex<double> **SS, **TT;
		SS = Tools::Get_2D_c_double(Np*2*pMax, Np*2*pMax);
		TT = Tools::Get_2D_c_double(Np*2*pMax, Np*2*pMax);
		// -----------------------------------------------------------------------------------------

		// start scan here - Cext, Cabs ------------------------------------------------------------
		complex<double> Cext(0.0, 0.0);
		complex<double> Cext1(0.0, 0.0);
		complex<double> Cext2(0.0, 0.0);
		complex<double> Cabs(0.0, 0.0);
		complex<double> Csca(0.0, 0.0);
		complex<double> Csca_validate(0.0, 0.0);

		complex<double> Qext(0.0, 0.0);
		complex<double> Qabs(0.0, 0.0);
		complex<double> Qsca(0.0, 0.0);
		// computed from normalised internal like & scattered coeffs
		complex<double> *Cabs_aux;
		Cabs_aux = new complex<double>[2*pMax];

		// change operating lam_0
		lam_0=lam_ip*1e-9;
		freq = consC/lam_0;
		omega=2*consPi*freq;
		waveK_0 = omega*sqrt(consEpsilon0*consMu0);

		cout<<"compute lam = "<<lam_0<<endl;

		// Sellmier dispersion relationship
		for (j=0; j<Np; j++){											// loop over all particles
			double lam_sl = lam_0*1e6;
			eps_i = ( 1 + 10.6684293*pow(lam_sl,2)/(pow(lam_sl,2)-pow(0.301516485,2)) + 0.003043475*pow(lam_sl,2)/(pow(lam_sl,2)-pow(1.13475115,2)) + 1.54133408*pow(lam_sl,2)/(pow(lam_sl,2)-pow(1104.0,2)) );
			eps_i = 1.61*1.61;
			particles[j].elmag.epsilon=consEpsilon0*eps_i;				// The electromagnetic properties of the scatterer
		}
	//	particles[1].vR = Spherical<double>(Tools::toSpherical(Cartesian<double>(0., 0., -rad_ip*1e-9)));	// The coordinates of the center of the scatterer

		if (l==0) cout<<"1 - Build system of equations 	: E = [SS].B"<<endl;
		time(&itime);                                       // calculate current time
		// 4.3 Build system of equations -----------------------------------------------------------
		for (i=0; i<Np; i++){											// loop over all particles
			for (j=0; j<Np; j++){										// loop over all particles

				// get local T-matrix ---------------------------------------------------------------
				getTLocal(omega, particles[j], nMax, T_matrix);

				if(i==j){

					// A - Beta terms ---------------------------------------------------------------
					R_rel = particles[i].vR - Origin;
					TC_beta.init(R_rel, waveK_0, BHreg_inc, nMax);
					TC_beta.populate();
					// ------------------------------------------------------------------------------

					// B - map beta terms obtained from Coupling.cpp to account for efficient multiplication
					for(ii=0; ii<pMax; ii++){
						for(jj=0; jj<pMax; jj++){
							// beta terms -----------------------------------------------------------
							// - Q11
							TC_betaa_t[ii][jj]					= TC_beta.dataApq[jj][ii];
							// - Q12
							TC_betaa_t[ii][jj+pMax]				= TC_beta.dataBpq[jj][ii];
							// - Q21
							TC_betaa_t[ii+pMax][jj]				= TC_beta.dataBpq[jj][ii];
							// - Q22
							TC_betaa_t[ii+pMax][jj+pMax]		= TC_beta.dataApq[jj][ii];
							// ----------------------------------------------------------------------
						}
					}
					// ------------------------------------------------------------------------------

					// C - incident coefficients translation ----------------------------------------
					// compute B0 - Implements the multiplication Y = alpha*A*X + beta*Y ------------
					Algebra::multiplyVectorMatrix(	TC_betaa_t, 		2*pMax,			2*pMax,
													Inc,				B0, 			1., 	0.);
					// populate the system translated incident coefficients -------------------------
					for(ii=0; ii<2*pMax; ii++){
						B[j*2*pMax+ii]=B0[ii];
					}
					// ------------------------------------------------------------------------------

					// D - populate global scattering matrix ----------------------------------------
					// identity matrix --------------------------------------------------------------
					for(ii=0; ii<2*pMax; ii++){
						for(jj=0; jj<2*pMax; jj++){
							// SS matrix ------------------------------------------------------------
							if(ii==jj)
								SS[ii+i*2*pMax][jj+j*2*pMax]	= complex<double>(1., 0.);	// diagonal terms
							else
								SS[ii+i*2*pMax][jj+j*2*pMax]	= complex<double>(0., 0.);	// everything else
							// SS matrix ------------------------------------------------------------

							// TT matrix ------------------------------------------------------------
							TT[ii+i*2*pMax][jj+j*2*pMax] 		= T_matrix[ii][jj];
							// TT matrix ------------------------------------------------------------
						}
					}
					// -----------------------------------------------------------------------------
				}	// if(i==j)

				else{
					// A - Alpha terms --------------------------------------------------------------
					R_rel = particles[i].vR - particles[j].vR;
					TC_alpha.init(R_rel, waveK_0, BHreg_sca, nMax);
					TC_alpha.populate();
					// ------------------------------------------------------------------------------

					// B - map aplha terms obtained from Coupling.cpp to account for efficient multiplication
					for(ii=0; ii<pMax; ii++){
						for(jj=0; jj<pMax; jj++){
							// alpha terms ----------------------------------------------------------
							// - Q11
							TC_alpha_t[ii][jj]					= TC_alpha.dataApq[jj][ii];
							// - Q12
							TC_alpha_t[ii][jj+pMax]				= TC_alpha.dataBpq[jj][ii];
							// - Q21
							TC_alpha_t[ii+pMax][jj]				= TC_alpha.dataBpq[jj][ii];
							// - Q22
							TC_alpha_t[ii+pMax][jj+pMax]		= TC_alpha.dataApq[jj][ii];
							// ----------------------------------------------------------------------
						}
					}

					// C - scattered coefficients translation ---------------------------------------
					// obtain multiplication : -1[T_][TC_alpha] -------------------------------------
					Algebra::multiplyMatrixMatrix(	TC_alpha_t, 		2*pMax,			2*pMax,
													T_matrix, 			2*pMax,			2*pMax,
													mT_TC_alpha, 		-1., 			0.);
					// ------------------------------------------------------------------------------

					// D - populate global scattering matrix ----------------------------------------
					// identity matrix --------------------------------------------------------------
					for(ii=0; ii<2*pMax; ii++){
						for(jj=0; jj<2*pMax; jj++){
							// SS matrix ------------------------------------------------------------
							SS[ii+i*2*pMax][jj+j*2*pMax]		= mT_TC_alpha[ii][jj];
							// TT matrix ------------------------------------------------------------
							TT[ii+i*2*pMax][jj+j*2*pMax] 		= complex<double>(0., 0.);
						}
					}
					// ------------------------------------------------------------------------------
				}	// else

			}
		}
		// op time ----------------------------------------------------------------------------------
		time(&ftime);                                       // calculate final time
		drift = difftime(ftime, itime);
		hours=int(drift/3600);
		mins =int(60*((drift/3600)-hours));
		secs =int(60*(60*((drift/3600)-hours)-mins));
		if (l==0) cout<<"1 - Build system of equations execution time\t"<<hours<<" hrs "<<mins<<" mins "<<secs<<" secs" << endl;
		//-------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------
		if (l==0) cout<<"2 - Solve system of equations 	: E = [SS].B"<<endl;
		// 4.4 excitation coeffs solution - solve E = [SS].B ----------------------------------------
		time(&itime);                                       // calculate initial time
		AlgebraS::solveMatrixVector(		SS, 		Np*2*pMax, 		Np*2*pMax, B, E);
		// op time --------------------
		time(&ftime);                                       // calculate final time
		drift = difftime(ftime, itime);
		hours=int(drift/3600);
		mins =int(60*((drift/3600)-hours));
		secs =int(60*(60*((drift/3600)-hours)-mins));
		if (l==0) cout<<"2 - Solve system of equations execution time\t"<<hours<<" hrs "<<mins<<" mins "<<secs<<" secs" << endl;
		//-------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------
		time(&itime);                                       // calculate initial time
		if (l==0) cout<<"3 - compute                    	: F = [TT].E"<<endl;
		// 4.5 matrix vector multiplication to finally obtain F : F = [TT].E ------------------------
		Algebra::multiplyVectorMatrix(	TT, 		Np*2*pMax,		Np*2*pMax,
											E,				F, 			1., 	0.);
		// op time --------------------
		time(&ftime);                                       // calculate final time
		drift = difftime(ftime, itime);
		hours=int(drift/3600);
		mins =int(60*((drift/3600)-hours));
		secs =int(60*(60*((drift/3600)-hours)-mins));
		if (l==0) cout<<"3 - compute system of equations execution time\t"<<hours<<" hrs "<<mins<<" mins "<<secs<<" secs" << endl;
		// ------------------------------------------------------------------------------------------

		// Calculate cross sections Cext, Cabs, Csca ------------------------------------------------------------------------------------------
		// 7 - computed from incident & scattered coeffs
		Cext=complex<double>(0.0, 0.0);
		Cext1=complex<double>(0.0, 0.0);
		Cext2=complex<double>(0.0, 0.0);
		Cabs=complex<double>(0.0, 0.0);
		Csca=complex<double>(0.0, 0.0);
		Csca_validate=complex<double>(0.0, 0.0);
		Qext=complex<double>(0.0, 0.0);
		Qabs=complex<double>(0.0, 0.0);
		Qsca=complex<double>(0.0, 0.0);
		std::vector<complex<double> > Cext_E(Cext_nMax+1);
		std::vector<complex<double> > Cext_H(Cext_nMax+1);
		int cext_opi=0;
		for(cext_opi=0; cext_opi<=Cext_nMax; cext_opi++){
			Cext_E[cext_opi]=complex<double>(0.0, 0.0);
			Cext_H[cext_opi]=complex<double>(0.0, 0.0);
		}
		complex<double> *F_translated;
		F_translated = new complex<double> [2*pMax];							// Scattered coeffs solution
		// 7.1 - Calculate abs cross section Cabs --------------------------------------------------------------------------------------------------
		time(&itime);                                       					// calculate initial time
		for(j=0; j<Np; j++)														// for each particle
		{
			// get Beta translation of scattered fields
			// A - Beta terms ---------------------------------------------------------------
			R_rel = particles[j].vR - Origin;			// Mackowski - translate from origin to particle
			TC_beta.init(R_rel, waveK_0, BHreg_inc, nMax);
			TC_beta.populate();
			// ------------------------------------------------------------------------------
			// B - map beta terms obtained from Coupling.cpp to account for efficient multiplication
			for(ii=0; ii<pMax; ii++){
				for(jj=0; jj<pMax; jj++){
					// beta terms -----------------------------------------------------------
					// - Q11
					TC_betaa_t[ii][jj]					= TC_beta.dataApq[jj][ii];
					// - Q12
					TC_betaa_t[ii][jj+pMax]				= TC_beta.dataBpq[jj][ii];
					// - Q21
					TC_betaa_t[ii+pMax][jj]				= TC_beta.dataBpq[jj][ii];
					// - Q22
					TC_betaa_t[ii+pMax][jj+pMax]		= TC_beta.dataApq[jj][ii];
					// ----------------------------------------------------------------------
				}
			}
			// ------------------------------------------------------------------------------
			// C - translate scattering coeffs
			complex<double> *F_local = new complex<double>[2*pMax];
			// translate scattered coeffs from particles to origion
			for(i=0; i<pMax; i++)				// Mackowski - translate excitation coeffs from origin to particle
			{
				F_local[i] 		= excitation.dataIncAp[i];
				F_local[i+pMax] = excitation.dataIncBp[i];
			}
			Algebra::multiplyVectorMatrix(TC_betaa_t, 2*pMax, 2*pMax, F_local, F_translated, consC1, consC0);
			delete [] F_local;

			// D - calculate Cext - Stout 2001 (62)
			for(p=0; p<pMax; p++)
			{
				i=p;
				// Stout
				// Mackowski
				if(j==0){	// particle 1 Cext
					Cext1 += real(	conj(F_translated[i+0*pMax]) * F[j*2*pMax+i+0*pMax]
						  +  conj(F_translated[i+1*pMax]) * F[j*2*pMax+i+1*pMax]);
				}
				if(j==1){	// particle 2 Cext
					Cext2 += real(	conj(F_translated[i+0*pMax]) * F[j*2*pMax+i+0*pMax]
						  +  conj(F_translated[i+1*pMax]) * F[j*2*pMax+i+1*pMax]);
				}
				// total Cext
				Cext += real(	conj(F_translated[i+0*pMax]) * F[j*2*pMax+i+0*pMax]
							 +  conj(F_translated[i+1*pMax]) * F[j*2*pMax+i+1*pMax]);

				// op individual Cext contributions -----------------------------------------------------------------------------------------------------
				// all contributions
				Cext_E[0]	+= - real(conj(F_translated[i+0*pMax]) * F[j*2*pMax+i+0*pMax]) / (waveK_0 * waveK_0)/radius/radius/1e-9/1e-9;
				Cext_H[0]	+= - real(conj(F_translated[i+1*pMax]) * F[j*2*pMax+i+1*pMax]) / (waveK_0 * waveK_0)/radius/radius/1e-9/1e-9;
				// op individual Cext contributions up to Cext_nMax
				for(cext_opi=1; cext_opi<=Cext_nMax; cext_opi++){
					if(p.first==cext_opi){
						Cext_E[cext_opi]	+= - real(conj(F_translated[i+0*pMax]) * F[j*2*pMax+i+0*pMax]) / (waveK_0 * waveK_0)/radius/radius/1e-9/1e-9;
						Cext_H[cext_opi]	+= - real(conj(F_translated[i+1*pMax]) * F[j*2*pMax+i+1*pMax]) / (waveK_0 * waveK_0)/radius/radius/1e-9/1e-9;
					}
				}
			}
		}
		Cext *= -1. / (waveK_0 * waveK_0);
		Cext1 *= -1. / (waveK_0 * waveK_0);
		Cext2 *= -1. / (waveK_0 * waveK_0);
		// -----------------------------------------------------------------------------------------------------------------------------------------
		// op time --------------------
		time(&ftime);                                       // calculate final time
		drift = difftime(ftime, itime);
		hours=int(drift/3600);
		mins =int(60*((drift/3600)-hours));
		secs =int(60*(60*((drift/3600)-hours)-mins));
		if (l==0) cout<<"4 - Cext TE/TM individual (n,m) execution time\t"<<hours<<" hrs "<<mins<<" mins "<<secs<<" secs" << endl;
		//-------------------------------------------------------------------------------------------

		// 7.2 - Calculate abs cross section Cabs --------------------------------------------------------------------------------------------------
		double temp1(0.), temp2(0.);
		for(j=0; j<Np; j++)
		{

			for(i=0; i<2.*pMax; i++){
				Cabs_aux[i]=complex<double>(0., 0.);
			}
			getCabsAux(omega, particles[j], nMax, Cabs_aux);
			// Stout 2001 (68)
			for(i=0; i<pMax; i++)
			{
				temp1 = abs(F[j*2*pMax+i+0*pMax]);
				temp1*=temp1;
				temp2 = abs(F[j*2*pMax+i+1*pMax]);
				temp2*=temp2;
				Cabs += 	temp1 * Cabs_aux[i+0*pMax]
						+ 	temp2 * Cabs_aux[i+1*pMax];
			}
		}
		Cabs *= 1. / (waveK_0 * waveK_0);
		// -----------------------------------------------------------------------------------------------------------------------------------------


		Csca_validate = Cext - Cabs;

		// normalise Cext, Cabs and Csca
		op_CextCabs<<lam_0<<" "<<Cext1.real()<<" "<<Cext2.real()<<" "<<Cext.real()<<" "<<Cabs.real();
		if (l!=Cext_steps-1) op_CextCabs<<endl;
		// op all Cext TE/TM contributions
		for(cext_opi=0; cext_opi<=Cext_nMax; cext_opi++){
			op_CextCabs_EH<<real(Cext_E[cext_opi])<<" "<<real(Cext_H[cext_opi])<<" ";
		}
		op_CextCabs_EH<<endl;
		// op individual Cext contributions
		op_CextCabs_EH_1<<lam_0*1e9<<" "<<real(Cext_E[1])<<" "<<real(Cext_H[1])<<endl;
		op_CextCabs_EH_2<<lam_0*1e9<<" "<<real(Cext_E[2])<<" "<<real(Cext_H[2])<<endl;
		op_CextCabs_EH_3<<lam_0*1e9<<" "<<real(Cext_E[3])<<" "<<real(Cext_H[3])<<endl;
		op_CextCabs_EH_4<<lam_0*1e9<<" "<<real(Cext_E[4])<<" "<<real(Cext_H[4])<<endl;
		op_CextCabs_EH_5<<lam_0*1e9<<" "<<real(Cext_E[5])<<" "<<real(Cext_H[5])<<endl;

		// at each scan delete all matrices -------------------------------------------------------
		delete[] B;
		delete[] E;
		delete[] B0;
		delete[] Cabs_aux;
		delete[] F_translated;
		for(int i=0; i<2*pMax; i++){
			delete[] T_matrix[i];
			delete[] TC_betaa_t[i];		delete[] TC_alpha_t[i];		delete[] mT_TC_alpha[i];
		}
			delete[] T_matrix;
			delete[] TC_betaa_t;		delete[] TC_alpha_t;		delete[] mT_TC_alpha;

		for(int i=0; i<Np*2*pMax; i++){
			delete[] SS[i];			delete[] TT[i];
		}
			delete[] SS;			delete[] TT;

	}
	// finish scan here -------------------------------------------------------------------------

	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
	// pre-plotting ----------------------------------------------------------------------------
	// 5 - obtain internal coeffs for each particle --------------------------------------------
	complex<double> *I_aux, *I_aux_j;
	I_aux = new complex<double> [2*Np*pMax];
	I_aux_j = new complex<double> [2*pMax];
	for (j=0; j<Np; j++){											// loop over all particles
		// Aux coeffs
		getIaux (omega, particles[j], nMax, I_aux_j);
		// internal coeffs
		for(i=0; i<p.max(nMax); i++){
			// particle 1
			I_aux[j*2*pMax+i+0*pMax] = I_aux_j[i+0*pMax] * F[j*2*pMax+i+0*pMax];
			I_aux[j*2*pMax+i+1*pMax] = I_aux_j[i+1*pMax] * F[j*2*pMax+i+1*pMax];
		}
	}
	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

	//	ofstream T1_op_real("T1_op_real");
	//	for(i=0; i<Np*2*pMax; i++){
	//		T1_op_real<<i<<"\t"<<real(BB[i])<<endl;
	//	}

	// produce scattered fields solution --------------------------------------------------------------------------------------------
	// do a 2D scan of EF profile ---------------------------------------------------------------------------------------------------
	double xmin=-4000.*1e-9, xmax=4000.*1e-9;
	double zmin=-4000.*1e-9, zmax=4000.*1e-9;
	double dx = (xmax-xmin)/steps;
	double dz = (zmax-zmin)/steps;
	Spherical<double> R_op, R_j;
	int Np_in=0;
	SphericalP<complex<double> > sum_up, sum_par, sum_sca, sum_inc, sum_ext, sum_int, sum_tot;

	// Cartesian 2D plot ----------------------------------------------------------------------------------------------------------------
	time(&itime);                                       // calculate initial time
	// start 2D scan
	// op inc
	ofstream E1_inc_real("E1_inc_real.txt");		ofstream E1_inc_imag("E1_inc_imag.txt");
	ofstream E2_inc_real("E2_inc_real.txt");		ofstream E2_inc_imag("E2_inc_imag.txt");
	ofstream E3_inc_real("E3_inc_real.txt");		ofstream E3_inc_imag("E3_inc_imag.txt");
	//op ext
	ofstream E1_real("E1_real.txt");		ofstream E1_imag("E1_imag.txt");
	ofstream E2_real("E2_real.txt");		ofstream E2_imag("E2_imag.txt");
	ofstream E3_real("E3_real.txt");		ofstream E3_imag("E3_imag.txt");

	// xy plane
	z=0.*1e-9;
	for(ii=0, x=xmin; ii<=steps; ii++, x+=dx){
		for(jj=0, y=zmin; jj<=steps; jj++, y+=dx){

	// xz plane

			// calculate incident, scattered, external & total fields -------------------------------------------------------------------
			sum_inc=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
			sum_sca=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
			sum_ext=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
			sum_tot=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));

			// Cartesian observation point ----------------------------------------------------------------------------------------------
			R_op=Tools::toSpherical(Cartesian<double>(x+1e-12, y+1e-12, z+1e-12));

			// 0 - does location of 'R_op' fall within any particle ---------------------------------------------------------------------
			Np_in=checkInner(R_op, particles, Np);

			// 1 - internal field -------------------------------------------------------------------------------------------------------
			// if inside any particle - contribution form one particle only (Np_in) -----------------------------------------------------
			sum_int=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
			if(Np_in!=-1){
				j=Np_in;														// assign j to Np_in
				// scattering filed Aux functions
				waveK_j = omega *sqrt(particles[j].elmag.epsilon*particles[j].elmag.epsilon_r * particles[j].elmag.mu*particles[j].elmag.mu_r);	// particle wavenumber
				R_j = Tools::toPoint(R_op, particles[j].vR);					// local observation point
				AuxCoefficients AuxCoeff_int;
				AuxCoeff_int.init(R_j, waveK_j, BHreg_inc, nMax);
				AuxCoeff_int.populate();
				// Calculate internal field ---------------------------------------------------------------------------------------------
				sum_int=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
				for(i=0; i<p.max(nMax); i++){
					p=i;
					if(p.first>=n_op_start && p.first<=n_op_final){
					// internal from scattered coeffs 	- particle j
					sum_int.rrr += I_aux[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].rrr + I_aux[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].rrr;
					sum_int.the += I_aux[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].the + I_aux[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].the;
					sum_int.phi += I_aux[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].phi + I_aux[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].phi;
					// only to check plane wave is uniformally extended to both particles
					}
				}
			}
			// if inside sphere ---------------------------------------------------------------------------------------------------------

			// 2 - external field -------------------------------------------------------------------------------------------------------
			// if outside all spheres ---------------------------------------------------------------------------------------------------
			else{
				// 2.1 scattered field from all particles -------------------------------------------------------------------------------
				for (j=0; j<Np; j++){											// scattered field contribution from all particles
					// scattering filed Aux functions
					R_j = Tools::toPoint(R_op, particles[j].vR);				// local observation point
					AuxCoefficients AuxCoeff_sca;
					AuxCoeff_sca.init(R_j, waveK_0, BHreg_sca, nMax);
					AuxCoeff_sca.populate();
					// Calculate scattered field ----------------------------------------------------------------------------------------
					sum_par=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
					for(i=0; i<p.max(nMax); i++){
						p=i;
						if(p.first>=n_op_start && p.first<=n_op_final){
						// sca from sphere 1
						sum_par.rrr += F[j*2*pMax+i+0*pMax]*AuxCoeff_sca.dataMp[i].rrr + F[j*2*pMax+i+1*pMax]*AuxCoeff_sca.dataNp[i].rrr;
						sum_par.the += F[j*2*pMax+i+0*pMax]*AuxCoeff_sca.dataMp[i].the + F[j*2*pMax+i+1*pMax]*AuxCoeff_sca.dataNp[i].the;
						sum_par.phi += F[j*2*pMax+i+0*pMax]*AuxCoeff_sca.dataMp[i].phi + F[j*2*pMax+i+1*pMax]*AuxCoeff_sca.dataNp[i].phi;
						}
					}
					// get total scattered - accumulation from all particles
					sum_sca = sum_sca + sum_par;							// total scattered field includes scattered scattered from all
				}	// for loop (j)
				// ----------------------------------------------------------------------------------------------------------------------

				// 2.2 global incident field --------------------------------------------------------------------------------------------
				AuxCoefficients AuxCoeff_inc;
				AuxCoeff_inc.init(R_op, waveK_0, BHreg_inc, nMax);
				AuxCoeff_inc.populate();
				// Calculate incident field ---------------------------------------------------------------------------------------------
				sum_inc=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
				for(i=0; i<p.max(nMax); i++){
					p=i;
					if(p.first>=n_op_start && p.first<=n_op_final){
	//				// local to particle
	//				sum_inc.rrr		+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].rrr + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].rrr;
	//				sum_inc.the 	+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].the + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].the;
	//				sum_inc.phi 	+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].phi + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].phi;
					// local to system origin
					sum_inc.rrr		+= 	Inc[i+0*pMax]*AuxCoeff_inc.dataMp[i].rrr + Inc[i+1*pMax]*AuxCoeff_inc.dataNp[i].rrr;
					sum_inc.the 	+= 	Inc[i+0*pMax]*AuxCoeff_inc.dataMp[i].the + Inc[i+1*pMax]*AuxCoeff_inc.dataNp[i].the;
					sum_inc.phi 	+= 	Inc[i+0*pMax]*AuxCoeff_inc.dataMp[i].phi + Inc[i+1*pMax]*AuxCoeff_inc.dataNp[i].phi;
					}
				}
				// ----------------------------------------------------------------------------------------------------------------------

				// 2.3 total external field = global incident field + scattered field from all particles --------------------------------
				sum_ext = sum_inc + sum_sca;								// total external field includes total scattered + incident

			}	// else if(!Np_in)
			// if outside all spheres ---------------------------------------------------------------------------------------------------

			// 3 - total = external (or) internal field ---------------------------------------------------------------------------------
			sum_tot = sum_ext + sum_int;						// total field will include either internal for a given particle, or external from all

			// inc
			E1_inc_real<<real(sum_inc.rrr)<<" ";		E1_inc_imag<<imag(sum_inc.rrr)<<" ";
			E2_inc_real<<real(sum_inc.the)<<" ";		E2_inc_imag<<imag(sum_inc.the)<<" ";
			E3_inc_real<<real(sum_inc.phi)<<" ";		E3_inc_imag<<imag(sum_inc.phi)<<" ";
			// tot
			E1_real<<real(sum_tot.rrr)<<" ";		E1_imag<<imag(sum_tot.rrr)<<" ";
			E2_real<<real(sum_tot.the)<<" ";		E2_imag<<imag(sum_tot.the)<<" ";
			E3_real<<real(sum_tot.phi)<<" ";		E3_imag<<imag(sum_tot.phi)<<" ";
			// tot mag
		}
		// inc
		E1_inc_real<<endl;		E1_inc_imag<<endl;
		E2_inc_real<<endl;		E2_inc_imag<<endl;
		E3_inc_real<<endl;		E3_inc_imag<<endl;
		// tot
		E1_real<<endl;		E1_imag<<endl;		//E1_mag<<endl;
		E2_real<<endl;		E2_imag<<endl;		//E2_mag<<endl;
		E3_real<<endl;		E3_imag<<endl;		//E3_mag<<endl;
	}
	// op time --------------------
	time(&ftime);                                       // calculate final time
	drift = difftime(ftime, itime);
	hours=int(drift/3600);
	mins =int(60*((drift/3600)-hours));
	secs =int(60*(60*((drift/3600)-hours)-mins));
	cout<<"5 - Fields profile execution time             \t"<<hours<<" hrs "<<mins<<" mins "<<secs<<" secs" << endl;
	// ---------------------------------------------------------------------------------------------------------------------------------- //


	cout<<"done";
	return 0;
}


