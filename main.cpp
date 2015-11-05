//OPTIMET 3D - main.cpp
//
//Restricted to testing for the moment.
//

#include "mpi.h"
//#include <iostream>
#include <fstream>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <time.h>
#include <complex>
#include "Spherical.h"
#include "AuxCoefficients.h"
#include "Bessel.h"
//#include "Reader.h"
#include "CompoundIterator.h"
#include "Legendre.h"
#include "Excitation.h"
#include "Coupling.h"
#include "SphericalP.h"
//#include "Result.h"
#include "Symbol.h"
#include "Tools.h"

#include "AJ_AuxFuns.h"
#include "PeriodicCoupling.h"
//#include "Simulation.h"
//#include "Solver.h"
#include "AlgebraS.h"
//#include "AlgebraP.h"

using namespace std;

void Pol2Cart_main(double rho_, double theta_, double &x_, double &y_)
{
	x_ = rho_ * cos(theta_);
	y_ = rho_ * sin(theta_);
}

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
//	Jrho_n.init(rho*r_0, 0, 0, nMax_);
////	Jrho_n.init(1., 0, 0, nMax_);
//	if(Jrho_n.populate())
//	{
//		cerr << "Error computing Bessel functions. Amos said: " << Jrho_n.ierr << "!";
//		return 1;
//	}
//
//	cout<<"0 - So far so good"<<endl;
//
//	H_n.init(r_0, 1, 0, nMax_);
////	H_n.init(1., 1, 0, nMax_);
//	if(H_n.populate())
//	{
//		cerr << "Error computing Hankel functions. Amos said: " << H_n.ierr << "!";
//	}
//
//	cout<<"1 - So far so good"<<endl;
//
//	Hrho_n.init(rho*r_0, 1, 0, nMax_);
////	Hrho_n.init(1., 1, 0, nMax_);
//	if(Hrho_n.populate())
//	{
//		cerr << "Error computing Hankel functions. Amos said: " << Hrho_n.ierr << "!";
//	}
//
//	cout<<"2 - So far so good"<<endl;

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
//				psi = r_0*J_n.data[p.first];
//				dpsi= r_0*J_n.ddata[p.first] + J_n.data[p.first];
//
//				ksi = r_0*H_n.data[p.first];
//				dksi= r_0*H_n.ddata[p.first] + H_n.data[p.first];
//
//				psirho = r_0*rho*Jrho_n.data[p.first];
//				dpsirho= r_0*rho*Jrho_n.ddata[p.first] + Jrho_n.data[p.first];
//
//				ksirho = r_0*rho*Hrho_n.data[p.first];
//				dksirho= r_0*rho*Hrho_n.ddata[p.first] + Hrho_n.data[p.first];
//
//				//TE Part
//				T_local_[p][q] = (psi/ksi) 		* (mu_sob*dpsi/psi - rho*dpsirho/psirho)
//												/ (rho*dpsirho/psirho - mu_sob*dksi/ksi);
//				//TM part
//				T_local_[(int)p+pMax][(int)q+qMax] = (psi/ksi)	* (mu_sob*dpsirho/psirho - rho*dpsi/psi)
//																/ (rho*dksi/ksi - mu_sob*dpsirho/psirho);

				//TE Part
				T_local_[p][q] 						= 1.;
				//TM part
				T_local_[(int)p+pMax][(int)q+qMax] 	= 1.;


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

//			// AJ - from excitation
//			ksi = r_0*H_n.data[p.first];
//			dksi= r_0*H_n.ddata[p.first] + H_n.data[p.first];
//			//TE Part
//			I_aux_[p] 			= 	(mu_j*rho*ksi*dpsi - mu_j*rho*psi*dksi) /
//									(mu_0*rho*ksi*dpsirho - mu_j*psirho*dksi);
//			//TM part
//			I_aux_[(int)p+pMax]	= 	(mu_j*rho*dksi*psi - mu_j*rho*dpsi*ksi) /
//									(mu_0*rho*dksi*psirho - mu_j*dpsirho*ksi);

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

//int getCabsAux (double omega_, Scatterer particle, int nMax_, complex<double> *Cabs_aux_){
//
//	complex<double> temp1(0., 0.), temp2(0., 0.);
//	complex<double> k_s = omega_ * sqrt(particle.elmag.epsilon * particle.elmag.mu);
//	complex<double> k_b = omega_ * sqrt(consEpsilon0 * consMu0);
//
//	complex<double> rho = k_s / k_b;
//	complex<double> r_0 = k_b * particle.radius;
//	complex<double> mu_j  = particle.elmag.mu;
//	complex<double> mu_0  = consMu0;
//	complex<double> mu_sob = particle.elmag.mu /consMu0;
//
//	complex<double> psi(0., 0.), ksi(0., 0.);
//	complex<double> dpsi(0., 0.), dksi(0., 0.);
//
//	complex<double> psirho(0., 0.);
//	complex<double> dpsirho(0., 0.);
//
//	Bessel J_n;
//	Bessel Jrho_n;
//
//	J_n.init(r_0, 0, 0, nMax_);
//	if(J_n.populate())
//	{
//		cerr << "Error computing Bessel functions. Amos said: " << J_n.ierr << "!";
//		return 1;
//	}
//
//	Jrho_n.init(rho*r_0, 0, 0, nMax_);
//	if(Jrho_n.populate())
//	{
//		cerr << "Error computing Bessel functions. Amos said: " << Jrho_n.ierr << "!";
//		return 1;
//	}
//
//	// AJ - direct
//	Bessel H_n;
//	H_n.init(r_0, 1, 0, nMax_);
//	if(H_n.populate())
//	{
//		cerr << "Error computing Hankel functions. Amos said: " << H_n.ierr << "!";
//	}
//
//	CompoundIterator p;
//
//	int pMax = p.max(nMax_);
//
//	/**
//	 * @todo Modify this for non-spherical objects (need matrix, also see Solver::solveInternal).
//	 */
//
//	for(p=0; (int)p<pMax; p++)
//	{
//			// obtain Riccati-Bessel functions
//			psi = r_0*J_n.data[p.first];
//			dpsi= r_0*J_n.ddata[p.first] + J_n.data[p.first];
//
//			psirho = r_0*rho*Jrho_n.data[p.first];
//			dpsirho= r_0*rho*Jrho_n.ddata[p.first] + Jrho_n.data[p.first];
//
//			// Stout 2002 - from scattered
//			//TE Part
//			temp1 = complex<double>(0., 1.)*rho*mu_0*conj(mu_j)*conj(psirho)*dpsirho;
//			temp2 = abs((mu_j*psirho*dpsi - mu_0*rho*dpsirho*psi));
//			temp2*=temp2;
//			Cabs_aux_[p] 			= real(temp1)/temp2;
//
//			//TM part
//			temp1 = complex<double>(0., 1.)*conj(rho)*mu_0*mu_j*conj(psirho)*dpsirho;
//			temp2 = abs((mu_0*rho*psirho*dpsi - mu_j*dpsirho*psi));
//			temp2*=temp2;
//			Cabs_aux_[(int)p+pMax]	= real(temp1)/temp2;
//	}
//	return 0;
//}

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

//	// from xml ----------------------------------------------------------------------------
//	Simulation simulation;
//	simulation.init(argv[1]);
//	simulation.run();
//	simulation.done();

	//************  MPI ***************************//
	int mpirank, mpisize;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
	MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
	//************  MPI ***************************//

	// Ahmed -------------------------------------------------------------------------------
	int nMax=20;							// max no of harmonics
	int Obs_part=0;							// cont check - choose observation particle
	int steps = 10;						// op steps for 2D plot
	int the_range = 10;					// cont check steps
	int phi_range = 10;					// cont check steps
	int Cext_nMax=1;						// individual TE/TM components to calculate

	CompoundIterator p, q;
	int pMax = p.max(nMax);

	// parallel partitioning
	int Np=4;
	double dblock=double((double(Np*Np))/512.);
	int N, M, Nb, Mb;
	N	= Np*2*pMax;					// no of rows in global matrix
	M	= Np*2*pMax;					// no of cols in global matrix
//	Nb	= int(dblock*double(2*pMax));	// no of rows in primitive block matrix - to be changed to 2*pMax
	//Mb	= int(dblock*double(2*pMax));	// no of cols in primitive block matrix - to be changed to 2*pMax
	Nb	= 2*2*pMax;						// no of rows in primitive block matrix - to be changed to 2*pMax
	Mb	= 2*2*pMax;						// no of cols in primitive block matrix - to be changed to 2*pMax
//	int procrows = 2;					// initialse Pr=procrows, Pc=proccols
	//int proccols = 2;					// initialse Pr=procrows, Pc=proccols
	int procrows = Np/2;
	int proccols = Np/2;								// initialse Pr=procrows, Pc=proccols

	if(mpirank==0){
		cout<<" Output partitioning parameters -------- "<<endl;
		cout<<" Number of particles	= "<<Np<<endl;
		cout<<" dblock factor		= "<<dblock<<endl;
		cout<<" Matrix size			= "<<Np*2*pMax<<endl;
		cout<<" I-matrix size		= "<<2*pMax<<endl;
		cout<<" block size rows Mb	= "<<Mb<<endl;
		cout<<" block size cols Mb	= "<<Nb<<endl;
		cout<<" Process Grid rows	= "<<procrows<<endl;
		cout<<" Process Grid cols	= "<<proccols<<endl;
		cout<<endl;
	}

	// scan parameters ---------------------------------------------------------------------
	double inc =1;
	double inc_interval=1000;
//	double lam_start=1300.+(inc-1.)*inc_interval;	// Cext scan start at
//	double lam_final=1300.+(inc)*inc_interval;		// Cext scan finish at
	double lam_start=1300.;							// Cext scan start at
	double lam_final=1300.;							// Cext scan finish at
	int Cext_steps=lam_final-lam_start+1;						// Cext steps - user defined
//	Cext_steps = 101;
	double lam_steps=((lam_final-lam_start)/double(Cext_steps-1));	// Cext scan steps
	// -------------------------------------------------------------------------------------

	// Translation Coefficients test - from Periodic_Coupling.cpp + Coupling.cpp
	int l, i, j, ii, jj;
	complex<double> c_temp(0., 0.);
//	CompoundIterator p, q;
//	int pMax = p.max(nMax);
	int BHreg_sca=0;					// compute regular values of BH - non-regular
	int BHreg_inc=1;					// compute regular values of BH - regular

	// input variables
	double lam_0=lam_start*1e-9;				// input wavelength
	double lam_ip=lam_0;				// default scan
	complex<double> eps_i=complex<double>(13, 0.);
	double freq = consC/lam_0;
	double omega=2*consPi*freq;
	complex<double> waveK_0 = omega*sqrt(consEpsilon0 * consMu0);	// background wave-number
	complex<double> waveK_j;
	double x(0.), y(0.), z(0.);

	if(mpirank==0) cout<<"AJ intake on OPTIMET-3D -------------------- "<<endl;
	if(mpirank==0) cout<<"Scan Cext_setps = "<<Cext_steps<<endl;

	// AJ - Chiral Auto build --------------------------------------------------
	//	Spiral structure building block
	double radius = 500.0;
	double sep =100.0;
	int arms =3;							// Number of arms in spiral structure
	int No = 2; 							// Number of particles on each arm including center
	//Distance based simulation
	double Theta = consPi/(No-1); 			// Calculate the separation angle
	double d = 2.0*radius + sep;			// Separation between two particles
	double Phi = (consPi-Theta)/2.0;		// remaining angles on equal sided triangle
	double R = d * sin(Phi)/sin(Theta);		// Determine R - given by triangle rule : d/sin(Theta) = R/sin(Phi)
	if(No==2)								// Special case where Theta==consPi
		R=d;

//	int Np =(No-1)*arms + 1;				// total number of particles in system
	//Np=4;
	if(mpirank==0) cout<<"Total number of particles = "<<Np<<endl;
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
//		Tools::Pol2Cart(R, i*Theta, x_, y_);
		Pol2Cart_main(R, i*Theta, x_, y_);
		x_arm[j] = x_ + R;
		y_arm[j] = y_;
		j++;
	}

	//Create vectors for r, theta, x and y, X and Y
	double X_= 0.;							// to store x-axis location of particles for all arm
	double Y_= 0.;							// to store y-axis location of particles for all arm
	eps_i=complex<double>(13, 0.);			// fixed here - Sellmier is used below
	// particle 1
	particles[0].vR = Spherical<double>(Tools::toSpherical(Cartesian<double>(0., 0., 0.*1e-9)));	// The coordinates of the center of the scatterer
	particles[0].elmag.epsilon=consEpsilon0*eps_i;						// The electromagnetic properties of the scatterer
	particles[0].elmag.mu=consMu0;										// The electromagnetic properties of the scatterer.
	particles[0].nMax=nMax;												// Maximum value of the n iterator.
	particles[0].radius=radius*1e-9;									// The radius of a sphere encompassing the scatterer.
	
	// chain of particles ----------------------------------------------------------------------
	sep =100;
	double radiusp = radius + sep;
	for(i=0; i<Np; i++){
		particles[i].vR = Spherical<double>(Tools::toSpherical(Cartesian<double>(0., 0., (i*2.0*radiusp)*1e-9 )));	// The coordinates of the center of the scatterer
		particles[i].elmag.epsilon=consEpsilon0*eps_i;						// The electromagnetic properties of the scatterer
		particles[i].elmag.mu=consMu0;										// The electromagnetic properties of the scatterer.
		particles[i].nMax=nMax;												// Maximum value of the n iterator.
		particles[i].radius=radius*1e-9;									// The radius of a sphere encompassing the scatterer.
	}
	// chain of particles ----------------------------------------------------------------------
	j=1;
//	for(int arm_i=0; arm_i<arms; arm_i++){
//		// Translate points around origin to finally define locations of all particles on arms
//		for(i=0; i<No-1; i++){
//			// obtain global coordinates
//			X_ = x_arm[i]*cos(double(arm_i)*Theta_rot) - y_arm[i]*sin(double(arm_i)*Theta_rot);
//			Y_ = x_arm[i]*sin(double(arm_i)*Theta_rot) + y_arm[i]*cos(double(arm_i)*Theta_rot);
//
//			// populate each particle parameters
//			eps_i=complex<double>(13, 0.);	// fixed here - Sellmier is used below
//			// x normal
////			particles[j].vR = Spherical<double>(Tools::toSpherical(Cartesian<double>(0., X_*1e-9, Y_*1e-9)));	// The coordinates of the center of the scatterer
//			// y normal
////			particles[j].vR = Spherical<double>(Tools::toSpherical(Cartesian<double>(X_*1e-9, 0., Y_*1e-9)));	// The coordinates of the center of the scatterer
//			// z normal
//			particles[j].vR = Spherical<double>(Tools::toSpherical(Cartesian<double>(X_*1e-9, Y_*1e-9, 0.)));	// The coordinates of the center of the scatterer
//			particles[j].elmag.epsilon=consEpsilon0*eps_i;						// The electromagnetic properties of the scatterer
//			particles[j].elmag.mu=consMu0;										// The electromagnetic properties of the scatterer.
//			particles[j].nMax=nMax;												// Maximum value of the n iterator.
//			particles[j].radius=radius*1e-9;									// The radius of a sphere encompassing the scatterer.
//			j++;
//		}
//	}
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
//	vKInc.the = consPi*(90./180.);
//	vKInc.phi = consPi*(90./180.);
	// x-axis propagation
//	vKInc.the = consPi*(90./180.);
//	vKInc.phi = consPi*(0.0001/180.);

	// Set Eaux Polarisation -----------------------------------------------------------------------------------------------------
	Spherical<double> vAux = Spherical<double>(0.0, vKInc.the, vKInc.phi);
	SphericalP<complex<double> > Eaux;
	// linear
//	Eaux = SphericalP<complex<double> >(complex<double>(0.0, 0.0),
//										complex<double>(0.0, 0.0),
//										complex<double>(1.0, 0.0));
	// circular
	double oSqrt2 = 1./sqrt(2.);
	Eaux = SphericalP<complex<double> >(complex<double>(0.0, 0.0),
										complex<double>(oSqrt2, 0.0),
										complex<double>(0.0, -oSqrt2));
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
	l=0;
	for(l=0, lam_ip=lam_start; l<Cext_steps; l++, lam_ip+=lam_steps){
	//	for(l=0, lam_ip=lam_start; l<1; l++, lam_ip+=lam_steps){

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
//		// System of equations vectors ------------------------------------------------------------
		complex<double> *B, *E;// *F;
		B = new complex<double> [Np*2*pMax];							// translated Incident coeffs
		E = new complex<double> [Np*2*pMax];							// Excitation coeffs solution
//		F = new complex<double> [Np*2*pMax];							// Scattered coeffs solution
		// -----------------------------------------------------------------------------------------
		// System of equations global matrices -----------------------------------------------------
		complex<double> **SS;// **TT;
		SS = Tools::Get_2D_c_double(Np*2*pMax, Np*2*pMax);
//		TT = Tools::Get_2D_c_double(Np*2*pMax, Np*2*pMax);
		// -----------------------------------------------------------------------------------------

		// start scan here - Cext, Cabs ------------------------------------------------------------
		complex<double> Cext(0.0, 0.0);
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

		if(mpirank==0) cout<<"compute lam = "<<lam_0<<endl;

		// Sellmier dispersion relationship
		for (j=0; j<Np; j++){											// loop over all particles
			double lam_sl = lam_0*1e6;
			eps_i = ( 1 + 10.6684293*pow(lam_sl,2)/(pow(lam_sl,2)-pow(0.301516485,2)) + 0.003043475*pow(lam_sl,2)/(pow(lam_sl,2)-pow(1.13475115,2)) + 1.54133408*pow(lam_sl,2)/(pow(lam_sl,2)-pow(1104.0,2)) );
//			eps_i = 12;
			particles[j].elmag.epsilon=consEpsilon0*eps_i;				// The electromagnetic properties of the scatterer
		}

		if (l==0 && mpirank==0) cout<<"1 - Build system of equations 	: E = [SS].B"<<endl;
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
//							TT[ii+i*2*pMax][jj+j*2*pMax] 		= T_matrix[ii][jj];
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
//							TT[ii+i*2*pMax][jj+j*2*pMax] 		= complex<double>(0., 0.);
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
		if (l==0 && mpirank==0) cout<<"1 - Build system of equations execution time\t"<<hours<<" hrs "<<mins<<" mins "<<secs<<" secs" << endl;
		//-------------------------------------------------------------------------------------------

		// ------------------------------------------------------------------------------------------
		if (l==0 && mpirank==0) cout<<"2 - Solve system of equations 	: E = [SS]^-1.B"<<endl;
		// 4.4 excitation coeffs solution - solve E = [SS].B ----------------------------------------
		time(&itime);                                       // calculate initial time
		// soleve in serial -------------------------------------------------------------
//		AlgebraS::solveMatrixVector(		SS, 		Np*2*pMax, 		Np*2*pMax, B, E);
		// ------------------------------------------------------------------------------

		// solve in parallel ------------------------------------------------------------
		MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
		MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
		bool mpiroot = (mpirank == 0);
		int iZERO = 0;			// used in function numroc_ (number of local rows or columns)
		// Global matrix A=[16x16]
		// A - We assign each subblock to contain 8x8 elements
//		int N, M, Nb, Mb;
//		N	= Np*2*pMax;		// no of rows in global matrix
//		M	= Np*2*pMax;		// no of cols in global matrix
//		Nb	= (Np/2)*2*pMax;				// no of rows in primitive block matrix - to be changed to 2*pMax
//		Mb	= (Np/2)*2*pMax;				// no of cols in primitive block matrix - to be changed to 2*pMax
		// map global 2D matrix into 1D matrix with transposition --------------
		complex<double> *A_glob; //*A_glob2;
		A_glob		= new complex<double>[N*M];
//		A_glob2		= new complex<double>[N*M];
		// map matrix (with transposition!)
		ii=0;
		// column major storage - i.e. transpose matrix to be solved
		for(j=0; j<M; j++){				// columns
			for(i=0; i<N; i++){			// rows
				A_glob[ii] = SS[i][j];
				ii++;
			}
		}
/*
		if(mpirank==0){
			cout << "Matrix A_glob:\n";
			for (int r = 0; r < N; ++r) {
				for (int c = 0; c < M; ++c) {
					cout <<" "<< *(A_glob + N*c + r) << " ";
				}
				cout << "\n";
			}
			cout << endl;
		}
*/
		// -------------------------------------------------------------------

		// B - Begin Cblas context -------------------------------------------
		// We assume that we have 4 processes and place them in a 2-by-2 grid
		int ctxt, myid, myrow, mycol, numproc;
//		int procrows = 2, proccols = 2;								// initialse Pr=procrows, Pc=proccols
//		int procrows = Np, proccols = Np;							// per bi-particle partitioning
		Cblacs_pinfo(&myid, &numproc);								// obtain myid=uniquely identifies each process, numproc=the number of processes available for BLACS use 
		Cblacs_get(0, 0, &ctxt);									// obtain a default system context
		Cblacs_gridinit(&ctxt, "Row-major", procrows, proccols);	// mapping the processes to the process grid
		Cblacs_pcoord(ctxt, myid, &myrow, &mycol);					// myid=The process number who's coordinates are to be determined
/*		
		if(mpirank==1){
			cout << "Check - Cblacs started:" << endl;
			cout<<"myid		= "<<myid<<endl;
			cout<<"numproc	= "<<numproc<<endl;
			cout<<"ctxt		= "<<ctxt<<endl;
			cout<<"myrow	= "<<myrow<<endl;
			cout<<"mycol	= "<<mycol<<endl;
		}
		// Print grid pattern 
//		if (myid == 0){
		if(mpirank==1) cout << "Processes grid pattern:" << endl;
		for (int r = 0; r < procrows; ++r) {
			for (int c = 0; c < proccols; ++c) {
				Cblacs_barrier(ctxt, "All");
				if (myrow == r && mycol == c) {
					cout << myid << " " << flush;
				}
			}
			Cblacs_barrier(ctxt, "All");
			if (myid == 0)
				cout << endl;
		}
		Cblacs_barrier(ctxt, "All");
//		}
*/
		// ------------------------------------------------------------------
		
		// C - Reserve space for local matrices -----------------------------
		// Number of rows and cols owned by the current process
		// The numroc utility will compute this for you. Just call it once for the rows and once for the columns. 
		// Then nrows*ncols will be the number of total local entries (i.e. the entries of A_loc).
		int nrows = numroc_(&N, &Nb, &myrow, &iZERO, &procrows);		// LOCr
		int ncols = numroc_(&M, &Mb, &mycol, &iZERO, &proccols);		// LOCc
		for (int id = 0; id < numproc; ++id) {
			Cblacs_barrier(ctxt, "All");
		}
		complex<double> *A_loc, *b_loc;
		b_loc = new complex<double>[nrows];			// corresponds to distributed solution
		A_loc = new complex<double>[nrows*ncols];	// distributed to global matrix A
		for (int iloc = 0; iloc < nrows*ncols; ++iloc) *(A_loc+iloc)=complex<double>(0., 0.);
		for (int iloc = 0; iloc < nrows; ++iloc)		b_loc[iloc]=complex<double>(0., 0.);
		// ------------------------------------------------------------------

		// D - Partitioning the global matrix -------------------------------
		// - Note that in the for loops we don’t just add 1 to the current row or column, 
		// but the corresponding block size, so at the beginning of every iteration A_glob[r, c] 
		// is the first entry in the new block that has to be sent.
		// - Using sendr and sendc we identify the process the block has to be sent to. 
		// - The recvr and recvc help us find where in the local matrix the first entry will be placed.
		// - A fundamental parameter for the matrices is the trailing dimension, 
		// i.e. the dinstance between the first element of two contiguos columns.
		// For the global matrix, the leading dimension is the number of rows, i.e. N; 
		// for the local matrix, this is the number of local-stored rows, i.e. nrows.
		// Partition matrix using send-recieve operations - to be changed to local assignment
		int sendr = 0, sendc = 0, recvr = 0, recvc = 0;
		for (int r = 0; r < N; r += Nb, sendr=(sendr+1)%procrows) {
			sendc = 0;
			// Number of rows to be sent
			// Is this the last row block?
			int nr = Nb;
			if (N-r < Nb)
				nr = N-r;
 
			for (int c = 0; c < M; c += Mb, sendc=(sendc+1)%proccols) {
				// Number of cols to be sent
				// Is this the last col block?
				int nc = Mb;
				if (M-c < Mb)
					nc = M-c;
 
				if (mpiroot) {
					// Send a nr-by-nc submatrix to process (sendr, sendc)
					Czgesd2d(ctxt, nr, nc, A_glob+N*c+r, N, sendr, sendc);
				}
 
				if (myrow == sendr && mycol == sendc) {
					// Receive the same data
					// The leading dimension of the local matrix is nrows!
					Czgerv2d(ctxt, nr, nc, A_loc+nrows*recvc+recvr, nrows, 0, 0);
					recvc = (recvc+nc)%ncols;
				}
 
			}
 
			if (myrow == sendr)
				recvr = (recvr+nr)%nrows;
		}
		// ------------------------------------------------------------------

		// partition RHS vector ---------------------------------------------
		// loop over primitive block rows
		for (int iloc = 0; iloc<Mb; ++iloc){
//			if(myrow==0 && mycol==0)
				b_loc[iloc] = iloc;
//				b_loc[iloc] = B[iloc];
//			if(myrow==0 && mycol==1)
				//b_loc[iloc] = B[iloc+pMax];
		}
		// ------------------------------------------------------------------

		// E - solve --------------------------------------------------------
		// according to AJ's implementaion ----------------------------------
		Cblacs_barrier(ctxt, "All");
		int desca[9], descb[9];
		int ZERO = 0;
		int ONE = 1;
		int IERR;

		int rows_A_ = N;
		int cols_A_ = M;
		int block_rows_A = Nb;
		int block_cols_A = Mb;
		int LLD = nrows;
		// - A general M_ by N_ distributed matrix is defined by its dimensions, 
		// - the size of the elementary MB_ by NB_ block used for its decomposition, 
		// - the coordinates of the process having in its local memory the first matrix entry (RSRC_,CSRC_), 
		// - the BLACS context (CTXT_) in which this matrix is defined. 
		// - Finally, a local leading dimension LLD_ is associated with the   
		// local memory address pointing to the data structure used for the local storage of this distributed matrix.
		descinit_(desca,	&rows_A_,	&cols_A_,	&block_rows_A, &block_cols_A, &ZERO, &ZERO, &ctxt, &LLD, &IERR);
		descinit_(descb,	&rows_A_,	&ONE,		&block_rows_A, &block_cols_A, &ZERO, &ZERO, &ctxt, &LLD, &IERR);
		int *ipiv = new int[nrows+block_rows_A+1];
		pzgesv(&block_rows_A, &ONE, A_loc, &ONE, &ONE, desca, ipiv, b_loc, &ONE, &ONE, descb, &IERR);
		Cblacs_barrier(ctxt, "All");
/*
		for(int sendr = 0; sendr < (int)sqrt(numproc); sendr++)
		{
            if(myrow == sendr)
            {
                Czgesd2d(ctxt, LLD, ONE, b_loc, 1, 0, 0);
//				Czgesd2d(context_, 1, block_cols_A, b, 1, 0, 0);			// CB
//				Czgesd2d(ctxt, nr, nc, A_glob+N*c+r, N, sendr, sendc);		// From above
            }

            if(myid == 0)
            {
				Czgerv2d(ctxt, LLD, 1, &E[sendr * block_cols_A], LLD, sendr, 0);
//                Czgerv2d(ctxt, 1, block_cols_A, &E[sendr * block_cols_A], 1, sendr, 0);	// CB
	//			Czgerv2d(ctxt, nr, nc, A_loc+nrows*recvc+recvr, nrows, 0, 0);				// From above
            }
		}
		//		int noRowBlocks_	= Np;
	//	int noColBlocks_	= Np;
		//int block_size		= 2*pMax;
//		AlgebraS::solveMatrixVectorBP_AJ(	SS, 		Np*2*pMax, 		Np*2*pMax, B, E, pMax, noRowBlocks_, noColBlocks_);
		// ------------------------------------------------------------------
*/
/*
		// F - reconstruct the solved solution ------------------------------
		// Gather matrix 
		sendr = 0;
		for (int r = 0; r < N; r += Nb, sendr=(sendr+1)%procrows) {
			sendc = 0;
			// Number of rows to be sent
			// Is this the last row block?
			int nr = Nb;
			if (N-r < Nb)
				nr = N-r;
 
			for (int c = 0; c < M; c += Mb, sendc=(sendc+1)%proccols) {
				// Number of cols to be sent
				// Is this the last col block?
				int nc = Mb;
				if (M-c < Mb)
					nc = M-c;
 
				if (myrow == sendr && mycol == sendc) {
					// Send a nr-by-nc submatrix to process (sendr, sendc)
					Czgesd2d(ctxt, nr, nc, A_loc+nrows*recvc+recvr, nrows, 0, 0);
					recvc = (recvc+nc)%ncols;
				}
 
				if (mpiroot) {
					// Receive the same data
					// The leading dimension of the local matrix is nrows!
					Czgerv2d(ctxt, nr, nc, A_glob2+N*c+r, N, sendr, sendc);
				}
 
			}
 
			if (myrow == sendr)
				recvr = (recvr+nr)%nrows;
		}
 /*
		// Print test matrix
		Cblacs_barrier(ctxt, "All");
		if (mpiroot) {
			cout << "Matrix A_glob re-constructed:\n";
			for (int r = 0; r < N; ++r) {
				for (int c = 0; c < M; ++c) {
					cout <<" "<< *(A_glob2+N*c+r) << " ";
				}
				cout << endl;
			}
		}
*/
		// ------------------------------------------------------------------
		// Release resources
		delete[] A_glob;
//		delete[] A_glob2;
		delete[] A_loc;

		Cblacs_gridexit(ctxt);
		// solve in parallel ------------------------------------------------------------

		// op time ----------------------------------------------------------------------
		time(&ftime);                                       // calculate final time
		drift = difftime(ftime, itime);
		hours=int(drift/3600);
		mins =int(60*((drift/3600)-hours));
		secs =int(60*(60*((drift/3600)-hours)-mins));
		if (l==0 && mpirank==0){
			cout<<"2 - Solve system of equations execution time\t"<<drift<<endl;
			cout<<"2 - Solve system of equations execution time\t"<<hours<<" hrs "<<mins<<" mins "<<secs<<" secs" << endl;
		}
		//-------------------------------------------------------------------------------------------
/*
		// ------------------------------------------------------------------------------------------
		time(&itime);                                       // calculate initial time
		if (l==0 && mpirank==0) cout<<"3 - compute                    	: F = [TT].E"<<endl;
		// 4.5 matrix vector multiplication to finally obtain F : F = [TT].E ------------------------
		Algebra::multiplyVectorMatrix(	TT, 		Np*2*pMax,		Np*2*pMax,
											E,				F, 			1., 	0.);
		// op time --------------------
		time(&ftime);                                       // calculate final time
		drift = difftime(ftime, itime);
		hours=int(drift/3600);
		mins =int(60*((drift/3600)-hours));
		secs =int(60*(60*((drift/3600)-hours)-mins));
		if (l==0 && mpirank==0) cout<<"3 - compute system of equations execution time\t"<<hours<<" hrs "<<mins<<" mins "<<secs<<" secs" << endl;
		// ------------------------------------------------------------------------------------------
*/
/*
		// Calculate cross sections Cext, Cabs, Csca ------------------------------------------------------------------------------------------
		// 7 - computed from incident & scattered coeffs
		Cext=complex<double>(0.0, 0.0);
		Cabs=complex<double>(0.0, 0.0);
		Csca=complex<double>(0.0, 0.0);
		Csca_validate=complex<double>(0.0, 0.0);
		Qext=complex<double>(0.0, 0.0);
		Qabs=complex<double>(0.0, 0.0);
		Qsca=complex<double>(0.0, 0.0);
		complex<double> Cext_E[Cext_nMax+1];
		complex<double> Cext_H[Cext_nMax+1];
		int cext_opi=0;
		for(cext_opi=0; cext_opi<=Cext_nMax; cext_opi++){
			Cext_E[cext_opi]=complex<double>(0.0, 0.0);
			Cext_H[cext_opi]=complex<double>(0.0, 0.0);
		}
		complex<double> *F_translated;
		F_translated = new complex<double> [2*pMax];							// Scattered coeffs solution
		// 7.1 - Calculate abs cross section Cabs --------------------------------------------------------------------------------------------------
		time(&itime);                                       					// calculate initial time
		for(j=0; j<Np; j++)
		{
			// get Beta translation of scattered fields
			// A - Beta terms ---------------------------------------------------------------
//			R_rel = Origin - particles[j].vR;
			R_rel = particles[j].vR - Origin;			// Mackowski
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
//			for(i=0; i<2*pMax; i++)
//			{
//				F_local[i] = F[i+j*2*pMax];
//			}
			for(i=0; i<pMax; i++)				// Mackowski
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
//				Cext += real(	conj(excitation.dataIncAp[i]) * F_translated[i+0*pMax]
//							 +  conj(excitation.dataIncBp[i]) * F_translated[i+1*pMax]);
				// Mackowski
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
		// -----------------------------------------------------------------------------------------------------------------------------------------
		// op time --------------------
		time(&ftime);                                       // calculate final time
		drift = difftime(ftime, itime);
		hours=int(drift/3600);
		mins =int(60*((drift/3600)-hours));
		secs =int(60*(60*((drift/3600)-hours)-mins));
		if (l==0 && mpirank==0) cout<<"4 - Cext TE/TM individual (n,m) execution time\t"<<hours<<" hrs "<<mins<<" mins "<<secs<<" secs" << endl;
		//-------------------------------------------------------------------------------------------
/*
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
*/
//		// 7.3 - Calculate sca cross section Csca --------------------------------------------------------------------------------------------------
//		for(j=0; j<Np; j++)
//		{
//			for(int ll=0; ll<Np; ll++)
//			{
//				// get Beta translation of scattered fields
//				// A - Beta terms ---------------------------------------------------------------
//				R_rel = particles[j].vR - particles[ll].vR;
//				TC_beta.init(R_rel, waveK_0, BHreg_inc, nMax);
//				TC_beta.populate();
//				// ------------------------------------------------------------------------------
//				// ------------------------------------------------------------------------------
//				// B - map beta terms obtained from Coupling.cpp to account for efficient multiplication
//				for(ii=0; ii<pMax; ii++){
//					for(jj=0; jj<pMax; jj++){
//						// beta terms -----------------------------------------------------------
//						// - Q11
//						TC_betaa_t[ii][jj]					= TC_beta.dataApq[jj][ii];
//						// - Q12
//						TC_betaa_t[ii][jj+pMax]				= TC_beta.dataBpq[jj][ii];
//						// - Q21
//						TC_betaa_t[ii+pMax][jj]				= TC_beta.dataBpq[jj][ii];
//						// - Q22
//						TC_betaa_t[ii+pMax][jj+pMax]		= TC_beta.dataApq[jj][ii];
//						// ----------------------------------------------------------------------
//					}
//				}
//				// ------------------------------------------------------------------------------
//				// C - translate scattering coeffs
//				complex<double> *F_local = new complex<double>[2*pMax];
//				for(i=0; i<2*pMax; i++)
//				{
//					F_local[i] = F[i+ll*2*pMax];
//				}
//				Algebra::multiplyVectorMatrix(TC_betaa_t, 2*pMax, 2*pMax, F_local, F_translated, consC1, consC0);
//				delete [] F_local;
//				// D - calculate Cext - Stout 2001 (58)
//				for(i=0; i<pMax; i++)
//				{
//					Csca += real(	conj(F[j*2*pMax+i+0*pMax]) * F_translated[i+0*pMax]
//								 +  conj(F[j*2*pMax+i+1*pMax]) * F_translated[i+1*pMax]);
//				}
//			}	// l loop
//		}	// j loop
//		Csca *= 1. / (waveK_0 * waveK_0);
//		// -----------------------------------------------------------------------------------------------------------------------------------------
/*
		Csca_validate = Cext - Cabs;

		// normalise Cext, Cabs and Csca
//		Qext=Cext/(particles[0].radius*particles[0].radius);
//		Qabs=Cabs/(particles[0].radius*particles[0].radius);
//		Qsca=Csca/(particles[0].radius*particles[0].radius);
//		op_CextCabs<<lam_0<<" "<<Cext.real()<<" "<<Cabs.real()<<" "<<Csca.real()<<endl;
		op_CextCabs<<lam_0<<" "<<Cext.real()<<" "<<Cabs.real();
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
*/
		// at each scan delete all matrices -------------------------------------------------------
		delete[] B;
		delete[] E;
		delete[] B0;
		delete[] Cabs_aux;
		//delete[] F_translated;
		for(int i=0; i<2*pMax; i++){
			delete[] T_matrix[i];
			delete[] TC_betaa_t[i];		delete[] TC_alpha_t[i];		delete[] mT_TC_alpha[i];
		}
			delete[] T_matrix;
			delete[] TC_betaa_t;		delete[] TC_alpha_t;		delete[] mT_TC_alpha;

		for(int i=0; i<Np*2*pMax; i++){
			delete[] SS[i];			//delete[] TT[i];
		}
			delete[] SS;			//delete[] TT;

	} // for (l scan)
	// finish scan here -------------------------------------------------------------------------
	MPI_Finalize();

	// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
	// pre-plotting ----------------------------------------------------------------------------
	// 5 - obtain internal coeffs for each particle --------------------------------------------
	complex<double> *I_aux, *I_aux_j;
	I_aux = new complex<double> [2*Np*pMax];
	I_aux_j = new complex<double> [2*pMax];
	for(i=0; i<p.max(nMax); i++){
		I_aux_j[i] = complex<double>(0., 0.);
	}
	for (j=0; j<Np; j++){											// loop over all particles
		// Aux coeffs
//		getIaux (omega, particles[j], nMax, I_aux_j);
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
//	double ymin=-500.*1e-9, ymax=500.*1e-9;
	double zmin=-4000.*1e-9, zmax=4000.*1e-9;
	double dx = (xmax-xmin)/steps;
//	double dy = 0.;
	double dz = (zmax-zmin)/steps;
	Spherical<double> R_op, R_j;
	int Np_in=0;
	SphericalP<complex<double> > sum_up, sum_par, sum_sca, sum_inc, sum_ext, sum_int, sum_tot;

//	// the-phi - 2D plot ----------------------------------------------------------------------------------------------------------------
//	time(&itime);                                       // calculate initial time
//	//op ext mag
//	ofstream E1_err_mag("E1_err_mag");
//	ofstream E2_err_mag("E2_err_mag");
//	ofstream E3_err_mag("E3_err_mag");
//	Spherical<double> R_surfS(0., 0., 0.);
//	Cartesian<double> R_surfC(0., 0., 0.), vR_Cart(0., 0., 0.), R_opC(0., 0., 0.);
//	for(ii=1; ii<the_range; ii++){
//		for(jj=1; jj<phi_range; jj++){
//
//			// calculate incident, scattered, external & total fields -------------------------------------------------------------------
//			sum_inc=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_sca=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_ext=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_tot=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//
//			// Spherical observation point ----------------------------------------------------------------------------------------------
//			// observation at particle Obs_part defined above
//			R_surfS.rrr = particles[Obs_part].radius;		// force observation particle - first particle for now
//			R_surfS.the = consPi*(double(ii)/180.);
//			R_surfS.phi = consPi*(double(jj)/180.);
//
//			// translate R_surfS local to Global
//			// 1 - cart location on R_surf
//			R_surfC = Tools::toCartesian(R_surfS);
//			// 2 - cart location at center of observation particle j
//			vR_Cart = Tools::toCartesian(particles[Obs_part].vR);
//			// 3 - final location
//			R_opC.x = vR_Cart.x + R_surfC.x;
//			R_opC.y = vR_Cart.y + R_surfC.y;
//			R_opC.z = vR_Cart.z + R_surfC.z;
//			// Spherical projection
//			R_op=Tools::toSpherical(Cartesian<double>(R_opC.x, R_opC.y, R_opC.z));
//
//			// 0 - define observation particle ------------------------------------------------------------------------------------------
//			Np_in=Obs_part;									// define observation particle
//
//			// 1 - internal field -------------------------------------------------------------------------------------------------------
//			// if inside any particle - contribution form one particle only (Np_in) -----------------------------------------------------
//			sum_int=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			if(Np_in!=-1){
//				j=Np_in;														// assign j to Np_in
//				// scattering filed Aux functions
//				waveK_j = omega *sqrt(particles[j].elmag.epsilon*particles[j].elmag.epsilon_r * particles[j].elmag.mu*particles[j].elmag.mu_r);	// particle wavenumber
//				R_j = Tools::toPoint(R_op, particles[j].vR);					// local observation point
//				AuxCoefficients AuxCoeff_int;
//				AuxCoeff_int.init(R_j, waveK_j, BHreg_inc, nMax);
//				AuxCoeff_int.populate();
//				// Calculate internal field ---------------------------------------------------------------------------------------------
//				sum_int=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//				for(i=0; i<p.max(nMax); i++){
//					// internal from scattered coeffs 	- particle j
//					sum_int.rrr += I_aux[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].rrr + I_aux[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].rrr;
//					sum_int.the += I_aux[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].the + I_aux[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].the;
//					sum_int.phi += I_aux[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].phi + I_aux[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].phi;
//					// only to check plane wave is uniformally extended to all particles
////					sum_int.rrr += BB[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].rrr + BB[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].rrr;
////					sum_int.the += BB[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].the + BB[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].the;
////					sum_int.phi += BB[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].phi + BB[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].phi;
//				}
//			}
//			// if inside sphere ---------------------------------------------------------------------------------------------------------
//
//			// 2 - external field -------------------------------------------------------------------------------------------------------
//			// if outside all spheres ---------------------------------------------------------------------------------------------------
////			else{
//			// 2.1 scattered field from all particles -------------------------------------------------------------------------------
//			for (j=0; j<Np; j++){											// scattered field contribution from all particles
//				// scattering filed Aux functions
//				R_j = Tools::toPoint(R_op, particles[j].vR);				// local observation point
//				AuxCoefficients AuxCoeff_sca;
//				AuxCoeff_sca.init(R_j, waveK_0, BHreg_sca, nMax);
//				AuxCoeff_sca.populate();
//				// Calculate scattered field ----------------------------------------------------------------------------------------
//				sum_par=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//				for(i=0; i<p.max(nMax); i++){
//					// sca from sphere 1
//					sum_par.rrr += F[j*2*pMax+i+0*pMax]*AuxCoeff_sca.dataMp[i].rrr + F[j*2*pMax+i+1*pMax]*AuxCoeff_sca.dataNp[i].rrr;
//					sum_par.the += F[j*2*pMax+i+0*pMax]*AuxCoeff_sca.dataMp[i].the + F[j*2*pMax+i+1*pMax]*AuxCoeff_sca.dataNp[i].the;
//					sum_par.phi += F[j*2*pMax+i+0*pMax]*AuxCoeff_sca.dataMp[i].phi + F[j*2*pMax+i+1*pMax]*AuxCoeff_sca.dataNp[i].phi;
//				}
//				// get total scattered - accumulation from all particles
//				sum_sca = sum_sca + sum_par;							// total scattered field includes scattered scattered from all
//			}	// for loop (j)
//			// ----------------------------------------------------------------------------------------------------------------------
//
//			// 2.2 global incident field --------------------------------------------------------------------------------------------
//			AuxCoefficients AuxCoeff_inc;
//			AuxCoeff_inc.init(R_op, waveK_0, BHreg_inc, nMax);
//			AuxCoeff_inc.populate();
//			// Calculate incident field ---------------------------------------------------------------------------------------------
//			sum_inc=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			for(i=0; i<p.max(nMax); i++){
////				// local to particle
////				sum_inc.rrr		+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].rrr + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].rrr;
////				sum_inc.the 	+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].the + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].the;
////				sum_inc.phi 	+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].phi + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].phi;
//				// local to system origin
//				sum_inc.rrr		+= 	Inc[i+0*pMax]*AuxCoeff_inc.dataMp[i].rrr + Inc[i+1*pMax]*AuxCoeff_inc.dataNp[i].rrr;
//				sum_inc.the 	+= 	Inc[i+0*pMax]*AuxCoeff_inc.dataMp[i].the + Inc[i+1*pMax]*AuxCoeff_inc.dataNp[i].the;
//				sum_inc.phi 	+= 	Inc[i+0*pMax]*AuxCoeff_inc.dataMp[i].phi + Inc[i+1*pMax]*AuxCoeff_inc.dataNp[i].phi;
//			}
//			// ----------------------------------------------------------------------------------------------------------------------
//
//			// 2.3 total external field = global incident field + scattered field from all particles --------------------------------
//			sum_ext = sum_inc + sum_sca;								// total external field includes total scattered + incident
////			}	// else if(!Np_in)
//			// if outside all spheres ---------------------------------------------------------------------------------------------------
//
//			// 3 - total = external (or) internal field ---------------------------------------------------------------------------------
//			sum_tot = sum_ext + sum_int;						// total field will include either internal for a given particle, or external from all
////			// from Cart to Spherical conversion
//			sum_ext = Tools::fromProjection(R_surfS, sum_ext);
//			sum_int = Tools::fromProjection(R_surfS, sum_int);
//			sum_int.rrr*=eps_i;					// continuity for D=epsi_r*E
////			sum_inc = Tools::fromProjection(R_surfS, sum_inc);		// projection here is irrelevant as it is a plane-wave
////			sum_tot = Tools::fromProjection(R_surfS, sum_tot);		// projection here illustrates filed profile with respect to a given coordinate system
//
//			// op -------------------------------------------------------------------
////			cout<<"computed "<< ii<<" out of a total of "<<180<<endl;
//			E1_err_mag<<abs(abs(sum_ext.rrr)-abs(sum_int.rrr))/abs(sum_int.rrr)<<" ";
//			E2_err_mag<<abs(abs(sum_ext.the)-abs(sum_int.the))/abs(sum_int.the)<<" ";
//			E3_err_mag<<abs(abs(sum_ext.phi)-abs(sum_int.phi))/abs(sum_int.phi)<<" ";
//		}
//		E1_err_mag<<endl;
//		E2_err_mag<<endl;
//		E3_err_mag<<endl;
//	}
//	// op time --------------------
//	time(&ftime);                                       // calculate final time
//	drift = difftime(ftime, itime);
//	hours=int(drift/3600);
//	mins =int(60*((drift/3600)-hours));
//	secs =int(60*(60*((drift/3600)-hours)-mins));
//	cout<<"4 - Fields continuity check execution time    \t"<<hours<<" hrs "<<mins<<" mins "<<secs<<" secs" << endl;
//	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------- //

	if(mpirank==0){
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
//	//op ext mag
//	ofstream E1_tot_mag("E1_tot_mag");
//	ofstream E2_tot_mag("E2_tot_mag");
//	ofstream E3_tot_mag("E3_tot_mag");

	// xy plane
	z=0.*1e-9;
	for(ii=0, x=xmin; ii<=steps; ii++, x+=dx){
		for(jj=0, y=zmin; jj<=steps; jj++, y+=dx){

	// xz plane
//	y=0.*1e-9;
////	x=0.*1e-9;
//	for(ii=0, x=xmin; ii<=steps; ii++, x+=dx){
//		for(jj=0, z=zmin; jj<=steps; jj++, z+=dz){

			// calculate incident, scattered, external & total fields -------------------------------------------------------------------
			sum_inc=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
			sum_sca=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
			sum_ext=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
			sum_tot=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));

			// Cartesian observation point ----------------------------------------------------------------------------------------------
			R_op=Tools::toSpherical(Cartesian<double>(x+1e-12, y+1e-12, z+1e-12));
//			R_op=Tools::toSpherical(Cartesian<double>(x, y, z));

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
					// internal from scattered coeffs 	- particle j
					sum_int.rrr += I_aux[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].rrr + I_aux[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].rrr;
					sum_int.the += I_aux[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].the + I_aux[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].the;
					sum_int.phi += I_aux[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].phi + I_aux[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].phi;
					// only to check plane wave is uniformally extended to both particles
//					sum_int.rrr += BB[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].rrr + BB[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].rrr;
//					sum_int.the += BB[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].the + BB[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].the;
//					sum_int.phi += BB[j*2*pMax+i+0*pMax]*AuxCoeff_int.dataMp[i].phi + BB[j*2*pMax+i+1*pMax]*AuxCoeff_int.dataNp[i].phi;
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
						// sca from sphere 1
						sum_par.rrr += F[j*2*pMax+i+0*pMax]*AuxCoeff_sca.dataMp[i].rrr + F[j*2*pMax+i+1*pMax]*AuxCoeff_sca.dataNp[i].rrr;
						sum_par.the += F[j*2*pMax+i+0*pMax]*AuxCoeff_sca.dataMp[i].the + F[j*2*pMax+i+1*pMax]*AuxCoeff_sca.dataNp[i].the;
						sum_par.phi += F[j*2*pMax+i+0*pMax]*AuxCoeff_sca.dataMp[i].phi + F[j*2*pMax+i+1*pMax]*AuxCoeff_sca.dataNp[i].phi;
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
	//				// local to particle
	//				sum_inc.rrr		+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].rrr + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].rrr;
	//				sum_inc.the 	+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].the + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].the;
	//				sum_inc.phi 	+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].phi + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].phi;
					// local to system origin
					sum_inc.rrr		+= 	Inc[i+0*pMax]*AuxCoeff_inc.dataMp[i].rrr + Inc[i+1*pMax]*AuxCoeff_inc.dataNp[i].rrr;
					sum_inc.the 	+= 	Inc[i+0*pMax]*AuxCoeff_inc.dataMp[i].the + Inc[i+1*pMax]*AuxCoeff_inc.dataNp[i].the;
					sum_inc.phi 	+= 	Inc[i+0*pMax]*AuxCoeff_inc.dataMp[i].phi + Inc[i+1*pMax]*AuxCoeff_inc.dataNp[i].phi;
				}
				// ----------------------------------------------------------------------------------------------------------------------

				// 2.3 total external field = global incident field + scattered field from all particles --------------------------------
				sum_ext = sum_inc + sum_sca;								// total external field includes total scattered + incident

			}	// else if(!Np_in)
			// if outside all spheres ---------------------------------------------------------------------------------------------------

			// 3 - total = external (or) internal field ---------------------------------------------------------------------------------
			sum_tot = sum_ext + sum_int;						// total field will include either internal for a given particle, or external from all
//			// from Cart to Spherical conversion
//			sum_inc = Tools::fromProjection(R_up, sum_inc);		// projection here is irrelevant as it is a plane-wave
//			sum_tot = Tools::fromProjection(R_up, sum_tot);		// projection here illustrates filed profile with respect to a given coordinate system

			// inc
			E1_inc_real<<real(sum_inc.rrr)<<" ";		E1_inc_imag<<imag(sum_inc.rrr)<<" ";
			E2_inc_real<<real(sum_inc.the)<<" ";		E2_inc_imag<<imag(sum_inc.the)<<" ";
			E3_inc_real<<real(sum_inc.phi)<<" ";		E3_inc_imag<<imag(sum_inc.phi)<<" ";
			// tot
//			cout<<imag(sum_tot.rrr)<<" ";
			E1_real<<real(sum_tot.rrr)<<" ";		E1_imag<<imag(sum_tot.rrr)<<" ";
			E2_real<<real(sum_tot.the)<<" ";		E2_imag<<imag(sum_tot.the)<<" ";
			E3_real<<real(sum_tot.phi)<<" ";		E3_imag<<imag(sum_tot.phi)<<" ";
			// tot mag
//			E1_tot_mag<<abs(sum_tot.rrr)<<" ";
//			E2_tot_mag<<abs(sum_tot.the)<<" ";
//			E3_tot_mag<<abs(sum_tot.phi)<<" ";
//			cout<<"computed "<< jj+ii*(steps+1)+1<<" out of a total of "<<(steps+1)*(steps+1)<<endl;
		}
		// inc
		E1_inc_real<<endl;		E1_inc_imag<<endl;
		E2_inc_real<<endl;		E2_inc_imag<<endl;
		E3_inc_real<<endl;		E3_inc_imag<<endl;
		// tot
//		cout<<endl;
		E1_real<<endl;		E1_imag<<endl;		//E1_mag<<endl;
		E2_real<<endl;		E2_imag<<endl;		//E2_mag<<endl;
		E3_real<<endl;		E3_imag<<endl;		//E3_mag<<endl;
	}
	// op time --------------------
//	E1_imag.flush();	E1_imag.close();	
	time(&ftime);                                       // calculate final time
	drift = difftime(ftime, itime);
	hours=int(drift/3600);
	mins =int(60*((drift/3600)-hours));
	secs =int(60*(60*((drift/3600)-hours)-mins));
	cout<<"5 - Fields profile execution time             \t"<<hours<<" hrs "<<mins<<" mins "<<secs<<" secs" << endl;
	// ---------------------------------------------------------------------------------------------------------------------------------- //
	}	// if(mpirank==0)

//	// at maximum point from above
//	for(ii=1; ii<2; ii++){
//		for(jj=1; jj<2; jj++){
//
//			ii = 89;
//			jj = 89;
//
//			// up, down & total points
//			sum_up=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_do=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_sca=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_inc=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_ext=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_int=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_tot=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//
//			R_op=Tools::toSpherical(Cartesian<double>(x+1e-12, y+1e-12, z+1e-12));
//			R_op.rrr = particles[0].radius;			// force to fall on radius
//			R_op.the = consPi*ii/180;			// force to fall on radius
//			R_op.phi = consPi*jj/180;			// force to fall on radius
//
////			// incident field Aux functions
////			AuxCoeff_inc.init(R_op, waveK_0, BHreg_inc, nMax);
////			AuxCoeff_inc.populate();
//
//			// check if point falls outside spheres
//			// inparticle0=checkInner(R_op, particles[0], dx);
//			 inparticle0_just_below=checkInner_just_below(R_op, particles[0], dx);
//			 inparticle0_just_above=checkInner_just_above(R_op, particles[0], dx);
//
////			if (inparticle0_just_above==1){
//
//			// outside sphere -----------------------------------------------------------------------------------------------------------
////			if(inparticle0==-1){
////				if(inparticle0_just_above==1){
//
//				// scattering filed Aux functions
//				//Translate to object 0
//				R_up = Tools::toPoint(R_op, particles[0].vR);
//				R_up.rrr = particles[0].radius;			// force to fall on radius
////					R_up.phi = consPi*150/180;			// force to fall on radius
//
//				// incident field Aux functions
//				AuxCoeff_inc.init(R_up, waveK_0, BHreg_inc, nMax);
//				AuxCoeff_inc.populate();
//
//				// outside sphere
//				AuxCoeff_Rup.init(R_up, waveK_0, BHreg_sca, nMax);
//				AuxCoeff_Rup.populate();
//
//				// inside sphere
//				AuxCoeff_Rdo.init(R_up, waveK_i, BHreg_inc, nMax);
//				AuxCoeff_Rdo.populate();
//
//				// inc, sca & int fields --------------------------------------------------------------------------------------------------
//				for(i=0; i<p.max(nMax); i++){
//
//					p=i;
//
//					// inc
//					sum_inc.rrr += excitation.dataIncAp[i]*AuxCoeff_inc.dataMp[i].rrr + excitation.dataIncBp[i]*AuxCoeff_inc.dataNp[i].rrr;
//					sum_inc.the += excitation.dataIncAp[i]*AuxCoeff_inc.dataMp[i].the + excitation.dataIncBp[i]*AuxCoeff_inc.dataNp[i].the;
//					sum_inc.phi += excitation.dataIncAp[i]*AuxCoeff_inc.dataMp[i].phi + excitation.dataIncBp[i]*AuxCoeff_inc.dataNp[i].phi;
//
//					// sca from sphere 1
//					sum_up.rrr += FF[i+0*p.max(nMax)]*AuxCoeff_Rup.dataMp[i].rrr + FF[i+1*p.max(nMax)]*AuxCoeff_Rup.dataNp[i].rrr;
//					sum_up.the += FF[i+0*p.max(nMax)]*AuxCoeff_Rup.dataMp[i].the + FF[i+1*p.max(nMax)]*AuxCoeff_Rup.dataNp[i].the;
//					sum_up.phi += FF[i+0*p.max(nMax)]*AuxCoeff_Rup.dataMp[i].phi + FF[i+1*p.max(nMax)]*AuxCoeff_Rup.dataNp[i].phi;
//
//					// sca from sphere 1
//					sum_do.rrr += I_aux_1[i+0*p.max(nMax)]*AuxCoeff_Rdo.dataMp[i].rrr + I_aux_1[i+1*p.max(nMax)]*AuxCoeff_Rdo.dataNp[i].rrr;
//					sum_do.the += I_aux_1[i+0*p.max(nMax)]*AuxCoeff_Rdo.dataMp[i].the + I_aux_1[i+1*p.max(nMax)]*AuxCoeff_Rdo.dataNp[i].the;
//					sum_do.phi += I_aux_1[i+0*p.max(nMax)]*AuxCoeff_Rdo.dataMp[i].phi + I_aux_1[i+1*p.max(nMax)]*AuxCoeff_Rdo.dataNp[i].phi;
//
//					// get external field
//					sum_ext = sum_inc + sum_up;
//
//					// from Cart to Spherical conversion
//					sum_ext = Tools::fromProjection(R_up, sum_ext);
//					sum_int = Tools::fromProjection(R_up, sum_do);
//					sum_int.rrr*=eps_i;
//
//					if(p.first == -p.second){
////						E1_tot_mag<<p.first<<" "<<abs(abs(sum_ext.rrr)-abs(sum_int.rrr))/abs(sum_int.rrr)<<endl;
////						E2_tot_mag<<p.first<<" "<<abs(abs(sum_ext.the)-abs(sum_int.the))/abs(sum_int.the)<<endl;
////						E3_tot_mag<<p.first<<" "<<abs(abs(sum_ext.phi)-abs(sum_int.phi))/abs(sum_int.phi)<<endl;
//						E1_tot_mag<<p.first<<" "<<abs(abs(sum_ext.rrr))<<" "<<abs(abs(sum_int.rrr))<<" "<<abs(abs(sum_ext.rrr)-abs(sum_int.rrr))/abs(sum_int.rrr)<<endl;
//						E2_tot_mag<<p.first<<" "<<abs(abs(sum_ext.the))<<" "<<abs(abs(sum_int.the))<<" "<<abs(abs(sum_ext.the)-abs(sum_int.the))/abs(sum_int.the)<<endl;
//						E3_tot_mag<<p.first<<" "<<abs(abs(sum_ext.phi))<<" "<<abs(abs(sum_int.phi))<<" "<<abs(abs(sum_ext.phi)-abs(sum_int.phi))/abs(sum_int.phi)<<endl;
//					}
//
//				}
//			    // else if inside sphere -----------------------------------------------------------------------------------------------------
//
//			// op
//			cout<<"computed "<< ii<<" out of a total of "<<90<<endl;
//
//		}
//		E1_tot_mag<<endl;
//		E2_tot_mag<<endl;
//		E3_tot_mag<<endl;
//	}
//	// ------------------------------------------------------------------------------------------------------------------------------------------------------------------- //

	if(mpirank==0) cout<<"done"<<endl;
	return 0;
}

//	// ***************************************************************************************************************************
//	// Mishchenko method - oneSphere system - matrix vector multiplication to obtain FF - scattering coeffs
//	complex<double> *B0,  *FF;
//	B0 = new complex<double> [2*p.max(nMax)];
//	FF = new complex<double> [2*p.max(nMax)];
//	for(i=0; i<2*p.max(nMax); i++){		// p
//		B0[i]=complex<double>(0., 0.);
//		FF[i]=complex<double>(0., 0.);
//	}
//	// Translate incident coeffs to coordinate system of particle -------------------------------------------------
////	for(i=0; i<p.max(nMax); i++){		// p
////		for(j=0; j<p.max(nMax); j++){	// q
////			B0[i+0*p.max(nMax)] += (TC_b0g.dataApq[j][i]*excitation.dataIncAp[j] + TC_b0g.dataBpq[j][i]*excitation.dataIncBp[j]);
////			B0[i+1*p.max(nMax)] += (TC_b0g.dataBpq[j][i]*excitation.dataIncAp[j] + TC_b0g.dataApq[j][i]*excitation.dataIncBp[j]);
////		}
////	}
//
//	// LAPACK multiplication ---------------------------------------------------------------------------------------
//	// obtain local excitation
//	Algebra::multiplyVectorMatrix(	TC_beta_1g, 2*p.max(nMax),	2*p.max(nMax),
//									Inc,		B0, 			1., 			0.);	// using transfer
//	// obtain local scattered
//	Algebra::multiplyVectorMatrix(	T_1, 		2*p.max(nMax),	2*p.max(nMax),
//									B0,			FF, 			1., 			0.);	// using transfer
//
//// obtain internal coeffs for each particle --------------------------------------------------------------------
//	complex<double> *I_aux_1;
//	I_aux_1 = new complex<double> [2*p.max(nMax)];
//	// Aux ceffs
//	getIaux (omega, particles[0], nMax, I_aux_1);
//	// internal coeffs
//	for(i=0; i<p.max(nMax); i++){
//		// particle 1
//		I_aux_1[i+0*p.max(nMax)] = I_aux_1[i+0*p.max(nMax)] * FF[i+0*p.max(nMax)];
//		I_aux_1[i+1*p.max(nMax)] = I_aux_1[i+1*p.max(nMax)] * FF[i+1*p.max(nMax)];
//	}
//
//	// 5 - produce scattered fields solution -----------------------------------------------------------------------
//	AuxCoefficients AuxCoeff_Rup, AuxCoeff_Rdo, AuxCoeff_inc;
//	int BHreg_sca=0;												// compute regular values of BH - non-regular
//	int BHreg_inc=1;												// compute regular values of BH - regular
//	// do a 2D scan of EF profile ----------------------------------------------------------------------------------
//	double xmin=-500.*1e-9, xmax=500.*1e-9;
////	double ymin=-500.*1e-9, ymax=500.*1e-9;
//	double zmin=-500.*1e-9, zmax=500.*1e-9;
//	int steps = 50;
//	double dx = (xmax-xmin)/steps;
////	double dy = 0.;
//	double dz = (zmax-zmin)/steps;
//	Spherical<double> R_op, R_up;
//	int inparticle0=0;
////	int inparticle0_just_below=0;
////	int inparticle0_just_above=0;
//	SphericalP<complex<double> > sum_up, sum_do, sum_sca, sum_inc, sum_ext, sum_int, sum_tot;
//	// start 2D scamn
//	// op inc
//	ofstream E1_inc_real("E1_inc_real");		ofstream E1_inc_imag("E1_inc_imag");
//	ofstream E2_inc_real("E2_inc_real");		ofstream E2_inc_imag("E2_inc_imag");
//	ofstream E3_inc_real("E3_inc_real");		ofstream E3_inc_imag("E3_inc_imag");
//	//op ext
//	ofstream E1_tot_real("E1_tot_real");		ofstream E1_tot_imag("E1_tot_imag");
//	ofstream E2_tot_real("E2_tot_real");		ofstream E2_tot_imag("E2_tot_imag");
//	ofstream E3_tot_real("E3_tot_real");		ofstream E3_tot_imag("E3_tot_imag");
//	//op ext mag
//	ofstream E1_tot_mag("E1_tot_mag");
//	ofstream E2_tot_mag("E2_tot_mag");
//	ofstream E3_tot_mag("E3_tot_mag");
//
//	Cartesian<double> Cart_R_op, Cart_R_temp;
//
//	// Cartesian 2D plot -----------------------------------------------------------------------------------------------------------
//	// xz plane
//	y=100.*1e-9;
//	y=0.*1e-9;
//	for(ii=0, x=xmin; ii<=steps; ii++, x+=dx){
//		for(jj=0, z=zmin; jj<=steps; jj++, z+=dz){
//
//	// yz plane
////	x=0.;
////	for(ii=0, y=xmin; ii<steps; ii++, y+=dx){
////		for(jj=0, z=zmin; jj<steps; jj++, z+=dz){
//
//			// up, down & total points
//			sum_up=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_sca=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_inc=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_ext=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_int=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//			sum_tot=SphericalP<complex<double> >(complex<double>(0.0, 0.0), complex<double>(0.0, 0.0), complex<double>(0.0, 0.0));
//
//			R_op=Tools::toSpherical(Cartesian<double>(x+1e-12, y+1e-12, z+1e-12));
////			R_op=Tools::toSpherical(Cartesian<double>(x, y, z));
//
//			// 1 - incident field ------------------------------------------------------------------------------------------------
//			// incident field Aux functions
////			R_up = Tools::toPoint(R_op, particles[0].vR);					// local observation point
//			AuxCoeff_inc.init(R_op, waveK_0, BHreg_inc, nMax);
//			AuxCoeff_inc.populate();
//			// Calculate incident field ------------------------------------------------------------------------------------------
//			for(i=0; i<p.max(nMax); i++){
////				// local to particle
////				sum_inc.rrr		+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].rrr + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].rrr;
////				sum_inc.the 	+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].the + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].the;
////				sum_inc.phi 	+= 	B0[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].phi + B0[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].phi;
//				// local to system origin
//				sum_inc.rrr		+= 	Inc[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].rrr + Inc[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].rrr;
//				sum_inc.the 	+= 	Inc[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].the + Inc[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].the;
//				sum_inc.phi 	+= 	Inc[i+0*p.max(nMax)]*AuxCoeff_inc.dataMp[i].phi + Inc[i+1*p.max(nMax)]*AuxCoeff_inc.dataNp[i].phi;
//			}
//			// get total sum
//			// --------------------------------------------------------------------------------------------------------------------
//
//			// 2 - scattered field ------------------------------------------------------------------------------------------------
//			// check if point falls outside spheres
//			inparticle0=checkInner(R_op, particles, Np);
//			// 2.1 - outside sphere -----------------------------------------------------------------------------------------------
//			if(inparticle0==-1){
//				// scattering filed Aux functions
//				R_up = Tools::toPoint(R_op, particles[0].vR);				// local observation point
//				AuxCoeff_Rup.init(R_up, waveK_0, BHreg_sca, nMax);
//				AuxCoeff_Rup.populate();
//
//				// Calculate scattered field --------------------------------------------------------------------------------------
//				for(i=0; i<p.max(nMax); i++){
//					// sca from sphere 1
//					sum_up.rrr += FF[i+0*p.max(nMax)]*AuxCoeff_Rup.dataMp[i].rrr + FF[i+1*p.max(nMax)]*AuxCoeff_Rup.dataNp[i].rrr;
//					sum_up.the += FF[i+0*p.max(nMax)]*AuxCoeff_Rup.dataMp[i].the + FF[i+1*p.max(nMax)]*AuxCoeff_Rup.dataNp[i].the;
//					sum_up.phi += FF[i+0*p.max(nMax)]*AuxCoeff_Rup.dataMp[i].phi + FF[i+1*p.max(nMax)]*AuxCoeff_Rup.dataNp[i].phi;
//				}
//				// get total sum
//				sum_sca = sum_up;
//				sum_ext = sum_inc + sum_up;		// total
//			}	// if(!inparticle)
//			// outside sphere ------------------------------------------------------------------------------------------------------------
//
//			// else if inside sphere -----------------------------------------------------------------------------------------------------
//			else if(inparticle0!=0){
//				// scattering filed Aux functions
//				waveK_j = omega *sqrt(particles[0].elmag.epsilon*particles[0].elmag.epsilon_r * particles[0].elmag.mu*particles[0].elmag.mu_r);	// particle wavenumber
//				R_up = Tools::toPoint(R_op, particles[0].vR);				// local observation point
//				AuxCoeff_Rup.init(R_up, waveK_j, BHreg_inc, nMax);
//				AuxCoeff_Rup.populate();
//
//				// inc, sca & int fields --------------------------------------------------------------------------------------------------
//				for(i=0; i<p.max(nMax); i++){
//					// int from scatterd coeffs 	- sphere 1
//					sum_up.rrr += I_aux_1[i+0*p.max(nMax)]*AuxCoeff_Rup.dataMp[i].rrr + I_aux_1[i+1*p.max(nMax)]*AuxCoeff_Rup.dataNp[i].rrr;
//					sum_up.the += I_aux_1[i+0*p.max(nMax)]*AuxCoeff_Rup.dataMp[i].the + I_aux_1[i+1*p.max(nMax)]*AuxCoeff_Rup.dataNp[i].the;
//					sum_up.phi += I_aux_1[i+0*p.max(nMax)]*AuxCoeff_Rup.dataMp[i].phi + I_aux_1[i+1*p.max(nMax)]*AuxCoeff_Rup.dataNp[i].phi;
//				}
//				// get total sum
//				sum_int = sum_up;			// must not add incoming wave here as internal field is inherently satisfied.
//			}
//			// else if inside sphere -----------------------------------------------------------------------------------------------------
//
//			// get total field -----------------------------------------------------------------------------------------------------------
//			sum_tot = sum_ext + sum_int;
////			// from Cart to Spherical conversion
////			sum_inc = Tools::fromProjection(R_up, sum_inc);		// projection here is irrelevant as it is a plane-wave
////			sum_tot = Tools::fromProjection(R_up, sum_tot);		// projection here illustrates filed profile with respect to a given coordinate system
//
//			// inc
//			E1_inc_real<<real(sum_inc.rrr)<<" ";		E1_inc_imag<<imag(sum_inc.rrr)<<" ";
//			E2_inc_real<<real(sum_inc.the)<<" ";		E2_inc_imag<<imag(sum_inc.the)<<" ";
//			E3_inc_real<<real(sum_inc.phi)<<" ";		E3_inc_imag<<imag(sum_inc.phi)<<" ";
//			// tot
//			E1_tot_real<<real(sum_tot.rrr)<<" ";		E1_tot_imag<<imag(sum_tot.rrr)<<" ";
//			E2_tot_real<<real(sum_tot.the)<<" ";		E2_tot_imag<<imag(sum_tot.the)<<" ";
//			E3_tot_real<<real(sum_tot.phi)<<" ";		E3_tot_imag<<imag(sum_tot.phi)<<" ";
//			// tot mag
//			E1_tot_mag<<abs(sum_tot.rrr)<<" ";
//			E2_tot_mag<<abs(sum_tot.the)<<" ";
//			E3_tot_mag<<abs(sum_tot.phi)<<" ";
//			cout<<"computed "<< jj+ii*(steps+1)+1<<" out of a total of "<<(steps+1)*(steps+1)<<endl;
//		}
//		// inc
//		E1_inc_real<<endl;		E1_inc_imag<<endl;
//		E2_inc_real<<endl;		E2_inc_imag<<endl;
//		E3_inc_real<<endl;		E3_inc_imag<<endl;
//		// tot
//		E1_tot_real<<endl;		E1_tot_imag<<endl;		E1_tot_mag<<endl;
//		E2_tot_real<<endl;		E2_tot_imag<<endl;		E2_tot_mag<<endl;
//		E3_tot_real<<endl;		E3_tot_imag<<endl;		E3_tot_mag<<endl;
//
//	}
//	// // oneSphere system ---------------------------------------------------------------------------------------------------------------------------------- //
//	// *****************************************************************************************************************************

// *********************************************************************************************************************************
// twoSpheres system - matrix vector multiplication to obtain FF - scattering coeffs
//	// Use this or the one below ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
//	// A - Original formulation --------------------------------------------------
//	// solve system of equations -------------------------------------------------
//	// 1 - auxiliary multiplication matrices
//	complex<double> **mT1_TC_alpha_01, **mT2_TC_alpha_10;
//	mT1_TC_alpha_01 = Tools::Get_2D_c_double(2*p.max(nMax), 2*p.max(nMax));
//	mT2_TC_alpha_10 = Tools::Get_2D_c_double(2*p.max(nMax), 2*p.max(nMax));
//	// obtain multiplication - -1[T_1][TC_alpha_01]
//	Algebra::multiplyMatrixMatrix(	T_1, 				2*p.max(nMax),	2*p.max(nMax),
//									TC_alpha_01, 		2*p.max(nMax),	2*p.max(nMax),
//									mT1_TC_alpha_01, 	-1., 			0.);
//	// obtain multiplication - -1[T_2][TC_alpha_10]
//	Algebra::multiplyMatrixMatrix(	T_2, 				2*p.max(nMax),	2*p.max(nMax),
//									TC_alpha_10, 		2*p.max(nMax),	2*p.max(nMax),
//									mT2_TC_alpha_10, 	-1., 			0.);
//	// 2 - global matrix --------------------------------------------------------
//	// populate SSSS ------------------------------------------------------------
//	complex<double> **SSSS;
//	SSSS = Tools::Get_2D_c_double(4*p.max(nMax), 4*p.max(nMax));
//	for(i=0; i<2*p.max(nMax); i++){
//		for(j=0; j<2*p.max(nMax); j++){
//			SSSS[i][j]	= complex<double>(0., 0.);
//		}
//	}
//	for(i=0; i<2*p.max(nMax); i++){
//		for(j=0; j<2*p.max(nMax); j++){
//			// - Q11
//			if(i==j)
//				SSSS[i][j]			 					=1.;
//			// - Q12
//			SSSS[i][j+2*p.max(nMax)]					=mT1_TC_alpha_01[i][j];
//			// - Q21
//			SSSS[i+2*p.max(nMax)][j]					=mT2_TC_alpha_10[i][j];
//			// - Q22
//			if(i+2*p.max(nMax)==j+2*p.max(nMax))
//				SSSS[i+2*p.max(nMax)][j+2*p.max(nMax)]	=1.;
//		}
//	}
//
//	// 3 - RHS vector ----------------------------------------------------------------
//	complex<double> *B0, *B1;
//	B0 = new complex<double> [2*p.max(nMax)];
//	B1 = new complex<double> [2*p.max(nMax)];
//	// excitation
//	complex<double> *Exc;
//	Exc = new complex<double> [2*p.max(nMax)];
//	SphericalP<complex<double> > Einc(1., 0., 0.);
//	Spherical<double> vKInc;
//	vKInc.rrr = real(waveK_0);
//	vKInc.the = consPi/2;
//	vKInc.phi = consPi/2;
//
//	Excitation excitation;
//	int Type=0;												// compute regular values of BH - regular
//	excitation.init(Type, Einc, vKInc, nMax);
//	excitation.populate();
//	for(i=0; i<p.max(nMax); i++){
//		Exc[i]=excitation.dataIncAp[i];
//		Exc[i+p.max(nMax)]=excitation.dataIncBp[i];
//	}
//
//	complex<double> *BB, *FF;
//	BB = new complex<double> [4*p.max(nMax)];
//	FF = new complex<double> [4*p.max(nMax)];
//	// auxiliary multiplication matrices
//	complex<double> **T1_TC_beta_0g, **T2_TC_beta_1g;
//	T1_TC_beta_0g = Tools::Get_2D_c_double(2*p.max(nMax), 2*p.max(nMax));
//	T2_TC_beta_1g = Tools::Get_2D_c_double(2*p.max(nMax), 2*p.max(nMax));
//	// obtain multiplication - -1[T_1][TC_alpha_01]
//	Algebra::multiplyMatrixMatrix(	T_1, 				2*p.max(nMax),	2*p.max(nMax),
//									TC_beta_0g, 		2*p.max(nMax),	2*p.max(nMax),
//									T1_TC_beta_0g, 		1., 			0.);
//	// obtain multiplication - -1[T_2][TC_alpha_10]
//	Algebra::multiplyMatrixMatrix(	T_2, 				2*p.max(nMax),	2*p.max(nMax),
//									TC_beta_1g, 		2*p.max(nMax),	2*p.max(nMax),
//									T2_TC_beta_1g, 		1., 			0.);
//	// compute B0 - Implements the multiplication Y = alpha*A*X + beta*Y
//	Algebra::multiplyVectorMatrix(	T1_TC_beta_0g, 		2*p.max(nMax),	2*p.max(nMax),
//									Exc,				B0, 			1., 	0.);
//	// compute B1 - Implements the multiplication Y = alpha*A*X + beta*Y
//	Algebra::multiplyVectorMatrix(	T2_TC_beta_1g, 		2*p.max(nMax),	2*p.max(nMax),
//									Exc,				B1, 			1., 	0.);
//	for(i=0; i<2*p.max(nMax); i++){
//		BB[i]=B0[i];
//		BB[i+2*p.max(nMax)]=B1[i];
//	}
//
//	// final solution - Solves the SSSS*FF = BB equation
//	Algebra::solveMatrixVector(SSSS, 4*p.max(nMax), 4*p.max(nMax), BB, FF);
//	// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
