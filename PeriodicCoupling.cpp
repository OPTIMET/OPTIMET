/*
 * Periodic_Couplin.cpp
 *
 *  Created on: June 27, 2013
 *      Author: uceealj
 */

#include "PeriodicCoupling.h"

#include "Tools.h"
#include "Legendre.h"
#include "Bessel.h"
#include "Symbol.h"
#include "CompoundIterator.h"
#include "constants.h"
#include "AJ_AuxFuns.h"

#include <cmath>

#include "gsl/gsl_sf_gamma.h"

double PeriodicCoupling::a_nm_p(double n, double m)
{
	return -sqrt( ((n+m+1.)*(n-m+1.))/((2.*n+1.)*(2.*n+3.)) );
}

double PeriodicCoupling::a_nm_m(double n, double m)
{
	return +sqrt( ((n+m)*(n-m))/((2.*n+1.)*(2.*n-1.)) );
}

double PeriodicCoupling::b_nm_p(double n, double m)
{
	return +sqrt( ((n+m+2.)*(n+m+1.))/((2.*n+1.)*(2.*n+3.)) );
}

double PeriodicCoupling::b_nm_m(double n, double m)
{
	return +sqrt( ((n-m)*(n-m-1.))/((2.*n+1.)*(2.*n-1.)) );
}

//------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------
void PeriodicCoupling::compute_AlBe_nmlk(Spherical<double> R, complex<double> waveK, int BHreg, int n_max, complex<double> ****AlBe_nmlk){

	// Translation coefficients calculation  -----------------------------------------
	// -------------------------------------------------------------------------------
	// temp variables
	complex<double> c_temp(0., 0.);
	double d_temp(0.), d_temp1(0.), d_temp2(0.);
	// working variables
	int i(0), j(0);
	int ii(0), jj(0);
	int n(0), m(0), l(0), k(0);					// duale subscript variables
	double d_n(0.), d_m(0.), d_l(0.), d_k(0.);	// duale subscript variables
	// auxiliary variables
	complex<double> AlBe2(0., 0.);				// auxiliary coefficients
	complex<double> AlBe3(0., 0.);				// auxiliary coefficients
	complex<double> AlBe4(0., 0.);				// auxiliary coefficients
	complex<double> AlBe2conj(0., 0.);				// auxiliary coefficients
	complex<double> AlBe3conj(0., 0.);				// auxiliary coefficients
	complex<double> AlBe4conj(0., 0.);				// auxiliary coefficients
	double A1(0.), A2(0.), A3(0.), A4(0.);		// auxiliary coefficients
	double B1(0.), B2(0.), B3(0.);				// auxiliary coefficients
	int k_mirror(0), m_mirror(0);				// for calculating (negative)m
	int p=(0); //q(0);							// compound like iterator
	CompoundIterator pl, ql; 					// Create a compound iterator

	// Assign corresponding input values ---------------------------------------------
	complex<double> wkR	= R.rrr*waveK;
	complex<double> wkRconj(0., 0.);
	complex<double> exp_iwk_theji(0., 0.);		// function of m & theji

	// prepare for computing transfer matrices ---------------------------------------
	// Matrices sizes ----------------------------------------------------------------
	// based on n_max
	int n_Matsize, m_Matsize;					// dependent on (n_max)
	n_Matsize = (n_max)+1;						// up to and including n_max		: indexed from 1
	m_Matsize = 2*(n_max)+1;					// up to and including m_max		: indexed from 0 + 1 for m==0
//	int p_max = n_max*n_max + 2*n_max;			// progression relationship - for n, m
	// based on n_max+e							// (e) is required from recurrence relations
	int e=7;									// required for extra values in 'AlBe_00lk'
	int n_Matsize1(0), m_Matsize1(0);			// dependent on (n_max+e)
	n_Matsize1 = (n_max+e)+1;					// up to and including (n_max+e)	: indexed from 1
	m_Matsize1 = 2*(n_max+e)+1;					// up to and including (n_max+e)	: indexed from 0 + 1 for m==0
//	int p_max1 = (n_max+e)*(n_max+e)+2*(n_max+e); // progression relationship - for n, m

	// prepare for calculating translation coefficients -----------------------------
	complex<double> *dataYp;
	dataYp = new complex<double>[pl.max(n_max+e)+1];
	compute_Yp(R, waveK, n_max+e, dataYp);

	// 1.2 Prepare for AlBe_nmlk[0][n_max+e][ii][jj] evaluation ---------------------
	Bessel HB, HBconj;									// create Bessel object
	// waveK
	HB.init(wkR, BHreg, 0, n_Matsize1);					// Initialize Bessel object : start from n==0 - not from n==1!!
	HB.populate();										// calculate Bessel object
	// conj(waveK)
	if(BHreg==0){										// regular case
		wkRconj = conj(wkR);
		HBconj.init(wkRconj, BHreg, 0, n_Matsize1);		// Initialize Bessel object : start from n==0 - not from n==1!!
		HBconj.populate();
	}
	else if(BHreg==1){									// Irregular case
		wkRconj = conj(wkR);
		HBconj.init(consCm1*wkRconj, BHreg, 0, n_Matsize1);					// Initialize Bessel object : start from n==0 - not from n==1!!
		HBconj.populate();
	}

	// Building blocks --------------------------------------------------------------
	// I - fundamental building blocks - Ynm ----------------------------------------
	complex<double> **Ynm;
	Ynm = Tools::Get_2D_c_double(n_Matsize1, m_Matsize1);
	for(i=0; i<n_Matsize1; i++){
		for(j=0; j<m_Matsize1; j++){
			Ynm[i][j]=complex<double>(0., 0.);
		}
	}

	// Obtain Y_nm(the, phi) --------------------------------------------------------
	complex<double> ALegendre00(1., 0.);								// to be changed
	for(i=0, n=0; i<n_Matsize1; i++, n++){								// start at n==0 up to and including n_max
		for(j=0, m=-(n_max+e); j<m_Matsize1; j++, m++){					// increment by the padded e value

			d_n=double(n);	d_m=double(m);

			// if n==m==0 -----------------------------------------------------------
			if(n==0 && m==0){

				d_temp1 = gsl_sf_fact(n-m);								// source of a possible overflow!
				d_temp2 = gsl_sf_fact(n+m);								// source of a possible overflow!
				d_temp  = ((2.*d_n+1.)*d_temp1)/(4*consPi*d_temp2);
				d_temp  = sqrt(d_temp);

				complex<double> exp_ik_phiji(cos(m*R.phi), sin(m*R.phi));
				Ynm[i][j]=d_temp*ALegendre00*exp_ik_phiji;				// eqn (A1)

			}

			// else, for all abs(m)<=n -----------------------------------------------
			else if(abs(m)<=n){

				pl.init(n,m);											// get vector location (p) corresponding to pl(n,m)
				p=pl;
				Ynm[i][j]=dataYp[p];									// eqn (A1) - spherical harmonic
				d_temp*=1.;

				//p++;													// increment p
			}

		}
	}
	// ------------------------------------------------------------------------------

	// II - Global matrix - Alpha-Beta (n, m, l, k) - AlBe_nmlk ---------------------
	complex<double> ****AlBe_nmlkconj;
	AlBe_nmlkconj 	= Tools::Get_4D_c_double(n_Matsize1, m_Matsize1, n_Matsize1, m_Matsize1);
	for(i=0; i<n_Matsize1; i++){
		for(j=0; j<m_Matsize1; j++){
			for(ii=0; ii<n_Matsize1; ii++){
				for(jj=0; jj<m_Matsize1; jj++){
					AlBe_nmlk[i][j][ii][jj]		= complex<double>(0., 0.);
					AlBe_nmlkconj[i][j][ii][jj]	= complex<double>(0., 0.);
				}
			}
		}
	}

	// A - calculate AlBe_00lk - fundamental building block -------------------------
	n=0;			m=0;						// n==m==0
	d_n=double(n);	d_m=double(m);				// n==m==0
	for(ii=0, l=0; ii<n_Matsize1; ii++, l++){
		for(jj=0, k=-(n_max+e); jj<m_Matsize1; jj++, k++){

			d_l=double(l);	d_k=double(k);

			// get mirror image in k
			k_mirror = n_max+e + (n_max+e-jj);

			if(abs(k)<=l){

				// AlBe_00lk[ii][jj] : Bessel/Hankel - Legendre - exp(imthe)
				d_temp=sqrt(4.*consPi)*pow(-1.,d_l+d_k);
				// wavek
				AlBe_nmlk[0][n_max+e][ii][jj]		= d_temp*Ynm[ii][k_mirror]*HB.data[l];		// eqn (C3)
				// waveconj
				AlBe_nmlkconj[0][n_max+e][ii][jj]	= d_temp*Ynm[ii][k_mirror]*HBconj.data[l];	// eqn (C3)

			}	// if(abs(k)<=l)

		}	// jj
	}	// ii

	// B.1 - calculate AlBe_nmlk - positive(m)==n ----------------------------------
	for (i = 1, n = 1; i < n_Matsize1 - e; i++, n++) { // start at n==1 up to and including n_max
		for (j = e, m = -n_max; j < m_Matsize1 - e; j++, m++) { // increment by the padded e value

			d_n = double(n);
			d_m = double(m);

			// if n==positive(m) ---------------------------------------------------
			if (n == m) {

				for (ii = 0, l = 0; ii < n_Matsize1 - 1; ii++, l++) { // start at n==1 up to and including n_max+e
					for (jj = 0, k = -(n_max + e - 0); jj < m_Matsize1 - 0; jj++, k++) { // increment by the padded e value

						d_l = double(l);
						d_k = double(k);

						if (abs(k) <= l) {
							// obtain three coefficients ----------------------------
							B1 = b_nm_p(d_n - 1., d_n - 1.);
							B2 = b_nm_p(d_l - 1., d_k - 1.);
							B3 = b_nm_m(d_l + 1., d_k - 1.);

							if (((ii - 1) < 0) || ((jj - 1) < 0)){
								AlBe2 		= complex<double>(0., 0.);
								AlBe2conj 	= complex<double>(0., 0.);
							}
							else{
								AlBe2 		= AlBe_nmlk[i - 1][j - 1][ii - 1][jj - 1];
								AlBe2conj 	= AlBe_nmlkconj[i - 1][j - 1][ii - 1][jj - 1];
							}
							//
							if (((ii + 1) > n_Matsize1-1) || ((jj - 1) < 0)){
								AlBe3 		= complex<double>(0., 0.);
								AlBe3conj 	= complex<double>(0., 0.);
							}
							else{
								AlBe3 		= AlBe_nmlk[i - 1][j - 1][ii + 1][jj - 1];
								AlBe3conj 	= AlBe_nmlkconj[i - 1][j - 1][ii + 1][jj - 1];
							}

							// evaluate current term --------------------------------
							// wavek
							AlBe_nmlk[i][j][ii][jj] = B2 * AlBe2
													+ B3 * AlBe3;
							AlBe_nmlk[i][j][ii][jj] /= B1;						// eqn (C1-B)
							// wavekconj
							AlBe_nmlkconj[i][j][ii][jj] = B2 * AlBe2conj
														+ B3 * AlBe3conj;
							AlBe_nmlkconj[i][j][ii][jj] /= B1;					// eqn (C1-B)
						}
					} // jj
				} // ii

			} // if(n>0 && n==m)

		} // j
	} // i

	// B.2 - calculate AlBe_nmlk - for all other positive(m) -----------------------
	for (i = 1, n = 1; i < n_Matsize1 - e; i++, n++) { // start at n==1 up to and including n_max
		for (j = e, m = -n_max; j < m_Matsize1 - e; j++, m++) { // increment by the padded e value

			d_n = double(n);
			d_m = double(m);

			if (abs(m) <= n) {

				// for all other positive m ------------------------------------------------
				if (m >= 0 && m != n) {

					for (ii = 0, l = 0; ii < n_Matsize1 - 1; ii++, l++) { // start at n==1 up to and including n_max+e
						for (jj = 0, k = -(n_max + e - 0); jj < m_Matsize1 - 0; jj++, k++) { // increment by the padded e value

							d_l = double(l);
							d_k = double(k);

							if (abs(k) <= l) {
								// obtain three coefficients -------------------------------
								A1 = a_nm_p(d_n - 1., d_m);
								A2 = a_nm_m(d_n - 1., d_m);
								A3 = a_nm_p(d_l - 1., d_k);
								A4 = a_nm_m(d_l + 1., d_k);

								// evaluate current term -----------------------------------
								if (((i - 2) < 0)){
									AlBe2 = complex<double> (0., 0.);
									AlBe2conj = complex<double> (0., 0.);
								}
								else{
									AlBe2 		= AlBe_nmlk[i - 2][j][ii][jj];
									AlBe2conj 	= AlBe_nmlkconj[i - 2][j][ii][jj];
								}
								//
								if (((ii - 1) < 0)){
									AlBe3 = 0.;
									AlBe3conj = 0.;
								}
								else{
									AlBe3 		= AlBe_nmlk[i - 1][j][ii - 1][jj];
									AlBe3conj 	= AlBe_nmlkconj[i - 1][j][ii - 1][jj];
								}
								//
								if (((ii + 1) > n_Matsize1 - 1)){
									AlBe4 = 0.;
									AlBe4conj = 0.;
								}
								else{
									AlBe4 		= AlBe_nmlk[i - 1][j][ii + 1][jj];
									AlBe4conj 	= AlBe_nmlkconj[i - 1][j][ii + 1][jj];
								}

								// wavek
								AlBe_nmlk[i][j][ii][jj] = -A2 * AlBe2
														+  A3 * AlBe3
														+  A4 * AlBe4;
								AlBe_nmlk[i][j][ii][jj] /= A1;					// eqn (C1-A)
								// wavek
								AlBe_nmlkconj[i][j][ii][jj] = -A2 * AlBe2conj
															+  A3 * AlBe3conj
															+  A4 * AlBe4conj;
								AlBe_nmlkconj[i][j][ii][jj] /= A1;					// eqn (C1-A)
							}
						} // jj
					} // ii

				} // if(m>=0 && m!=n)

			} // if(m<=n)

		} // j
	} // i

	// C - calculate AlBe_nmlk - for all other negative(m) --------------------------------
	for (i = 1, n = 1; i < n_Matsize1 - e; i++, n++) { // start at n==1 up to and including n_max
		for (j = e, m = -n_max; j < m_Matsize1 - e; j++, m++) { // increment by the padded e value

			d_n = double(n);
			d_m = double(m);

			if (abs(m) <= n) {

				// for all other negative m, including n==abs(m) ---------------------------
				if (m < 0) {

					for (ii = 0, l = 0; ii < n_Matsize1 - 1; ii++, l++) { // start at n==1 up to and including n_max+e
						for (jj = 0, k = -(n_max + e - 0); jj < m_Matsize1 - 0; jj++, k++) { // increment by the padded e value

							d_l = double(l);
							d_k = double(k);

							if (abs(k) <= l) {

								// get mirror image in m -----------------------------------
								m_mirror = n_max + e + (n_max + e - j);
								// get mirror image in k
								k_mirror = n_max + e + (n_max + e - jj);

								if(BHreg==0){															// regular 		- Bessel
									// obtain mirror coefficient -------------------------------
									d_temp = pow(-1., d_k + d_m);										// eqn (C4-A)
									AlBe_nmlk[i][j][ii][jj] = d_temp * conj(AlBe_nmlkconj[i][m_mirror][ii][k_mirror]);
								}
								else if(BHreg==1){														// iregular 	- Hankel
									// obtain mirror coefficient -------------------------------
									d_temp = pow(-1., d_n + d_l + d_k + d_m);										// eqn (C4-B)
									AlBe_nmlk[i][j][ii][jj] = 1.*d_temp * conj(AlBe_nmlkconj[i][m_mirror][ii][k_mirror]);
								}

							}
						} // jj
					} // ii

				} // if(m<0)

			} // if(m<=n)

		} // j
	} // i

	// C - delete all auxiliary matrices ---------------------------------------------------
	for(i=0; i<n_Matsize1; i++){
		for(j=0; j<m_Matsize1; j++){
			for(ii=0; ii<n_Matsize1; ii++){
				delete [] AlBe_nmlkconj[i][j][ii];
	}	}	}

	for(i=0; i<n_Matsize1; i++){
		for(j=0; j<m_Matsize1; j++){
				delete [] AlBe_nmlkconj[i][j];
	}	}

	for(i=0; i<n_Matsize1; i++){
				delete [] AlBe_nmlkconj[i];
				delete [] Ynm[i];
	}

	delete [] AlBe_nmlkconj;
	delete [] Ynm;
	delete [] dataYp;			// AJ - this needs to be deleted - currently producing an error - AJ //

}
// ---------------------------------------------------------------------------------------



// ---------------------------------------------------------------------------------------
// compute Bnm function - obtained with the aid of Wigner 3j symboles
complex<double> PeriodicCoupling::compute_Bnm(int n, int n1, int n2, int m, int m1, int m2){

	return sqrt( ((2*n+1) * (2*n1+1) * (2*n2+1)) / (4*consPi) )*
			Symbol::Wigner3j(n, n1, n2,  0,  0,  0)*
			Symbol::Wigner3j(n, n1, n2,  m,  m1, m2);

}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// compute alpha_nm function - MODINOS 1987 eq (6a)
double PeriodicCoupling::compute_alpha_nm(int n, int m){

	return 0.5*sqrt(double((n-m)*(n+m+1)));

}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// compute beta_nm function - MODINOS 1987 eq (6b)
double PeriodicCoupling::compute_beta_nm(int n, int m){

	return 0.5*sqrt(double((n+m)*(n-m+1)));

}
// ---------------------------------------------------------------------------------------

// ---------------------------------------------------------------------------------------
// compute Znmlk(-Rn) function - this is obtained from the scalar addition-translation theorem
// Znmlk(-Rn) = G_nmlk *sum(exp(iK||.Rn))	MODINOS 1987 eq (6b)
// This is computed for the special case of a rectangular lattice only for now
// a1 = x-axis unit size of the 2D lattice
// a2 = y-axis unit size of the 2D lattice
// Rmax is a truncation limit of on the reciprocal lattice - user defined for now
int PeriodicCoupling::compute_Znmlk(double a1, double a2, complex<double> waveK, int Rmax, int n_max, complex<double> ****Z_nmlk){

	int ll(0);
	int BHreg(1);				// indicates Hankel / non-regular

	double b1(0.), b2(0.);
//	double Rn_mag(0.);
	Spherical<double> Rs;
	Cartesian<double> Rn;
	Cartesian<double> waveKc;	// conversion of vKinc into Cartesian coords wrt to the 2D lattice - make use of components that are parallel to the 2D lattice only
	double dotproduct2D(0.);
	complex<double> exp_iKpRn_p(0., 0.);
	complex<double> exp_iKpRn_n(0., 0.);

	// Obtain Cartesian values of spherical vKinc
	// start: AJ - this is temporary - AJ //
	Spherical<double> vKinc;
	double lam =1.;
	vKinc = Spherical<double>(2*consPi / lam, (90/180.)*consPi, (90/180.)*consPi);
	// end:   AJ - this is temporary - AJ //
	waveKc = Tools::toCartesian(vKinc);		// (waveKp_x, waveKp_y, waveKp_z) == (K||, Kz)

	// declare variable and matrices for computing Gnmlk==AlBe_nmlk ------------------------
	// Transfer function AlBe --------------------------------------------------------------
	int i, j, n, m, ii, jj, l, k;
	int e=1;									// required for extra values in 'AlBe_00lk'
	int n_Matsize1(0), m_Matsize1(0);			// dependent on (n_max+e)
	n_Matsize1 = (n_max+e)+1;					// up to and including (n_max+e)	: indexed from 1
	m_Matsize1 = 2*(n_max+e)+1;					// up to and including (n_max+e)	: indexed from 0 + 1 for m==0
	complex<double> ****Gp_nmlk, ****Gn_nmlk; 	// positive Rn, negative Rn
	Gp_nmlk = Tools::Get_4D_c_double(n_Matsize1, m_Matsize1, n_Matsize1, m_Matsize1);	// positive Rn
	Gn_nmlk = Tools::Get_4D_c_double(n_Matsize1, m_Matsize1, n_Matsize1, m_Matsize1);	// negative Rn
	for (i = 0; i < n_Matsize1; i++) {
		for (j = 0; j < m_Matsize1; j++) {
			for (ii = 0; ii < n_Matsize1; ii++) {
				for (jj = 0; jj < m_Matsize1; jj++) {
						Gp_nmlk[i][j][ii][jj] =complex<double>(0., 0.);
						Gn_nmlk[i][j][ii][jj] =complex<double>(0., 0.);
	} } } }

	// Determine reciprocal lattice constants -----------------------------------------------------
	// Assuming cubic lattice for now
	// start: AJ - this is temporary - AJ //
	b1 = 2.*consPi/a1;			// rectangular lattice only
	b2 = 2.*consPi/a2;			// rectangular lattice only
	// end:   AJ - this is temporary - AJ //

	// By definition: loop over all lattice constants - excluding the self-term Rn=0	MODINOS 1987 eq (6b)
	for(ll=1; ll<=Rmax; ll++){

		Rn.x 			= double(ll)*b1;
		Rn.y 			= double(ll)*b2;
		Rn.z 			= 0.;

//		// get magnitude of Rn
//		Rn_mag = sqrt( Rn.x*Rn.x + Rn.y*Rn.y + Rn.z*Rn.z );
//		// Automatic check if truncation limit is reached - for later - ------------------------
//		if(Rn_mag - Rmax <0.){

		// positive lattice members of Rn(+l) --------------------------------------------------
		Rn.x 			= double(ll)*b1;
		Rn.y 			= double(ll)*b2;
		Rn.z 			= 0.;
		dotproduct2D 	= waveKc.x*Rn.x + waveKc.y*Rn.y;
		exp_iKpRn_p		= complex<double> ( cos(dotproduct2D), sin(dotproduct2D) );
		// obtain corresponding Gnmlk(-Rn)==AlBe_nmlk(-Rn)
		// negate relative vector rn=r-Rn, r==Origin -> rn=-Rn
		Rn.x*=-1;		Rn.y*=-1;		Rn.z*=-1;
		BHreg = 1;		// Hankel / non-regular
		Rs = Tools::toSpherical(Rn);		// AJ - to be double checked - AJ //
		compute_AlBe_nmlk(Rs, waveK, BHreg, n_max, Gp_nmlk);
		// -------------------------------------------------------------------------------------

		// positive lattice members of Rn(+l) --------------------------------------------------
		Rn.x 			= double(-ll)*b1;
		Rn.y 			= double(-ll)*b2;
		Rn.z 			= 0.;
		dotproduct2D	= waveKc.x*Rn.x + waveKc.y*Rn.y;
		exp_iKpRn_n		= complex<double> ( cos(dotproduct2D), sin(dotproduct2D) );
		// obtain corresponding Gnmlk(-Rn)==AlBe_nmlk(-Rn)
		// negate relative vector rn=r-Rn, r==Origin -> rn=-Rn
		Rn.x*=-1;		Rn.y*=-1;		Rn.z*=-1;
		BHreg = 1;		// Hankel / non-regular
		Rs = Tools::toSpherical(Rn);		// AJ - to be double checked - AJ //
		compute_AlBe_nmlk(Rs, waveK, BHreg, n_max, Gn_nmlk);
		// ------------------------------------------------------------------------------------

		// update Z_nmlk:	Z_nmlk+= G_nmlk(-Rn) *sum(exp(iK||.Rn)) -----------------------
		for (i = 1, n = 1; i < n_Matsize1 - e; i++, n++) { // start at n==1 up to and including n_max
			for (j = e, m = -n_max; j < m_Matsize1 - e; j++, m++) { // increment by the padded e value
				if (abs(m) <= n) {
					for (ii = 0, l = 0; ii < n_Matsize1 - 1; ii++, l++) { // start at n==1 up to and including n_max+e
						for (jj = 0, k = -(n_max + e - 0); jj < m_Matsize1 - 0; jj++, k++) { // increment by the padded e value
							if (abs(k) <= l) {
								// update Z_nmlk
								// from positive integers Rn(+l)
								Z_nmlk[i][j][ii][jj] += exp_iKpRn_n*Gp_nmlk[i][j][ii][jj];
								// from positive integers Rn(-l)
								Z_nmlk[i][j][ii][jj] += exp_iKpRn_n*Gp_nmlk[i][j][ii][jj];
							} // if(abs(k)<=l)
						} // jj
					} // ii
				} // if(abs(m)<=n)
			} // j
		} // i
		// --------------------------------------------------------------------------------------
//		}
//		// otherwise, exit routine --------------------------------------
//		else{
//			break;
//		}

	}

	// delete all other matrices ----------------------------------------
	for(i=0; i<n_Matsize1; i++){
		for(j=0; j<m_Matsize1; j++){
			for(ii=0; ii<n_Matsize1; ii++){
				delete [] Gp_nmlk[i][j][ii];
				delete [] Gn_nmlk[i][j][ii];
	}	}	}

	for(i=0; i<n_Matsize1; i++){
		for(j=0; j<m_Matsize1; j++){
				delete [] Gp_nmlk[i][j];
				delete [] Gn_nmlk[i][j];
	}	}

	for(i=0; i<n_Matsize1; i++){
				delete [] Gp_nmlk[i];
				delete [] Gn_nmlk[i];
	}

	delete [] Gp_nmlk;
	delete [] Gn_nmlk;

	return 0;

}
// ---------------------------------------------------------------------------------------

//// ---------------------------------------------------------------------------------------
// compute Omega_nmlk function - this is identically the scalar addition-translation theorem
int PeriodicCoupling::compute_OMEGAnmlk(complex<double> waveK, int n_max, complex<double> **dataOMEGA_1pq, complex<double> **dataOMEGA_2pq){

	// first of all increment n_max by 1 to allow for all required values of Z_nmlk
	n_max+=1;

	// Other quantities
	double alpha_nm_I(0.), alpha_nm_II(0.);
	double beta_nm_I(0.), beta_nm_II(0.);
	complex<double> B_nm_I(0., 0.), B_nm_II(0., 0.);
	complex<double> c_temp_I(0., 0.), c_temp_II(0., 0.), c_temp_III(0., 0.);

	// 1- compute Z_nmlk -------------------------------------------------------------------
	int i, j, n, m, ii, jj, l, k;
	int e=1;									// required for extra values in 'AlBe_00lk'
	int n_Matsize1(0), m_Matsize1(0);			// dependent on (n_max+e)
	n_Matsize1 = (n_max+e)+1;					// up to and including (n_max+e)	: indexed from 1
	m_Matsize1 = 2*(n_max+e)+1;					// up to and including (n_max+e)	: indexed from 0 + 1 for m==0
	complex<double> ****Z_nmlk, ****OMEGA_1_nmlk, ****OMEGA_2_nmlk;
	Z_nmlk = Tools::Get_4D_c_double(n_Matsize1, m_Matsize1, n_Matsize1, m_Matsize1);
	for (i = 0; i < n_Matsize1; i++) {
		for (j = 0; j < m_Matsize1; j++) {
			for (ii = 0; ii < n_Matsize1; ii++) {
				for (jj = 0; jj < m_Matsize1; jj++) {
						Z_nmlk[i][j][ii][jj] = complex<double>(0., 0.);
						OMEGA_1_nmlk[i][j][ii][jj] = complex<double>(0., 0.);
						OMEGA_2_nmlk[i][j][ii][jj] = complex<double>(0., 0.);
	} } } }

	// 2- obtain Z_ nmlk -------------------------------------------------------------------
	// This is computed for the special case of a rectangular lattice only for now
	double a1 = 1.;	// x-axis unit size of the 2D lattice
	double a2 = 1.;	// y-axis unit size of the 2D lattice
	int Rmax = 100;	// Rmax is a truncation limit of on the reciprocal lattice - user defined for now
	compute_Znmlk(a1, a2, waveK, Rmax, n_max, Z_nmlk);

	// 3- obtain OMEGA_1_nmlk & OMEGA_2_nmlk from Z_nmlk -----------------------------------
	for (i = 1, n = 1; i < n_Matsize1 - e; i++, n++) { // start at n==1 up to and including n_max
		for (j = e, m = -n_max; j < m_Matsize1 - e; j++, m++) { // increment by the padded e value
			if (abs(m) <= n) {
				for (ii = 1, l = 1; ii < n_Matsize1 - 1; ii++, l++) { // start at n==1 up to and including n_max+e
					for (jj = e, k = -n_max; jj < m_Matsize1 - 0; jj++, k++) { // increment by the padded e value
						if (abs(k) <= l) {

							// A- obtain OMEGA_1_nmlk -------------------------------------
							// I- evaluate all other quantities ---------------------------
							alpha_nm_I = compute_alpha_nm(l, k-1);
							beta_nm_I = compute_alpha_nm(n, m);
							// check if any of the above are ZERO by default!
							if(Tools::equalDoubles(alpha_nm_I, 0.) || Tools::equalDoubles(beta_nm_I, 0.))
								c_temp_I = complex<double>(0., 0.);
							else
								c_temp_I = 2.*alpha_nm_I*beta_nm_I*Z_nmlk[i][j-1][ii][jj-1];
							// II- evaluate all other quantities
							alpha_nm_II = compute_alpha_nm(n, m);
							beta_nm_II = compute_alpha_nm(l, k+1);
							// check if any of the above are ZERO by default!
							if(Tools::equalDoubles(alpha_nm_II, 0.) || Tools::equalDoubles(beta_nm_II, 0.))
								c_temp_II = complex<double>(0., 0.);
							else
								c_temp_II = 2.*alpha_nm_II*beta_nm_II*Z_nmlk[i][j+1][ii][jj+1];
							// III- evaluate all other quantities
								c_temp_III = double(m)*double(k)*Z_nmlk[i][j][ii][jj];
							// obtain OMEGA_1_nmlk --------------------------------------
							OMEGA_1_nmlk[i][j][ii][jj] = 	1./sqrt(double(n*(n+1)*l*(l+1)))
															*(c_temp_I + c_temp_II + c_temp_III);
							// ----------------------------------------------------------

							// B- obtain OMEGA_2_nmlk -----------------------------------
							// I- evaluate all other quantities -------------------------
							alpha_nm_I = compute_alpha_nm(n, m);
							B_nm_I = compute_Bnm(l-1, k+1, 1, -1, l, k);
							// check if any of the above are ZERO by default!
							if(Tools::equalDoubles(alpha_nm_I, 0.) || (Tools::equalDoubles(B_nm_I.real(), 0.) && (Tools::equalDoubles(B_nm_I.imag(), 0.))))
								c_temp_I = complex<double>(0., 0.);
							else
								c_temp_I = sqrt(8.*consPi/3.)*pow(-1.,k)*alpha_nm_I*B_nm_I*Z_nmlk[i][j+1][ii-1][jj+1];
							// II- evaluate all other quantities
							beta_nm_II = compute_alpha_nm(n, m);
							B_nm_I = compute_Bnm(l-1, k-1, 1, 1, l, k);
							// check if any of the above are ZERO by default!
							if(Tools::equalDoubles(beta_nm_II, 0.) || (Tools::equalDoubles(B_nm_II.real(), 0.) && (Tools::equalDoubles(B_nm_II.imag(), 0.))))
								c_temp_II = complex<double>(0., 0.);
							else
								c_temp_I =-sqrt(8.*consPi/3.)*pow(-1.,k)*beta_nm_I*B_nm_I*Z_nmlk[i][j-1][ii-1][jj-1];
							// III- evaluate all other quantities
								c_temp_III = double(m)*Z_nmlk[i][j][ii-1][jj]*sqrt((l+k)*(l-k)/(2.*l-1)/(2.*l+1));
							// obtain OMEGA_1_nmlk --------------------------------------
							OMEGA_2_nmlk[i][j][ii][jj] = 	1./sqrt(double(n*(n+1)*l*(l+1)))
															*(2.*l+1)/waveK
															*(c_temp_I + c_temp_II + c_temp_III);
							// ----------------------------------------------------------

						} // if(abs(k)<=l)
					} // jj
				} // ii
			} // if(abs(m)<=n)
		} // j
	} // i

	// 4 - ---------------------------------------------------------------------------------
	// store OMEGA_1_nmlk & OMEGA_2_nmlk into a compound iterator like 2D matrix ---------------------
	int n_check(0), m_check(0), l_check(0), k_check(0);
	CompoundIterator pl, ql; 					// Create a compound iterator
	int p, q;
	for(i=1, n=1; i<n_Matsize1-e; i++, n++){							// start at n==1 up to and including n_max
		for(j=e, m=-n_max; j<m_Matsize1-e; j++, m++){					// increment by the padded e value

			if(abs(m)<=n){

				// double check current n & m -----------------------------------------------
				pl.init(n,m);
				p=pl;
				n_check = int(sqrt(double(p)+1.));
				m_check = -(p+1)+n*(n+1);
				//q=0;

				// for all m values --------------------------------------------------------
				for(ii=1, l=1; ii<n_Matsize1-e; ii++, l++){
					for(jj=e, k=-n_max; jj<m_Matsize1-e; jj++, k++){

						if(abs(k)<=l){
							// double check current l & k -----------------------------------
							ql.init(l,k);
							q=ql;
							l_check = int(sqrt(double(q)+1.));
							k_check = -(q+1)+l*(l+1);
							// map 4D matrix --> 2D matrix ----------------------------------
							dataOMEGA_1pq[p][q]=OMEGA_1_nmlk[i][j][ii][jj];
							dataOMEGA_2pq[p][q]=OMEGA_2_nmlk[i][j][ii][jj];
							//q++;										// increment q
						}

					}	// jj
				}	// ii

				//p++;													// increment p

			}	// if(abs(m)<=n)

		}	// j
	}	// i

	// C - delete all other matrices --------------------------------------------------------
	for(i=0; i<n_Matsize1; i++){
			for(j=0; j<m_Matsize1; j++){
				for(ii=0; ii<n_Matsize1; ii++){
					delete [] Z_nmlk[i][j][ii];
					delete [] OMEGA_1_nmlk[i][j][ii];
					delete [] OMEGA_2_nmlk[i][j][ii];
		}	}	}

		for(i=0; i<n_Matsize1; i++){
			for(j=0; j<m_Matsize1; j++){
					delete [] Z_nmlk[i][j];
					delete [] OMEGA_1_nmlk[i][j];
					delete [] OMEGA_2_nmlk[i][j];
		}	}

		for(i=0; i<n_Matsize1; i++){
					delete [] Z_nmlk[i];
					delete [] OMEGA_1_nmlk[i];
					delete [] OMEGA_2_nmlk[i];
		}

		delete [] Z_nmlk;
		delete [] OMEGA_1_nmlk;
		delete [] OMEGA_2_nmlk;

	return 0;

}
// ---------------------------------------------------------------------------------------

int PeriodicCoupling::compute_Ap_vkg(int n_max, Cartesian<double> vecKg, complex<double> waveK, double Ao, Cartesian<complex<double> > *dataAp_vkg){

	int n, m;
	CompoundIterator p, i, j;

	// Other quantities
	Spherical<double> R;
	double alpha_nm(0.), beta_nm(0.);
	complex<double> c_temp(0., 0.), c_temp_I(0., 0.), c_temp_II(0., 0.), c_temp_III(0., 0.);

	// compute spherical harmonics for n_max+1
	complex<double> *dataYp_R;
	dataYp_R = new complex<double> [p.max(n_max+1)+1];		// extra element for (n,m)=(0,0)

	R = Tools::toSpherical(vecKg);							// AJ - to be double checked - AJ //
	compute_Yp(R, waveK, n_max+1, dataYp_R);

	// obtain dataAp_vKg --------------------------------------------------------------
	for(p=0; p<p.max(n_max); p++){

		// Aux values
		n = abs(p.first);
		m = abs(p.second);
		alpha_nm 	= compute_alpha_nm(n, m);
		beta_nm  	= compute_alpha_nm(n, m);
		i.init(n, m+1);
		j.init(n, m-1);
		c_temp = 2.*consPi*pow(complex<double>(0, 1),-n)
				/(waveK*Ao*vecKg.z*sqrt(n*(n+1)));

		// I --------------------------------------------
		c_temp_I		= alpha_nm*dataYp_R[i] + beta_nm*dataYp_R[j];

		// II -------------------------------------------
		c_temp_II 		= alpha_nm*dataYp_R[i] - beta_nm*dataYp_R[j];
		c_temp_II 		*=complex<double>(0, -1);

		// III ------------------------------------------
		c_temp_III 		= double(p.second)*dataYp_R[p];

		// allocate values of dataAp_vkg from I, II, III
		dataAp_vkg[p].x	= c_temp*c_temp_I;
		dataAp_vkg[p].y	= c_temp*c_temp_II;
		dataAp_vkg[p].z	= c_temp*c_temp_III;
	}


//	return(		complex<double>(1,  0)*(compute_alpha_nm(n, m) + compute_beta_nm(n, m)),
//				complex<double>(0, -1)*(compute_alpha_nm(n, m) - compute_beta_nm(n, m)),
//				complex<double>(1,  0)*m);
	return 0;
}
