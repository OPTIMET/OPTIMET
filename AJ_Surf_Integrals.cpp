/*
 * AJ_Surf_Integrals.cpp
 *
 *  Created on: Mar 4, 2013
 *      Author: A Al-Jarro
 */

#include <iostream>
#include <math.h>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <complex>
#include "Spherical.h"
//#include "AuxCoefficients.h"
#include "Bessel.h"
//#include "Reader.h"
#include "CompoundIterator.h"
#include "Legendre.h"
#include "Excitation.h"
#include "Coupling.h"
#include "SphericalP.h"
#include "Solver.h"

//#include "AJ_Wigner.h"
//#include "AJ_AuxFuns.h"
//#include "AJ_SphWaves.h"
//#include "AJ_TransCoeff.h"
#include "gsl/gsl_sf_legendre.h"
//#include "AJ_Y00Trans.h"
#include "AJ_Surf_Integrals.h"

#define PI 3.14159265358979323846
#define EPS 3.0e-11 							// EPS is the relative precision.

int VIGdVIG(int nMax, int m_, Spherical<double> R, double *Wigner, double *dWigner){
        /*-------------------------------------------------------------------------------*/
                /* PURPOSE: Evaluate Wigner d function & its derivative : vig_d  and d_vig_d --- */
                /* VIG_d(n)             = vig_d  (l=0, m, n, the) -------------------------------------- */
                /* d_VIG_d(n)   = d_vig_d(l=0, m, n, the) -------------------------------------- */
                /*-------------------------------------------------------------------------------*/

                // variables ---------------------------------------------------------------------
                // loop variables ----------------------------------------------------------------
                int i(0), ii(0);
                int check_m_negative(0);                                                // flag for m<0
                // temp variables ---------------------------------------------------------------
                double d_temp(0.);                                                              // double temp
                double d_temp1(0.), d_temp2(0.), d_temp3(0.), d_temp4(0.), d_temp5(0.), d_temp6(0.);
                // Wigner function indices -------------------------------------------------------
                int n_min(0);                                                                   // vector spherical function index              : n>=1, n_min=max(|l|,|m|)
                int l(0), l_abs(0);                                                             // wigner function index
                int m_abs(0);                                                                   // vector spherical function index
                // Wigner function variables -----------------------------------------------------
                //double vig_d(0.);                                                             // wigner function value
                double vig_the=R.the;                                                           // wigner function argument                             : [0<=the<=PI]
                double vig_x(0.);                                                               // wigner function auxiliary variable   : x=cos(the)
                double vig_exy_lm(0.);                                                  // wigner function auxiliary variable   : exy_lm=eq(B.16)                       - book reference
                // Wigner recursive relationship terms - to be used for finally obtaining vig_d -
                double vig_d_n_min(0.);                                                 // wigner function recursive term               : d_n_min=eq(B.24)                      - book reference
                // -------------------------------------------------------------------------------

                // ------------------ prepare for VIG functions evaluations ----------------------
                // -------------------------------------------------------------------------------
                for(i=0; i<=nMax; i++){                                                 // rest
                        Wigner[i]=0.;
                        dWigner[i]=0.;
                }
                // prepare for vig_d evaluation --------------------------------------------------
                l_abs = abs(l);
                m_abs = abs(m_);
                // -------------------------------------------------------------------------------

                // -------------------------------------------------------------------------------
                // I - determine n_min : n_min=max(|l|,|m|) --------------------------------------
        //      if(l_abs>n_min)
          //    n_min = l_abs;
        //      if(m_abs>n_min)
                        n_min = m_abs;                                                          // this is always true, since l==0
                // -------------------------------------------------------------------------------


                // -------------------------------------------------------------------------------
                // II - obtain vig_d_n_min ---------------------------------------------------
                // 0 - set : vig_the==vig_the
        //      if(m_>=0)                                                                       // solve directly eqs(B.22-B24)
        //              vig_the = vig_the;                                              // keep as it is
                if(m_<0){                                                                       // solve using symmetry relation eq(B.7)
                        vig_the = consPi - vig_the;                             // modify value of vig_x   : eq(B.7)
                        m_=abs(m_);                                                             // obtain solution for |m| : eq(B.7)
                        check_m_negative=1;                                             // check m<0 : to use in IV below
                }
                vig_x = cos(vig_the);                                           // calculate in radians

                // A - evaluate vig_exy_lm                                      // wigner function recursive term in d_n_min=eq(B.23)   : exy_lm=eq(B.16)       - book reference
                if(m_>=l)
                        vig_exy_lm      = 1.;                                           // this is always true, since l==0 & m>=0
                else if(m_<l)
                        vig_exy_lm      = pow(-1., double(l-m_));

                // B - evaluate all other terms in d_n_min = eqn(B.24) -----------------------
                d_temp1 = pow(2., -n_min);
                //
                d_temp2 = sqrt(double(gsl_sf_fact(2*n_min)));
                d_temp3 = sqrt(double(gsl_sf_fact(abs(l-m_))));
                d_temp4 = sqrt(double(gsl_sf_fact(abs(l+m_))));
                //
                d_temp5 = pow(1. - vig_x, abs(double(l-m_))/2.);
                d_temp6 = pow(1. + vig_x, abs(double(l+m_))/2.);

                // C - evaluate vig_d_n_min --------------------------------------------------
                // vig_d_n_min = vig_exy_lm*d_temp1*(d_temp2/d_temp3/d_temp4)*d_temp5*d_temp6;
                vig_d_n_min = vig_exy_lm;
                vig_d_n_min*= d_temp1;
                vig_d_n_min*= d_temp2;
                vig_d_n_min/= d_temp3;
                vig_d_n_min/= d_temp4;
                vig_d_n_min*= d_temp5;
                vig_d_n_min*= d_temp6;
                // ----------------------------------------------------------------------------

                // ----------------------------------------------------------------------------
                // III - calculate VIG_d & d_VIG_d --------------------------------------------
                // III.1 - check for singularity cases - i.e (m==0) ---------------------------
                if(m_==0){                                                                      // singularity check
                        // A - if (i<n_min), then set VIG_d to zero -----------------------------------
        //              for(i=0; i<n_min; i++){                                 // for loop
                //              VIG_d[i]  =0.;                                          // == zero by default
                        //}                                                                             // end for loop
                        // B - obtain the readily availble VIG_d[n_min] value ------------------------
                        // vig_d ---------------------------------------------------------------------
                        Wigner[n_min] = vig_d_n_min;                                    // by definition
                        // C - obtain all other values in VIG_d[i] from recursive relationship --------
                        //   - Based on 'n_min' value and 'n' value; total recursive steps == n - n_min
                        // obtain first term in recursive relationship eq(B.22) - from special case
                        Wigner[n_min+1] = vig_x;
                        // d_vig_d for the previous from current : d_VIG_d[i] = 0.*VIG_d[i-1] + 0.*VIG_d[i] + ...*VIG_d[i+1]    - eqn(B.26)
                        ii=n_min;
                        d_temp1 = sqrt((double(ii)+1.)*(double(ii)+1.) - double(l*l));
                        d_temp2 = sqrt((double(ii)+1.)*(double(ii)+1.) - double(m_*m_));
                        d_temp3 = (2.*double(ii)+1.);
                        d_temp4 = (double(ii)+1.);
                        dWigner[ii] = double(ii)*d_temp1*d_temp2/(d_temp3*d_temp4)*Wigner[ii+1];
                        dWigner[ii]/=sin(vig_the);                              // == zero by definition
                        // obtain all successive terms in eq(B.22) - directrly
                        for(i=n_min+2; i<=nMax+1; i++){                 // for loop
                                ii=i-1;
                                d_temp1 = sqrt((double(ii)+1.)*(double(ii)+1.) - double(l*l));
                                d_temp2 = sqrt((double(ii)+1.)*(double(ii)+1.) - double(m_*m_));
                                d_temp3 = (2.*double(ii)+1.)*( ii*(ii+1.)*vig_x - l*m_ );
                                d_temp4 = (double(ii)+1.);
                                d_temp5 = sqrt(double(ii*ii) - double(l*l));
                                d_temp6 = sqrt(double(ii*ii) - double(m_*m_));
                                // evaluate from above
                                // vig_d        : VIG_d[i] = ...*VIG_d[i-1] + ...*VIG_d[i-2]                                            - eqn(B.22)
                                d_temp  = (1./double(ii)/d_temp1/d_temp2)*( (d_temp3)*Wigner[i-1]
                                                - (d_temp4*d_temp5*d_temp6)*Wigner[i-2] );
                                if(i<=nMax)
                                        Wigner[i] = d_temp;
                                // d_vig_d      : d_VIG_d[i-1] = ...*VIG_d[i-2] + ...*VIG_d[i-1] + ...*VIG_d[i]         - eqn(B.26)
                                dWigner[i-1] =  - d_temp4*d_temp5*d_temp6/(double(ii)*(2.*double(ii)+1.))*Wigner[i-2]
                                                                - 0.                            // since l==0
                                                                + double(ii)*d_temp1*d_temp2/(d_temp4*(2.*double(ii)+1.))*d_temp;
                                dWigner[i-1]/=sin(vig_the);
                        }                                                                               // end for loop
                }                                                                                       // end singularity check

                // III.2 - else if (m!=0), no special case ----------------------------------------
                else{                                                                           // non-singular case : if (m!=0)
                        // A - if (i<n_min), then set VIG_d to zero -----------------------------------
        //              for(i=0; i<n_min; i++){                                 // for loop
                //              VIG_d[i]  =0.;                                          // == zero by default
                        //      d_VIG_d[i]=0.;                                          // == zero by default
        //              }                                                                               // end for loop
                        // B - obtain the readily availble VIG_d[n_min] value ------------------------
                        // vig_d ---------------------------------------------------------------------
                        Wigner[n_min] = vig_d_n_min;                            // by definition
                        // d_vig_d for the previous from current : d_VIG_d[i] = 0.*VIG_d[i-1] + 0.*VIG_d[i] + ...*VIG_d[i+1]    - eqn(B.26)
                        ii=n_min-1;
                        d_temp1 = sqrt((double(ii)+1.)*(double(ii)+1.) - double(l*l));
                        d_temp2 = sqrt((double(ii)+1.)*(double(ii)+1.) - double(m_*m_));
                        d_temp3 = (2.*double(ii)+1.);
                        d_temp4 = (double(ii)+1.);
                        dWigner[ii] = double(ii)*d_temp1*d_temp2/(d_temp3*d_temp4)*Wigner[n_min];
                        dWigner[ii]/=sin(vig_the);                              // == zero by definition
                        // C - obtain all other values in VIG_d[i] from recursive relationship --------
                        //   - Based on 'n_min' value and 'n' value; total recursive steps == n - n_min
                        for(i=n_min+1; i<=nMax+1; i++){                 // for loop
                                ii=i-1;
                                d_temp1 = sqrt((double(ii)+1.)*(double(ii)+1.) - double(l*l));
                                d_temp2 = sqrt((double(ii)+1.)*(double(ii)+1.) - double(m_*m_));
                                d_temp3 = (2.*double(ii)+1.)*( ii*(ii+1.)*vig_x - l*m_ );
                                d_temp4 = (double(ii)+1.);
                                d_temp5 = sqrt(double(ii*ii) - double(l*l));
                                d_temp6 = sqrt(double(ii*ii) - double(m_*m_));
                                // evaluate from terms above
                                d_temp  = (1./double(ii)/d_temp1/d_temp2)*( (d_temp3)*Wigner[i-1]
                                                - (d_temp4*d_temp5*d_temp6)*Wigner[i-2] );
                                // populate VIG_d[i] array
                                if(i<=nMax)
                                        Wigner[i] = d_temp;
                                // d_vig_d for the previous from current : d_VIG_d[i-1] = ...*VIG_d[i-2] + ...*VIG_d[i-1] + ...*VIG_d[i]        - eqn(B.26)
                                dWigner[ii] =   - d_temp4*d_temp5*d_temp6/(double(ii)*(2.*double(ii)+1.))*Wigner[i-2]
                                                                - 0.                                    // since l==0
                                                                + double(ii)*d_temp1*d_temp2/(d_temp4*(2.*double(ii)+1.))*d_temp;
                                dWigner[ii]/=sin(vig_the);
                        }                                                                                       // end : for loop
                }                                                                                               // end : non-singular case : if (m!=0)
                // -------------------------------------------------------------------------------

                // IV - if (m<0) : apply symmetry property eq(B.7) to eqs(B.22-B.24) -------------
                if(check_m_negative==1){                                                                                // solve using symmetry relation eq(B.7)
                        for(i=0; i<=nMax; i++){                                         // for loop
                                d_temp1    = 1./(pow(-1., double(i+l)));// this can be simplified, since l==0
                                Wigner[i]  *= d_temp1;                                  // obtain final VIG_d
                                dWigner[i]*= -d_temp1;                                  // obtain final VIG_d
                        }                                                                                       // end : for loop
                }                                                                                               // end : if (m<0) condition

                // ------------------------------------------------------------------------------

                return 0;

}

// from NR
void gauleg(double x1, double x2, double x[], double w[], int n){
// Given the lower and upper limits of integration x1 and x2, and given n, this routine returns
// arrays x[1..n] and w[1..n] of length n, containing the abscissas and weights of the Gauss-
// Legendre n-point quadrature formula.

	int m,j,i;
	double z1,z,xm,xl,pp,p3,p2,p1; 				// High precision is a good idea for this routine.
	m=(n+1)/2; 									// The roots are symmetric in the interval, so
	xm=0.5*(x2+x1); 							// we only have to find half of them.
	xl=0.5*(x2-x1);
	for (i=1;i<=m;i++){ 						// Loop over the desired roots.
		z=cos(3.141592654*(i-0.25)/(n+0.5));
		// Starting with the above approximation to the ith root, we enter the main loop of
		// refinement by Newton’s method.
		do {
			p1=1.0;
			p2=0.0;
			for (j=1;j<=n;j++){ 				// Loop up the recurrence relation to get the
				p3=p2; 							// Legendre polynomial evaluated at z.
				p2=p1;
				p1=((2.0*j-1.0)*z*p2-(double(j)-1.0)*p3)/j;
			}
			// p1 is now the desired Legendre polynomial. We next compute pp, its derivative,
			// by a standard relation involving also p2, the polynomial of one lower order.
			pp=n*(z*p1-p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp; 						// Newton’s method.
		} while (fabs(z-z1) > EPS);
		x[i]=xm-xl*z; 							// Scale the root to the desired interval,
		x[n+1-i]=xm+xl*z; 						// and put in its symmetric counterpart.
		w[i]=2.0*xl/((1.0-z*z)*pp*pp); 			// Compute the weight
		w[n+1-i]=w[i]; 							// and its symmetric counterpart.
	}
}
// ------------------------------------------------------------------------------------------

// ellipse **********************************************************************************
// obtain r wrt theta : r(the) --------------------------------------------------------------
double ellipe_rthe( double a, double b, double costhe){

	double value(0.);
	value = a*b;
	value*= pow((b*b-a*a)*costhe + a*a, -1./2.);
	return value;
}
// ------------------------------------------------------------------------------------------

// obtain the derivative of ellipe_rthe wrt theta : dr_dthe_r(the)---------------------------
double ellipe_drthe_dthe( double a, double b, double costhe){

	double value(0.);
	value = 0.5*a*b;
	value*= pow((b*b-a*a)*costhe + a*a, -3./2.);
	value*= 2.*(b*b-a*a)*costhe*pow(1. - costhe*costhe, 1./2.);
	return value;
}
// ------------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------------
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

// ellipsoid ********************************************************************************
// obtain r wrt theta and phi : r(the, phi)--------------------------------------------------
double ellipsoid_rthephi( double a, double b, double c, double the, double phi){

	double value(0.);
	double sinthe = sin(the);
	double costhe = cos(the);
	double sinphi = sin(phi);
	double cosphi = cos(phi);
	// evaluate r(the, phi)
	value = ( (sinthe*sinthe * cosphi*cosphi) / (a*a)
			+ (sinthe*sinthe * sinphi*sinphi) / (b*b)
			+ (costhe*costhe)/(c*c)						);
	value = pow(value, -1./2.);

	return value;
}
// ------------------------------------------------------------------------------------------

// obtain derivative of ellipsoid_rthephi wrt theta : dr_dthe_r(the, phi) -------------------
double ellipsoid_drthephi_dthe( double a, double b, double c, double the, double phi){

	double value(0.);
	double sinthe = sin(the);
	double costhe = cos(the);
	double sinphi = sin(phi);
	double cosphi = cos(phi);
	// evaluate dr_dthe_r(the, phi)
	value = ( (sinthe*sinthe * cosphi*cosphi) / (a*a)
			+ (sinthe*sinthe * sinphi*sinphi) / (b*b)
			+ (costhe*costhe)/(c*c)						);
	value = 2.*pow(value, 3./2.);
	value =-( (2.*sinthe*costhe * cosphi*cosphi) / (a*a)
			+ (2.*sinthe*costhe * sinphi*sinphi) / (b*b)
			- (2.*sinthe*costhe)/(c*c)					)/value;
	return value;
}
// ------------------------------------------------------------------------------------------

// obtain derivative of ellipsoid_rthephi wrt phi : dr_dthe_r(the, phi) ---------------------
double ellipsoid_drthephi_dphi( double a, double b, double c, double the, double phi){

	double value(0.);
	double sinthe = sin(the);
	double costhe = cos(the);
	double sinphi = sin(phi);
	double cosphi = cos(phi);
	// evaluate dr_dphi_r(the, phi)
	value = ( (sinthe*sinthe * cosphi*cosphi) / (a*a)
			+ (sinthe*sinthe * sinphi*sinphi) / (b*b)
			+ (costhe*costhe)/(c*c)						);
	value = 2.*pow(value, 3./2.);
	value =-(-(2.*sinthe*sinthe * sinphi*cosphi) / (a*a)
			+ (2.*sinthe*sinthe * sinphi*cosphi) / (b*b))/value;
	return value;
}
// ------------------------------------------------------------------------------------------
// ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// sphere ***********************************************************************************
Qclass calcTsphere(int n, int m, int np, int mp, complex<double> wavek1, complex<double> wavek2, double a){

	// Mie-Lorentz : equations 5.219 - 5.220

	// declare class ---------------------------------------------------------------
	Qclass computeJs,  computeTs;
	computeJs.Q11 = complex<double>(0., 0.);
	computeJs.Q12 = complex<double>(0., 0.);
	computeJs.Q21 = complex<double>(0., 0.);
	computeJs.Q22 = complex<double>(0., 0.);

	computeTs.Q11 = complex<double>(0., 0.);
	computeTs.Q12 = complex<double>(0., 0.);
	computeTs.Q21 = complex<double>(0., 0.);
	computeTs.Q22 = complex<double>(0., 0.);

	// start calculations -----------------------------------------------------------
	if(m == mp && n==np){
		// Aux variables ------------------------------------------------------------
		complex<double> wk1a(0., 0.), wk2a(0., 0.), mk2k1(0., 0.), x(0., 0.), mx(0., 0.);
		wk1a 	= wavek1*a;
		wk2a 	= wavek2*a;
		mk2k1	= wavek2/wavek1;
		x 		= wk1a;
		mx		= mk2k1*wk1a;
		complex<double> exymx(0., 0.), 	dexymx(0., 0.);
		complex<double> exyx(0., 0.),  	dexyx(0., 0.);
		complex<double> scyx(0., 0.),  	dscyx(0., 0.);
		complex<double> an(0., 0.), 	bn(0., 0.);
		// --------------------------------------------------------------------------

		// Hankel-Bessel evaluations ------------------------------------------------
		Bessel H, B, Bmx;
		int Hreg, Breg;
		complex<double> Hvalue(0., 0.);
		complex<double> dHvalue(0.,0.);
		complex<double> Bvalue(0., 0.);
		complex<double> dBvalue(0.,0.);
		complex<double> Bmxvalue(0., 0.);
		complex<double> dBmxvalue(0.,0.);
		// obtain Hvalue, B, bmx
		// H
		Hreg = 1;
		H.init(x, Hreg, 0, n);						// Initialize Bessel object : start from n==0
		H.populate();								// calculate Bessel object
		Hvalue  = H.data[n];						// obtain required value
		dHvalue = H.ddata[n];						// obtain required value
		// B
		Breg = 0;
		B.init(x, Breg, 0, n);						// Initialize Bessel object : start from n==0
		B.populate();								// calculate Bessel object
		Bvalue  = B.data[n];						// obtain required value
		dBvalue = B.ddata[n];						// obtain required value
		// B
		Breg = 0;
		Bmx.init(mx, Breg, 0, n);					// Initialize Bessel object : start from n==0
		Bmx.populate();								// calculate Bessel object
		Bmxvalue  = Bmx.data[n];					// obtain required value
		dBmxvalue = Bmx.ddata[n];					// obtain required value
		// ----------------------------------------------------------------------------

		// obtain exy, dexy, scy, dscy ------------------------------------------------
		scyx  	= x*Hvalue;
		dscyx 	= x*dHvalue+Hvalue;

		exyx  	= x*Bvalue;
		dexyx 	= x*dBvalue+Bvalue;

		exymx  	= mx*Bmxvalue;
		dexymx 	= mx*dBmxvalue+Bmxvalue;
		// ----------------------------------------------------------------------------

//		// obtain Js ------------------------------------------------------------------
//		computeJs.Q11 = complex<double>(0., 0.);
//		computeJs.Q12 = 1./wavek1/wavek1/mk2k1*dscyx*exymx;
//		computeJs.Q21 = -1./wavek1/wavek1/mk2k1*scyx*dexymx;
//		computeJs.Q22 = complex<double>(0., 0.);
//		// ----------------------------------------------------------------------------

		// obtain an, bn --------------------------------------------------------------
		// an
		an = 	(mk2k1*exymx*dexyx - exyx*dexymx) /
				(mk2k1*exymx*dscyx - scyx*dexymx);
		// bn
		bn = 	(mk2k1*exyx*dexymx - exymx*dexyx) /
				(mk2k1*scyx*dexymx - exymx*dscyx);
		// -----------------------------------------------------------------------------

		// obtain Ts -------------------------------------------------------------------
		computeTs.Q11 = complex<double>(-1., 0.)*bn;
		computeTs.Q12 = complex<double>(0., 0.);
		computeTs.Q21 = complex<double>(0., 0.);
		computeTs.Q22 = complex<double>(-1., 0.)*an;
		// -----------------------------------------------------------------------------
	}

//	return computeJs;
	return computeTs;

}

Qclass calcQ_rthe(int n, int m, int np, int mp, int N, int BHreg, complex<double> wavek1, complex<double> wavek2, double a, double b){

// np==n_prime == l in AlBe_nmlk
// mp==m_prime == k in AlBe_nmlk
// rotationally symmetric particle : equations 5.196 - 5.199 (Mishchenko book)
// rthe for spheroid

	// a is the projection on the z-axis by definition - circular cross section along z
	// b is the projection on the y-axis by definition - y==x
	// This routine is based on
	// a == semi-major axis
	// b == semi-minor axis
	// check if a != semi-major axis. If yes, switch a and b values
	if(b>a){
		double tempb = b;
		b = a;						// semi-minor axis
		a = tempb;					// semi-major axis
	}

	// declare class ---------------------------------------------------------------
	Qclass computeJ, computeQ;
	computeJ.Q11 = complex<double>(0., 0.);
	computeJ.Q12 = complex<double>(0., 0.);
	computeJ.Q21 = complex<double>(0., 0.);
	computeJ.Q22 = complex<double>(0., 0.);

	computeQ.Q11 = complex<double>(0., 0.);
	computeQ.Q12 = complex<double>(0., 0.);
	computeQ.Q21 = complex<double>(0., 0.);
	computeQ.Q22 = complex<double>(0., 0.);

	// start calculations -----------------------------------------------------------
	if(m == mp){
		// Aux variables ------------------------------------------------------------
		int i(0);
		double tempd(0.);
		double d_m(0.), d_n(0.), d_np(0.);
		complex<double> J11sum(0.), J12sum(0.);
		complex<double> J21sum(0.), J22sum(0.);
		complex<double> wk1rthe(0., 0.), wk2rthe(0., 0.);
		Spherical<double> R;
		double RR(0.);
		double costhe(0.), rthe0(0.), drdthe0(0.);
		// --------------------------------------------------------------------------

		// Hankel-Bessel evaluations ------------------------------------------------
		Bessel H, B;
		int Hreg, Breg;
		complex<double> Hvalue(0., 0.);
		complex<double> Bvalue(0., 0.);
		complex<double> dHvalue(0., 0.);
		complex<double> dBvalue(0., 0.);
		Breg = 0;
		if(BHreg){
			Hreg = 0;
		}
		else if(!BHreg){
			Hreg = 1;
		}
		// ----------------------------------------------------------------------------

		// Wigner function declarations -----------------------------------------------
		double *Wignern, *dWignern;
		Wignern	 = new double[n+1];
		dWignern = new double[n+1];
		double *Wignernp,*dWignernp;
		Wignernp = new double[np+1];
		dWignernp= new double[np+1];
		double PImn(0.), TAmn(0.), PImnp(0.), TAmnp(0.);
		// ----------------------------------------------------------------------------

		// evaluate J11 ---------------------------------------------------------------
		d_m = double(m);
		d_n = double(n);
		d_np= double(np);

		tempd = 2.*d_n+1.;
		tempd*= 2.*d_np+1.;
		tempd/= d_n*(d_n+1.);
		tempd/= d_np*(d_np+1.);
		tempd = sqrt(tempd);

		// parameters for the GL routine------------------------------------------------
	//	int N=10;									// number of samples in GL
		double *x, *w;
		x = new double [N+1];
		w = new double [N+1];
		// obtain Gauss-Legendre weights & abscissas -----------------------------------
		gauleg(-1, 1, x, w, N);

		// compute integral for Jsum(s)
		for(i=1; i<=N; i++){

			// get the value of (the) : abscissas --------------------------------------
			costhe = x[i];

			// obtain rthe, drdthe -----------------------------------------------------
			rthe0 	= ellipe_rthe(a, b, costhe);
			wk1rthe = wavek1*rthe0;
			wk2rthe = wavek2*rthe0;
			RR 		= rthe0*rthe0;
			drdthe0 = ellipe_drthe_dthe(a, b, costhe);

			// obtain Hvalue and Bvalue ------------------------------------------------
			// Hvalue
			H.init(wk1rthe, Hreg, 0, n);			// Initialize Bessel object : start from n==0
			H.populate();								// calculate Bessel object
			Hvalue = H.data[n];							// obtain required value
			dHvalue= H.ddata[n];						// obtain required value
			// Bvalue
			B.init(wk2rthe, Breg, 0, np);			// Initialize Bessel object : start from n==0
			B.populate();								// calculate Bessel object
			Bvalue = B.data[np];						// obtain required value
			dBvalue= B.ddata[np];						// obtain required value

			// obtain PImn, TAmn, PImnp, TAmnp -----------------------------------------
			// wrt n
			R.the = acos(x[i]);
			VIGdVIG(n, m, R, Wignern, dWignern);
			PImn  = d_m/sin(R.the)*Wignern[n];
			TAmn  = dWignern[n];
			// wrt np
			VIGdVIG(np, m, R, Wignernp, dWignernp);
			PImnp = d_m/sin(R.the)*Wignernp[np];
			TAmnp = dWignernp[n];

			// obtain J11, J12, J21, J22 ----------------------------------------------
			// evaluate J11sum
			J11sum += w[i]*RR*Hvalue*Bvalue*(PImn*TAmnp + TAmn*PImnp);

			// evaluate J12sum
			J12sum += w[i]*RR*Bvalue*(	(1./wk1rthe)*(wk1rthe*dHvalue+Hvalue)*(PImn*PImnp+TAmn*TAmnp)
										+ drdthe0*d_n*(d_n+1.)*Hvalue/wk1rthe*Wignern[n]*TAmnp	);

			// evaluate J12sum
			J21sum += w[i]*RR*Hvalue*(	(1./wk2rthe)*(wk2rthe*dBvalue+Bvalue)*(PImn*PImnp+TAmn*TAmnp)
										+ drdthe0*d_np*(d_np+1.)*Bvalue/wk2rthe*TAmn*Wignernp[n]);

			// evaluate J12sum
			J22sum += w[i]*RR*(			(1./wk1rthe)*(wk1rthe*dHvalue+Hvalue)*(1./wk2rthe)*(wk2rthe*dBvalue+Bvalue)*(PImn*TAmnp+TAmn*PImnp)
										+ drdthe0*(d_n*(d_n+1.)*Hvalue/wk1rthe*(1./wk2rthe)*(wk2rthe*dBvalue+Bvalue)
										+ (1./wk1rthe)*(wk1rthe*dHvalue+Hvalue)*d_np*(d_np+1.)*Bvalue/wk2rthe)*PImn*Wignernp[n]	);
		}

		// delete arrays
		delete [] Wignern;
		delete [] dWignern;
		delete [] Wignernp;
		delete [] dWignernp;

		// once the integral from GL is evaluated, obtain J
		// J11
		computeJ.Q11 = complex<double>( 0., -0.5*tempd)*J11sum;

		// J12
		computeJ.Q12 = complex<double>( 0.5*tempd,  0.)*J12sum;

		// J21
		computeJ.Q21 = complex<double>(-0.5*tempd,  0.)*J21sum;

		// J22
		computeJ.Q22 = complex<double>( 0., -0.5*tempd)*J22sum;
		// -----------------------------------------------------------------------------

		// obtain Q from J -------------------------------------------------------------
		computeQ.Q11 = complex<double>(0., -1.)*wavek1*wavek2*computeJ.Q21 + complex<double>(0., -1.)*wavek1*wavek1*computeJ.Q12;
		computeQ.Q12 = complex<double>(0., -1.)*wavek1*wavek2*computeJ.Q11 + complex<double>(0., -1.)*wavek1*wavek1*computeJ.Q22;
		computeQ.Q21 = complex<double>(0., -1.)*wavek1*wavek2*computeJ.Q22 + complex<double>(0., -1.)*wavek1*wavek1*computeJ.Q11;
		computeQ.Q22 = complex<double>(0., -1.)*wavek1*wavek2*computeJ.Q12 + complex<double>(0., -1.)*wavek1*wavek1*computeJ.Q21;
		// -----------------------------------------------------------------------------
	}

//	return computeJ;
	return computeQ;

}
/*
Qclass calcQ_rthephi(int n, int m, int np, int mp, int N, int BHreg, complex<double> wavek1, complex<double> wavek2, double a, double b, double c){

// arbitrarily shaped particle : equation 5.184
// rthephi for the general shape of a spheroid

	// a is the projection on the x-axis by definition
	// b is the projection on the y-axis by definition
	// c is the projection on the z-axis by definition

	// declare class ---------------------------------------------------------------
	Qclass computeJ, computeQ;
	computeJ.Q11 = complex<double>(0., 0.);
	computeJ.Q12 = complex<double>(0., 0.);
	computeJ.Q21 = complex<double>(0., 0.);
	computeJ.Q22 = complex<double>(0., 0.);

	computeQ.Q11 = complex<double>(0., 0.);
	computeQ.Q12 = complex<double>(0., 0.);
	computeQ.Q21 = complex<double>(0., 0.);
	computeQ.Q22 = complex<double>(0., 0.);

	// start calculations -----------------------------------------------------------

	// Aux variables ------------------------------------------------------------
	int i(0), j(0);
	double tempd(0.);
	complex<double> tempc11(0., 0.), tempc12(0., 0.), tempc21(0., 0.), tempc22(0., 0.);
	double d_m=double(m);
	complex<double> J11sum(0., 0.), J12sum(0., 0.);
	complex<double> J21sum(0., 0.), J22sum(0., 0.);
	complex<double> wk1rthe(0., 0.), wk2rthe(0., 0.);
	Spherical<double> R;
	double RR(0.);
	double the(0.), phi(0.), rthephi(0.), drdthe(0.), drdphi(0.);
	// --------------------------------------------------------------------------

	// Hankel-Bessel evaluations ------------------------------------------------
	int Breg(0), Hreg(0);
	// This is to be passed again - therefore keep to the logic convention!!! - this is a cause of confusion!
	Breg = 1;
	if(BHreg){
		Hreg = 1;
	}
	else if(!BHreg){
		Hreg = 0;
	}
	// ----------------------------------------------------------------------------

	// Spherical waves M & N declarations -----------------------------------------
	Spherical<complex<double> > Mnvalue, Nnvalue;
	Spherical<complex<double> > *Mn;			// M function arrays
	Spherical<complex<double> > *Nn;			// N function arrays
	Mn = new Spherical<complex<double> > [n+1];
	Nn = new Spherical<complex<double> > [n+1];
	Spherical<complex<double> > Mnpvalue, Nnpvalue;
	Spherical<complex<double> > *Mnp;			// M function arrays
	Spherical<complex<double> > *Nnp;			// N function arrays
	Mnp= new Spherical<complex<double> > [np+1];
	Nnp= new Spherical<complex<double> > [np+1];
	// ----------------------------------------------------------------------------

	// obtain Gauss-Legendre weights & abscissas -----------------------------------
	// the	: 0 <= the <= PI
	double bma2_the = PI/2.;
	double *xthe, *wthe;
	xthe = new double [N+1];
	wthe = new double [N+1];
	gauleg(0, PI, xthe, wthe, N);

	// phi	: 0 <= phi <= 2PI
	double bma2_phi = PI;
	double *xphi, *wphi;
	xphi = new double [2*N+1];
	wphi = new double [2*N+1];
	gauleg(0, 2*PI, xphi, wphi, 2*N);

	// compute integral for Jsum(s)
	for(i=1; i<=N; i++){							// the loop
		for(j=1; j<=2*N; j++){						// phi loop

			// get value of (the) and (phi) : abscissas --------------------------------
			// the = bma2_the*xthe[i] + bma2_the;		// GL normalization
			// phi = bma2_phi*xphi[j] + bma2_phi;		// GL normalization

			the = xthe[i];					// GL normalization
			phi = xphi[j];					// GL normalization

			// obtain rthephi, drdthe, drdphi ------------------------------------------
			rthephi = ellipsoid_rthephi(a, b, c, the, phi);
			wk1rthe = wavek1*rthephi;
			wk2rthe = wavek2*rthephi;
			RR 		= rthephi*rthephi;
			drdthe  = ellipsoid_drthephi_dthe(a, b, c, the, phi);
			drdphi  = ellipsoid_drthephi_dphi(a, b, c, the, phi);

			// obtain Mn, Nn, Mp, Np arrays --------------------------------------------
			// update R value in the and phi
			R.rrr = rthephi;
			R.the = the;
			R.phi = phi;
			// evaluate GL samples on the and phi
			// Mn, Nn 	- dependent on 'Hreg'
			compute_MnSurf(n, -m, Hreg, R, wavek1, Mn);
			compute_NnSurf(n, -m, Hreg, R, wavek1, Nn);
			Mnvalue = Mn[n];							// obtain required value
			Nnvalue = Nn[n];							// obtain required value

			// Mnp, Nnp	- dependent on 'Breg' - always regular
			compute_MnSurf(np, mp, Breg, R, wavek2, Mnp);
			compute_NnSurf(np, mp, Breg, R, wavek2, Nnp);
			Mnpvalue = Mnp[n];							// obtain required value
			Nnpvalue = Nnp[n];							// obtain required value

			// obtain J11, J12, J21, J22 ----------------------------------------------
			// evaluate J11sum
			tempc11 = RR*sin(the)*(  Mnpvalue.the*Mnvalue.phi
									-Mnpvalue.phi*Mnvalue.the );
			tempc11*= wthe[i]*wphi[j];								// GL normalization
			J11sum += tempc11;

			// evaluate J12sum
			tempc12 = RR*sin(the)*( Mnpvalue.the*Nnvalue.phi - Mnpvalue.phi*Nnvalue.the )
					 -R.rrr*drdthe*(Mnpvalue.phi*Nnvalue.rrr)*sin(the)
					 +R.rrr*drdphi*(Mnpvalue.the*Nnvalue.rrr);
			tempc12*= wthe[i]*wphi[j];								// GL-normalization
			J12sum += tempc12;

			// evaluate J21sum
			tempc21 = RR*sin(the)*( Nnpvalue.the*Mnvalue.phi - Nnpvalue.phi*Mnvalue.the )
					 +R.rrr*drdthe*(Nnpvalue.rrr*Mnvalue.phi)*sin(the)
					 -R.rrr*drdphi*(Nnpvalue.rrr*Mnvalue.the);
			tempc21*= wthe[i]*wphi[j];								// GL-normalization
			J21sum += tempc21;

			// evaluate J22sum
			tempc22 = RR*sin(the)*( Nnpvalue.the*Nnvalue.phi - Nnpvalue.phi*Nnvalue.the)
					 -R.rrr*drdthe*(Nnpvalue.phi*Nnvalue.rrr - Nnpvalue.rrr*Nnvalue.phi)*sin(the)
					 -R.rrr*drdphi*(Nnpvalue.rrr*Nnvalue.the - Nnpvalue.the*Mnvalue.rrr);
			tempc22*= wthe[i]*wphi[j];								// GL-normalization
			J22sum += tempc22;

			// evaluate surface area of particle - check GL
			tempd+=RR*sin(the)*wthe[i]*wphi[j];
		}	// j loop (phi)
	}	// i loops (the)

	// once the integral from GL is evaluated, obtain final J ----------------------
	// J11
	computeJ.Q11 = pow(-1, d_m)*J11sum;
	// J12
	computeJ.Q12 = pow(-1, d_m)*J12sum;
	// J21
	computeJ.Q21 = pow(-1, d_m)*J21sum;
	// J22
	computeJ.Q22 = pow(-1, d_m)*J22sum;
	// -----------------------------------------------------------------------------

	// obtain Q from J -------------------------------------------------------------
	computeQ.Q11 = complex<double>(0., -1.)*wavek1*wavek2*computeJ.Q21 + complex<double>(0., -1.)*wavek1*wavek1*computeJ.Q12;
	computeQ.Q12 = complex<double>(0., -1.)*wavek1*wavek2*computeJ.Q11 + complex<double>(0., -1.)*wavek1*wavek1*computeJ.Q22;
	computeQ.Q21 = complex<double>(0., -1.)*wavek1*wavek2*computeJ.Q22 + complex<double>(0., -1.)*wavek1*wavek1*computeJ.Q11;
	computeQ.Q22 = complex<double>(0., -1.)*wavek1*wavek2*computeJ.Q12 + complex<double>(0., -1.)*wavek1*wavek1*computeJ.Q21;
	// -----------------------------------------------------------------------------

	return computeJ;
//	return computeQ;

}
*/

