
#include "AJ_AuxFuns.h"

#include "AuxCoefficients.h"
#include "Bessel.h"
#include "CompoundIterator.h"
#include "Spherical.h"
#include "constants.h"

#include <cmath>
#include <iostream>
#include <cstdlib>


// -----------------------------------------------------------------------
// computes all corresponding Spherical Harmonics values given (nMax, m)
int compute_YJn_m(Spherical<double> R, std::complex<double> waveK, int BHreg, int nMax, int m, std::complex<double> *YJnm){


  int i(0), n(0);
  double d_n(0.);
  double d_temp(0.);

  double *dn;
  dn = new double[nMax+1];

  double *Wigner, *dWigner;
  Wigner  = new double[nMax+1];
  dWigner = new double[nMax+1];

  AuxCoefficients auxCoefficients;
  auxCoefficients.compute_dn(nMax, dn);
  auxCoefficients.VIGdVIG(nMax, m, R, Wigner, dWigner);

  // Initialize and populate Bessel object
  Bessel besselH_R;
  besselH_R.init(R.rrr*waveK, BHreg, 0, nMax);
  besselH_R.populate();

  double dm = pow(-1., double(m));          // Legendre to Wigner function
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
int compute_Yn_m(Spherical<double> R, std::complex<double>, int nMax, int m, std::complex<double> *Ynm){


  int i(0), n(0);
  double d_n(0.);
  double d_temp(0.);

  double *dn;
  dn = new double[nMax+1];

  double *Wigner, *dWigner;
  Wigner  = new double[nMax+1];
  dWigner = new double[nMax+1];

  AuxCoefficients auxCoefficients;
  auxCoefficients.compute_dn(nMax, dn);
  auxCoefficients.VIGdVIG(nMax, m, R, Wigner, dWigner);

  double dm = pow(-1., double(m));          // Legendre to Wigner function
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
