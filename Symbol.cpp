#include "Symbol.h"

#include <cmath>
#include "constants.h"
#include "Bessel.h"
#include "CompoundIterator.h"
#include "gsl/gsl_sf_coupling.h"

using std::sqrt;
using std::pow;

double Symbol::Wigner3j(int j1, int j2, int j3, int m1, int m2, int m3)
{
  return gsl_sf_coupling_3j(2*j1, 2*j2, 2*j3, 2*m1, 2*m2, 2*m3);
}

double Symbol::Wigner6j(int j1, int j2, int j3, int j4, int j5, int j6)
{
  return gsl_sf_coupling_6j(2*j1, 2*j2, 2*j3, 2*j4, 2*j5, 2*j6);
}

double Symbol::Wigner9j(int j11, int j12, int j13, int j21, int j22, int j23,
    int j31, int j32, int j33)
{
  return gsl_sf_coupling_9j(2*j11, 2*j12, 2*j13, 2*j21, 2*j22, 2*j23, 2*j31, 2*j32, 2*j33);
}

double Symbol::CleGor(int j, int m, int j1, int m1, int j2, int m2)
{
  return pow(-1.0, m+j1-j2) * sqrt(2.0 * j + 1.0) * Wigner3j(j1, j2, j, m1, m2, -m);
}

// AJ --------------------------------------------------------------------------------
// C numbers
double Symbol::C_01m1(int J1, int M1, int J2, int M2, int J, int M)
{
  return sqrt(3./2./consPi)*CleGor(J, M, J1, M1, J2, M2)*(

      sqrt((J1+1.)*(J2)*(2.*J1-1.)*(2.*J2-1.))*
      Wigner9j(J1, J1-1, 1, J2, J2-1, 1, J, J, 1)*
      CleGor(J, 0, J1-1, 0, J2-1, 0)

      - sqrt((J1+1.)*(J2+1.)*(2.*J1-1.)*(2.*J2+3.))*
      Wigner9j(J1, J1-1, 1, J2, J2+1, 1, J, J, 1)*
      CleGor(J, 0, J1-1, 0, J2+1, 0)

      + sqrt((J1)*(J2)*(2.*J1+3.)*(2.*J2-1.))*
      Wigner9j(J1, J1+1, 1, J2, J2-1, 1, J, J, 1)*
      CleGor(J, 0, J1+1, 0, J2-1, 0)

      - sqrt((J1)*(J2+1)*(2.*J1+3.)*(2.*J2+1.))*
      Wigner9j(J1, J1+1, 1, J2, J2+1, 1, J, J, 1)*
      CleGor(J, 0, J1+1, 0, J2+1, 0)              );
}

double Symbol::C_10m1(int J1, int M1, int J2, int M2, int J, int M)
{
  return sqrt(3./2./consPi)*(2.*J1+1.)*CleGor(J, M, J1, M1, J2, M2)*(

      sqrt((J2)*(2.*J1-1.))*
      Wigner9j(J1, J1, 1, J2, J2-1, 1, J, J+1, 1)*
      CleGor(J+1, 0, J1, 0, J2-1, 0)*
      sqrt((J)/(2.*J+1.))

      - sqrt((J2+1.)*(2.*J2+3.))*
      Wigner9j(J1, J1, 1, J2, J2+1, 1, J, J+1, 1)*
      CleGor(J+1, 0, J1, 0, J2+1, 0)*
      sqrt((J)/(2.*J+1.))

      + sqrt((J2)*(2.*J2-1.))*
      Wigner9j(J1, J1, 1, J2, J2-1, 1, J, J-1, 1)*
      CleGor(J-1, 0, J1, 0, J2-1, 0)*
      sqrt((J+1.)/(2.*J+1.))

      - sqrt((J2+1)*(2.*J1+3.))*
      Wigner9j(J1, J1, 1, J2, J2+1, 1, J, J-1, 1)*
      CleGor(J-1, 0, J1, 0, J2+1, 0)*
      sqrt((J+1)/(2.*J+1.))                   );
}

double Symbol::C_00m1(int J1, int M1, int J2, int M2, int J, int M)
{
  return sqrt(3./2./consPi)*(2.*J1+1.)*CleGor(J, M, J1, M1, J2, M2)*(

      sqrt((J2)*(2.*J1-1.))*
      Wigner9j(J1, J1, 1, J2, J2-1, 1, J, J, 1)*
      CleGor(J, 0, J1, 0, J2-1, 0)

      - sqrt((J2+1.)*(2.*J2+3.))*
      Wigner9j(J1, J1, 1, J2, J2+1, 1, J, J, 1)*
      CleGor(J+1, 0, J1, 0, J2+1, 0)              );
}

double Symbol::C_11m1(int J1, int M1, int J2, int M2, int J, int M)
{
  return sqrt(3./2./consPi)*CleGor(J, M, J1, M1, J2, M2)*(

      sqrt((J1+1.)*(J2)*(2.*J1-1.)*(2.*J2-1.))*
      Wigner9j(J1, J1-1, 1, J2, J2-1, 1, J, J+1, 1)*
      CleGor(J+1, 0, J1-1, 0, J2-1, 0)*
      sqrt((J)/(2.*J+1.))

      - sqrt((J1+1.)*(J2+1)*(2.*J1-1.)*(2.*J2+3.))*
      Wigner9j(J1, J1-1, 1, J2, J2+1, 1, J, J+1, 1)*
      CleGor(J+1, 0, J1-1, 0, J2+1, 0)*
      sqrt((J)/(2.*J+1.))

      + sqrt((J1)*(J2)*(2.*J1+3.)*(2.*J2-1.))*
      Wigner9j(J1, J1+1, 1, J2, J2-1, 1, J, J+1, 1)*
      CleGor(J+1, 0, J1+1, 0, J2-1, 0)*
      sqrt((J)/(2.*J+1.))

      - sqrt((J1)*(J2+1)*(2.*J1+3.)*(2.*J2+3.))*
      Wigner9j(J1, J1+1, 1, J2, J2+1, 1, J, J+1, 1)*
      CleGor(J+1, 0, J1+1, 0, J2+1, 0)*
      sqrt((J)/(2.*J+1.))

      + sqrt((J1+1.)*(J2)*(2.*J1-1.)*(2.*J2-1.))*
      Wigner9j(J1, J1-1, 1, J2, J2-1, 1, J, J-1, 1)*
      CleGor(J-1, 0, J1-1, 0, J2-1, 0)*
      sqrt((J+1)/(2.*J+1.))

      - sqrt((J1+1.)*(J2+1)*(2.*J1-1.)*(2.*J2+3.))*
      Wigner9j(J1, J1-1, 1, J2, J2+1, 1, J, J-1, 1)*
      CleGor(J-1, 0, J1-1, 0, J2+1, 0)*
      sqrt((J+1)/(2.*J+1.))

      + sqrt((J1)*(J2)*(2.*J1+3.)*(2.*J2-1.))*
      Wigner9j(J1, J1+1, 1, J2, J2-1, 1, J, J-1, 1)*
      CleGor(J-1, 0, J1+1, 0, J2-1, 0)*
      sqrt((J+1)/(2.*J+1.))

      - sqrt((J1)*(J2+1)*(2.*J1+3.)*(2.*J2+3.))*
      Wigner9j(J1, J1+1, 1, J2, J2+1, 1, J, J-1, 1)*
      CleGor(J-1, 0, J1+1, 0, J2+1, 0)*
      sqrt((J+1)/(2.*J+1.))                   );
}

// A numbers
complex<double> Symbol::A_0(int, int n, double R, complex<double> waveK_i, complex<double> cmn_1)
{
  int BHreg =1;
  Bessel besselH;
  besselH.init(R*waveK_i, BHreg, 0, n);
  besselH.populate();

  return besselH.data[n]*cmn_1;

}

complex<double> Symbol::A_1(int, int n, double R, complex<double> waveK_i, complex<double> dmn_1)
{
  int BHreg =1;
  Bessel besselH;
  besselH.init(R*waveK_i, BHreg, 0, n);
  besselH.populate();
  // to be double checked - d/dr(J(kr)) --------------------
  return complex<double>(0., 1.)*(1./waveK_i)*
      (waveK_i*besselH.ddata[n]+besselH.data[n]/R)*dmn_1;

}

complex<double> Symbol::A_m1(int, int n, double R, complex<double> waveK_i, complex<double> dmn_1)
{
  int BHreg =1;
  Bessel besselH;
  besselH.init(R*waveK_i, BHreg, 0, n);
  besselH.populate();

  return complex<double>(0., 1.)*sqrt(n*(n+1.))*(1./waveK_i/R)*
      (besselH.data[n])*dmn_1;

}

double Symbol::W(int L1, int J1, int M1, int L2, int J2, int M2, int L, int M)
{
  return pow(-1., J2+L1+L)*
      sqrt((2.*J1+1.)*(2.*J2+1.)*(2.*L1+1.)*(2.*L2+1.)/(4/consPi/(2.*L+1)))*
      Wigner6j(L1, L2, L, J2, J1, 1)*
      CleGor(L, 0, L1, 0, L2, 0)*
      CleGor(L, M, J1, M1, J2, M2);

}

// In what follows, need to devise CompoundIterator to start from n=0.
// otherwise, use double for loops

complex<double> Symbol::up_mn(int m, int n, int nMax_, complex<double> cmn_1, complex<double> dmn_1, double omega_, Scatterer *object_, ElectroMagnetic bground_)
{
  // Basic relations ------------------------------------
  double R        = object_->radius;
  complex<double> mu_0    = consMu0;
  complex<double> mu_b    = bground_.mu;
  complex<double> eps_b   = bground_.epsilon;
  // SH Basic relations ---------------------------------
  complex<double> b_non = object_->elmag.b_SH;

  // Auxiliary variables --------------------------------------------------------------------------------
  complex<double> waveK_j1 = (omega_ / 2.0) * sqrt(object_->elmag.epsilon * object_->elmag.mu);

  int m1(0), n1(0), m2(0), n2(0);
  CompoundIterator p;
  CompoundIterator q;
  complex<double> sum(0., 0.);

  for (p = CompoundIterator(nMax_, nMax_); p < p.max(nMax_); p++) {
    n1 = p.first;
    m1 = p.second;
    for (q = CompoundIterator(nMax_, nMax_); q < q.max(nMax_); q++) {
      n2 = q.first;
      m2 = q.second;
      sum +=  A_1(m1, n1, R, waveK_j1, dmn_1)*
          A_m1(m2, n2, R, waveK_j1, dmn_1)*
          C_10m1(n1, m1, n2, m2, n, m)

        + A_0(m1, n1, R, waveK_j1, cmn_1)*
          A_m1(m2, n2, R, waveK_j1, dmn_1)*
          C_11m1(n1, m1, n2, m2, n, m);
    }
  }
  return sum*complex<double>(0., 1.)*(-b_non/2.)*sqrt(mu_b/eps_b)/sqrt(mu_0/eps_b);
}

complex<double> Symbol::vp_mn(int m, int n, int nMax_, complex<double> cmn_1, complex<double> dmn_1, double omega_, Scatterer *object_, ElectroMagnetic bground_)
{
  // Basic relations ------------------------------------
  double R        = object_->radius;
  complex<double> mu_0    = consMu0;
  complex<double> eps_0   = consEpsilon0;
  complex<double> mu_b    = bground_.mu;
  complex<double> eps_b   = bground_.epsilon;
  // SH Basic relations ---------------------------------
  complex<double> b_non = object_->elmag.b_SH;

  // Auxiliary variables --------------------------------------------------------------------------------
  complex<double> waveK_j1  = (omega_ / 2.0) * sqrt(object_->elmag.epsilon * object_->elmag.mu);

  int m1(0), n1(0), m2(0), n2(0);
  CompoundIterator p;
  CompoundIterator q;
  complex<double> sum(0., 0.);
  for (p = CompoundIterator(nMax_, nMax_); p < p.max(nMax_); p++) {
    n1 = p.first;
    m1 = p.second;
    for (q = CompoundIterator(nMax_, nMax_); q < q.max(nMax_); q++) {
      n2 = q.first;
      m2 = q.second;
      sum +=  A_1(m1, n1, R, waveK_j1, dmn_1)*
          A_m1(m2, n2, R, waveK_j1, dmn_1)*
          C_00m1(n1, m1, n2, m2, n, m)

        + A_0(m1, n1, R, waveK_j1, cmn_1)*
          A_m1(m2, n2, R, waveK_j1, dmn_1)*
          C_01m1(n1, m1, n2, m2, n, m);
    }
  }
  return sum*complex<double>(-2., 0.)*(-b_non/2.)*sqrt(mu_b/eps_b)/sqrt(mu_0/eps_0);
}

complex<double> Symbol::upp_mn(int m, int n, int nMax_, complex<double> cmn_1, complex<double> dmn_1, double omega_, Scatterer *object_, ElectroMagnetic)
{

  // Basic relations ------------------------------------
  double R        = object_->radius;
  complex<double> eps_0   = consEpsilon0;
  complex<double> mu_j    = object_->elmag.mu;
  complex<double> eps_j   = object_->elmag.epsilon;
  // SH Basic relations ---------------------------------
  complex<double> a_non = object_->elmag.a_SH;
  complex<double> d_non = object_->elmag.d_SH;
  complex<double> eps_j2  = object_->elmag.epsilon_SH;

  // Auxiliary variables --------------------------------------------------------------------------------
  complex<double> waveK_01  = (omega_ / 2.0) * sqrt(eps_j * mu_j);
  complex<double> waveK_j1  = (omega_ / 2.0) * sqrt(eps_j * mu_j);


  int m1(0), n1(0), m2(0), n2(0);
  CompoundIterator p;
  CompoundIterator q;
  complex<double> gmn(0., 0.);
  complex<double> fmn(0., 0.);

  for (p = CompoundIterator(nMax_, nMax_); p < p.max(nMax_); p++) {
    n1 = p.first;
    m1 = p.second;
    for (q = CompoundIterator(nMax_, nMax_); q < q.max(nMax_); q++) {
      n2 = q.first;
      m2 = q.second;

      gmn +=  A_m1(m1, n1, R, waveK_j1, dmn_1)*
          A_m1(m2, n2, R, waveK_j1, dmn_1)*
          ( sqrt((n1)  /(2.*n1+1.))*sqrt((n2)  /(2.*n2+1.))*W(n1-1, n1, m1, n2-1, n2, m2, n, m)
          + sqrt((n1+1)/(2.*n1+1.))*sqrt((n2+1)/(2.*n2+1.))*W(n1+1, n1, m1, n2+1, n2, m2, n, m)
          - sqrt((n1)  /(2.*n1+1.))*sqrt((n2+1)/(2.*n2+1.))*W(n1-1, n1, m1, n2+1, n2, m2, n, m)
          - sqrt((n1+1)/(2.*n1+1.))*sqrt((n2)  /(2.*n2+1.))*W(n1+1, n1, m1, n2-1, n2, m2, n, m));

      fmn +=  A_1(m1, n1, R, waveK_j1, dmn_1)*
          A_1(m2, n2, R, waveK_j1, dmn_1)*
          ( sqrt((n1+1)/(2.*n1+1.))*sqrt((n2+1)/(2.*n2+1.))*W(n1-1, n1, m1, n2-1, n2, m2, n, m)
          + sqrt((n1)  /(2.*n1+1.))*sqrt((n2)  /(2.*n2+1.))*W(n1+1, n1, m1, n2+1, n2, m2, n, m)
          - sqrt((n1+1)/(2.*n1+1.))*sqrt((n2)  /(2.*n2+1.))*W(n1-1, n1, m1, n2+1, n2, m2, n, m)
          - sqrt((n1)  /(2.*n1+1.))*sqrt((n2+1)/(2.*n2+1.))*W(n1+1, n1, m1, n2-1, n2, m2, n, m)
        + A_0(m1, n1, R, waveK_j1, cmn_1)*
          A_0(m2, n2, R, waveK_j1, cmn_1)*
          W(n1, n1, m1, n2, n2, m2, n, m)
        + A_1(m1, n1, R, waveK_j1, dmn_1)*
          A_0(m2, n2, R, waveK_j1, cmn_1)*
          ( sqrt((n1+1)/(2.*n1+1.))*W(n1-1, n1, m1, n2, n2, m2, n, m)
          + sqrt((n1)  /(2.*n1+1.))*W(n1+1, n1, m1, n2, n2, m2, n, m))
        + A_0(m1, n1, R, waveK_j1, cmn_1)*
          A_1(m2, n2, R, waveK_j1, dmn_1)*
          ( sqrt((n1+1)/(2.*n1+1.))*W(n1, n1, m1, n2-1, n2, m2, n, m)
          + sqrt((n2)  /(2.*n1+1.))*W(n1, n1, m1, n2+1, n2, m2, n, m))            );
    }
  }
  return    complex<double>(0., 1.)*(-a_non/4.)*sqrt(n*(n+1))*gmn/waveK_01/R
      + complex<double>(0., 1.)*(-d_non/8.)*(eps_0/eps_j2)*sqrt(n*(n+1))*(gmn+fmn)/waveK_01/R;
}


// obtain the scattering and internal coefficients at (2w)


// ------------------------------------------------------------------------------------------------
