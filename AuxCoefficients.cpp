#include "AuxCoefficients.h"

#include "constants.h"
#include "Tools.h"
#include "Legendre.h"
#include "Bessel.h"
#include "CompoundIterator.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "gsl/gsl_sf_gamma.h"

int AuxCoefficients::compute_dn(int nMax, double *dn) {

  double d_n(0.);

  dn[0] = -1000; // this should return infinity!
  for (int i = 1; i <= nMax; i++) {
    d_n = double(i);
    dn[i] = std::sqrt((2. * d_n + 1.) / (4 * consPi * d_n * (d_n + 1.)));
  }

  return 0;
}

int AuxCoefficients::compute_Pn(int nMax, double *Wigner,
                                SphericalP<std::complex<double>> *Pn) {
  /*-------------------------------------------------------------------------------*/
  /* PURPOSE: Evaluate P_nm function from vig_d  and d_vig_d
   * --------------------- */
  /*-------------------------------------------------------------------------------*/
  for (int i = 0; i <= nMax; i++) {
    //    Pn[i].rrr = exp_imphi * Wigner[i];
    Pn[i].rrr = Wigner[i];
    Pn[i].the = 0.;
    Pn[i].phi = 0.;
  }

  return 0;
} // end P_func

int AuxCoefficients::compute_Pp(Spherical<double> R, int nMax,
                                SphericalP<std::complex<double>> *dataPp) {

  double *dn;
  dn = new double[nMax + 1];

  // Wigner d function test
  double *Wigner, *dWigner;
  Wigner = new double[nMax + 1];
  dWigner = new double[nMax + 1];

  // Vector spherical functions
  SphericalP<std::complex<double>> *Pn; // P function arrays
  Pn = new SphericalP<std::complex<double>>[nMax + 1];

  CompoundIterator p;
  CompoundIterator q;

  double d_n(0.0);
  double d_temp(0.0);

  for (q = CompoundIterator(nMax, nMax); q < q.max(nMax); q++) {

    // prepare for spherical functions calculation
    VIGdVIG(nMax, q.second, R, Wigner, dWigner);
    compute_dn(nMax, dn);

    // call vector spherical functions
    compute_Pn(nMax, Wigner, Pn);

    for (int n = abs(q.second); n <= nMax; n++) {
      if (n != 0) {
        d_n = double(n);

        double dm = pow(-1., double(q.second)); // Legendre to Wigner function
        std::complex<double> exp_imphi = exp(consCi * (double)q.second * R.phi);

        d_temp = dm * dn[n] * std::sqrt(d_n * (d_n + 1.));

        p.init(n, q.second);

        // dataPp[p] = Tools::toProjection(R, Pn[n] * (exp_imphi * d_temp));
        dataPp[p] = Pn[n] * (exp_imphi * d_temp);
      }
    }
  }

  delete[] Pn;
  delete[] Wigner;
  delete[] dWigner;
  delete[] dn;

  return 0;
}
// ----------------------------------------------------------------------------------
int AuxCoefficients::compute_Cn(int nMax, int m_, Spherical<double> R,
                                double *Wigner, double *dWigner,
                                SphericalP<std::complex<double>> *Cn) {
  /*-------------------------------------------------------------------------------*/
  /* PURPOSE: Evaluate P_nm function from vig_d  and d_vig_d
   * --------------------- */
  /*-------------------------------------------------------------------------------*/

  int i(0);
  double A(0.), B(0.);

  for (i = 0; i <= nMax; i++) {
    if (m_ == 0)
      A = 0.0;
    else if (std::abs(R.the) < 1e-10 ||
             (std::abs(R.the) - consPi + 1e-10) > 0.0)
      A = double(m_) / cos(R.the) * dWigner[i];
    else
      A = double(m_) / sin(R.the) * Wigner[i];

    B = dWigner[i];
    Cn[i].rrr = std::complex<double>(0., 0);
    Cn[i].the = std::complex<double>(0., A);
    Cn[i].phi = std::complex<double>(-B, 0.);
  }

  return 0;
} // end P_func
// ---------------------------------------------------------------------

// ----------------------------------------------------------------------------------
int AuxCoefficients::compute_Bn(int nMax, int m_, Spherical<double> R,
                                double *Wigner, double *dWigner,
                                SphericalP<std::complex<double>> *Bn) {
  /*-------------------------------------------------------------------------------*/
  /* PURPOSE: Evaluate P_nm function from vig_d  and d_vig_d
   * --------------------- */
  /*-------------------------------------------------------------------------------*/

  int i(0);
  double A(0.), B(0.);

  for (i = 0; i <= nMax; i++) {

    if (m_ == 0)
      A = 0.0;
    else if (std::abs(R.the) < 1e-10 ||
             (std::abs(R.the) - consPi + 1e-10) > 0.0)
      A = double(m_) / cos(R.the) * dWigner[i];
    else
      A = double(m_) / sin(R.the) * Wigner[i];

    B = dWigner[i];

    Bn[i].rrr = std::complex<double>(0., 0);
    Bn[i].the = std::complex<double>(B, 0.);
    Bn[i].phi = std::complex<double>(0., A);
  }

  return 0;
} // end P_func
// ---------------------------------------------------------------------

// ----------------------------------------------------------------------------------
int AuxCoefficients::compute_Mn(int nMax, int m_, Spherical<double> R,
                                std::complex<double> waveK, double *dn,
                                SphericalP<std::complex<double>> *Cn,
                                SphericalP<std::complex<double>> *Mn,
                                int BHreg) {
  /*-------------------------------------------------------------------------------*/
  /* PURPOSE: Evaluate M_nm function from c_sph C_nm[] and bh functions
   * ---------- */
  /* bh_indx == 1 -> use bessj
   * --------------------------------------------------- */
  /* bh_indx == 2 -> use bessy
   * --------------------------------------------------- */
  /* bh_indx == 3 -> use bessh1
   * -------------------------------------------------- */
  /* bh_indx == 4 -> use bessh2
   * -------------------------------------------------- */
  /* M_nm(kR) = (-1)^m * dn * bh(kR) * C_nm * exp(i*m*phi)  - eq(7) ref paper --
   */
  /*-------------------------------------------------------------------------------*/

  int i(0);
  std::complex<double> c_temp(0.0, 0.0);

  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) =
      optimet::bessel(R.rrr * waveK, (optimet::BESSEL_TYPE)BHreg, 0, nMax);

  double dm = pow(-1., double(m_)); // Legendre to Wigner function
  std::complex<double> exp_imphi(cos(double(m_) * R.phi),
                                 sin(double(m_) * R.phi));

  for (i = 0; i <= nMax; i++) {

    c_temp = dm * dn[i] * exp_imphi;

    Mn[i].rrr = c_temp * data[i] * Cn[i].rrr;
    Mn[i].the = c_temp * data[i] * Cn[i].the;
    Mn[i].phi = c_temp * data[i] * Cn[i].phi;
  }

  return 0;
} // end P_func
// ---------------------------------------------------------------------

// ----------------------------------------------------------------------------------
int AuxCoefficients::compute_Nn(int nMax, int m_, Spherical<double> R,
                                std::complex<double> waveK, double *dn,
                                SphericalP<std::complex<double>> *Pn,
                                SphericalP<std::complex<double>> *Bn,
                                SphericalP<std::complex<double>> *Nn,
                                int BHreg) {
  /*-------------------------------------------------------------------------------*/
  /* PURPOSE: Evaluate N_nm function from P_nm & B_nm[] and bh functions
   * --------- */
  /* bh_indx == 1 -> use bessj
   * --------------------------------------------------- */
  /* bh_indx == 2 -> use bessy
   * --------------------------------------------------- */
  /* bh_indx == 3 -> use bessh1
   * -------------------------------------------------- */
  /* bh_indx == 4 -> use bessh2
   * -------------------------------------------------- */
  /* N_nm(kR) = (-1)^m * dn * exp(i*m*phi) *     -- eq(8) ref paper --       */
  /*        [ n*(n+1)/(kR)bh(kR)*P_nm + 1./(kR)*(kR*derbh(kR)+bh(kR))*B_nm ] */
  /*-------------------------------------------------------------------------------*/

  int i(0);
  double d_n(0.);

  std::complex<double> Kr = waveK * R.rrr;
  std::complex<double> c_temp(0., 0.);

  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) =
      optimet::bessel(R.rrr * waveK, (optimet::BESSEL_TYPE)BHreg, 0, nMax);

  double dm = pow(-1., double(m_)); // Legendre to Wigner function
  std::complex<double> exp_imphi(cos(double(m_) * R.phi),
                                 sin(double(m_) * R.phi));

  Nn[0].rrr = std::complex<double>(0., 0.);
  Nn[0].the = std::complex<double>(0., 0.);
  Nn[0].phi = std::complex<double>(0., 0.);

  for (i = 1; i <= nMax; i++) {
    d_n = double(i);
    Nn[i].rrr = (1. / Kr) * dm * dn[i] *
                ((d_n * (d_n + 1.) * data[i] * Pn[i].rrr) +
                 ((Kr * ddata[i] + data[i]) * Bn[i].rrr)) *
                exp_imphi;
    Nn[i].the = (1. / Kr) * dm * dn[i] *
                ((d_n * (d_n + 1.) * data[i] * Pn[i].the) +
                 ((Kr * ddata[i] + data[i]) * Bn[i].the)) *
                exp_imphi;
    Nn[i].phi = (1. / Kr) * dm * dn[i] *
                ((d_n * (d_n + 1.) * data[i] * Pn[i].phi) +
                 ((Kr * ddata[i] + data[i]) * Bn[i].phi)) *
                exp_imphi;
  }
  return 0;
} // end P_func
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------------------
int AuxCoefficients::compute_MpNp(Spherical<double> R,
                                  std::complex<double> waveK, int BHreg,
                                  int nMax,
                                  SphericalP<std::complex<double>> *dataMp,
                                  SphericalP<std::complex<double>> *dataNp) {

  //  double d_temp=0.;

  // Wigner d function test
  // --------------------------------------------------------------
  double *Wigner, *dWigner;
  Wigner = new double[nMax + 1];
  dWigner = new double[nMax + 1];

  // Vector spherical functions
  // ----------------------------------------------------------
  SphericalP<std::complex<double>> *Pn; // P function arrays
  SphericalP<std::complex<double>> *Cn; // C function arrays
  SphericalP<std::complex<double>> *Bn; // B function arrays
  Pn = new SphericalP<std::complex<double>>[nMax + 1];
  Cn = new SphericalP<std::complex<double>>[nMax + 1];
  Bn = new SphericalP<std::complex<double>>[nMax + 1];

  // Vector spherical Waves
  // --------------------------------------------------------------
  SphericalP<std::complex<double>> *Mn; // M function arrays
  SphericalP<std::complex<double>> *Nn; // N function arrays
  Mn = new SphericalP<std::complex<double>>[nMax + 1];
  Nn = new SphericalP<std::complex<double>>[nMax + 1];
  // -------------------------------------------------------------------------------------

  CompoundIterator p;
  CompoundIterator q;

  for (q = CompoundIterator(nMax, nMax); q < q.max(nMax); q++) {

    VIGdVIG(nMax, q.second, R, Wigner, dWigner);
    compute_dn(nMax, dn);

    // call vector spherical functions
    // -----------------------------------------------------
    compute_Pn(nMax, Wigner, Pn);
    compute_Cn(nMax, q.second, R, Wigner, dWigner, Cn);
    compute_Bn(nMax, q.second, R, Wigner, dWigner, Bn);

    // call vector spherical waves
    // ---------------------------------------------------------
    compute_Mn(nMax, q.second, R, waveK, dn, Cn, Mn, BHreg);
    compute_Nn(nMax, q.second, R, waveK, dn, Pn, Bn, Nn, BHreg);

    for (int n = std::abs(q.second); n <= nMax; n++) {
      if (n != 0) {

        p.init(n, q.second);

        dataBp[p] = Tools::toProjection(R, Bn[n]);
        dataCp[p] = Tools::toProjection(R, Cn[n]);
        dataMp[p] = Tools::toProjection(R, Mn[n]);
        dataNp[p] = Tools::toProjection(R, Nn[n]);
      }
    }
  }

  delete[] Pn;
  delete[] Cn;
  delete[] Bn;
  delete[] Mn;
  delete[] Nn;

  delete[] Wigner;
  delete[] dWigner;

  return 0;
}

int AuxCoefficients::VIGdVIG(int nMax, int m_, Spherical<double> R,
                             double *Wigner, double *dWigner) {
  /*-------------------------------------------------------------------------------*/
  /* PURPOSE: Evaluate Wigner d function & its derivative : vig_d  and d_vig_d
   * --- */
  /* VIG_d(n)             = vig_d  (l=0, m, n, the)
   * -------------------------------------- */
  /* d_VIG_d(n)   = d_vig_d(l=0, m, n, the)
   * -------------------------------------- */
  /*-------------------------------------------------------------------------------*/

  // variables
  // ---------------------------------------------------------------------
  // loop variables
  // ----------------------------------------------------------------
  int i(0), ii(0);
  int check_m_negative(0); // flag for m<0
  // temp variables
  // ---------------------------------------------------------------
  double d_temp(0.); // double temp
  double d_temp1(0.), d_temp2(0.), d_temp3(0.), d_temp4(0.), d_temp5(0.),
      d_temp6(0.);
  // Wigner function indices
  // -------------------------------------------------------
  int n_min(0);       // vector spherical function index              : n>=1,
                      // n_min=max(|l|,|m|)
  int l(0), l_abs(0); // wigner function index
  int m_abs(0);       // vector spherical function index
                      // Wigner function variables
                      // -----------------------------------------------------
  // double vig_d(0.); // wigner function value
  double vig_the = R.the; // wigner function argument : [0<=the<=PI]
  double vig_x(0.);       // wigner function auxiliary variable   : x=cos(the)
  double vig_exy_lm(0.);  // wigner function auxiliary variable   :
                          // exy_lm=eq(B.16)                       - book
                          // reference
  // Wigner recursive relationship terms - to be used for finally obtaining
  // vig_d -
  double vig_d_n_min(0.); // wigner function recursive term               :
                          // d_n_min=eq(B.24)                      - book
                          // reference
  // -------------------------------------------------------------------------------

  // ------------------ prepare for VIG functions evaluations
  // ----------------------
  // -------------------------------------------------------------------------------
  for (i = 0; i <= nMax; i++) { // rest
    Wigner[i] = 0.;
    dWigner[i] = 0.;
  }
  // prepare for vig_d evaluation
  // --------------------------------------------------
  l_abs = std::abs(l);
  m_abs = std::abs(m_);
  // -------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------
  // I - determine n_min : n_min=max(|l|,|m|)
  // --------------------------------------
  //      if(l_abs>n_min)
  //    n_min = l_abs;
  //      if(m_abs>n_min)
  n_min = m_abs; // this is always true, since l==0
                 // -------------------------------------------------------------------------------

  // -------------------------------------------------------------------------------
  // II - obtain vig_d_n_min ---------------------------------------------------
  // 0 - set : vig_the==vig_the
  //      if(m_>=0) // solve directly eqs(B.22-B24)
  //              vig_the = vig_the; // keep as it is
  if (m_ < 0) {                 // solve using symmetry relation eq(B.7)
    vig_the = consPi - vig_the; // modify value of vig_x   : eq(B.7)
    m_ = std::abs(m_);          // obtain solution for |m| : eq(B.7)
    check_m_negative = 1;       // check m<0 : to use in IV below
  }
  vig_x = cos(vig_the); // calculate in radians

  // A - evaluate vig_exy_lm                                      // wigner
  // function recursive term in d_n_min=eq(B.23)   : exy_lm=eq(B.16)       -
  // book reference
  if (m_ >= l)
    vig_exy_lm = 1.; // this is always true, since l==0 & m>=0
  else if (m_ < l)
    vig_exy_lm = pow(-1., double(l - m_));

  // B - evaluate all other terms in d_n_min = eqn(B.24) -----------------------
  d_temp1 = pow(2., -n_min);
  //
  d_temp2 = std::sqrt(double(gsl_sf_fact(2 * n_min)));
  d_temp3 = std::sqrt(double(gsl_sf_fact(std::abs(l - m_))));
  d_temp4 = std::sqrt(double(gsl_sf_fact(std::abs(l + m_))));
  //
  d_temp5 = pow(1. - vig_x, std::abs(double(l - m_)) / 2.);
  d_temp6 = pow(1. + vig_x, std::abs(double(l + m_)) / 2.);

  // C - evaluate vig_d_n_min --------------------------------------------------
  // vig_d_n_min = vig_exy_lm*d_temp1*(d_temp2/d_temp3/d_temp4)*d_temp5*d_temp6;
  vig_d_n_min = vig_exy_lm;
  vig_d_n_min *= d_temp1;
  vig_d_n_min *= d_temp2;
  vig_d_n_min /= d_temp3;
  vig_d_n_min /= d_temp4;
  vig_d_n_min *= d_temp5;
  vig_d_n_min *= d_temp6;
  // ----------------------------------------------------------------------------

  // ----------------------------------------------------------------------------
  // III - calculate VIG_d & d_VIG_d
  // --------------------------------------------
  // III.1 - check for singularity cases - i.e (m==0)
  // ---------------------------
  if (m_ == 0) { // singularity check
                 // A - if (i<n_min), then set VIG_d to zero
                 // -----------------------------------
    //              for(i=0; i<n_min; i++){                                 //
    //              for loop
    //              VIG_d[i]  =0.;                                          //
    //              == zero by default
    //} // end for loop
    // B - obtain the readily availble VIG_d[n_min] value
    // ------------------------
    // vig_d
    // ---------------------------------------------------------------------
    Wigner[n_min] = vig_d_n_min; // by definition
    // C - obtain all other values in VIG_d[i] from recursive relationship
    // --------
    //   - Based on 'n_min' value and 'n' value; total recursive steps == n -
    //   n_min
    // obtain first term in recursive relationship eq(B.22) - from special case
    Wigner[n_min + 1] = vig_x;
    // d_vig_d for the previous from current : d_VIG_d[i] = 0.*VIG_d[i-1] +
    // 0.*VIG_d[i] + ...*VIG_d[i+1]    - eqn(B.26)
    ii = n_min;
    d_temp1 = std::sqrt((double(ii) + 1.) * (double(ii) + 1.) - double(l * l));
    d_temp2 =
        std::sqrt((double(ii) + 1.) * (double(ii) + 1.) - double(m_ * m_));
    d_temp3 = (2. * double(ii) + 1.);
    d_temp4 = (double(ii) + 1.);
    dWigner[ii] =
        double(ii) * d_temp1 * d_temp2 / (d_temp3 * d_temp4) * Wigner[ii + 1];
    dWigner[ii] /= sin(vig_the); // == zero by definition
    // obtain all successive terms in eq(B.22) - directrly
    for (i = n_min + 2; i <= nMax + 1; i++) { // for loop
      ii = i - 1;
      d_temp1 =
          std::sqrt((double(ii) + 1.) * (double(ii) + 1.) - double(l * l));
      d_temp2 =
          std::sqrt((double(ii) + 1.) * (double(ii) + 1.) - double(m_ * m_));
      d_temp3 = (2. * double(ii) + 1.) * (ii * (ii + 1.) * vig_x - l * m_);
      d_temp4 = (double(ii) + 1.);
      d_temp5 = std::sqrt(double(ii * ii) - double(l * l));
      d_temp6 = std::sqrt(double(ii * ii) - double(m_ * m_));
      // evaluate from above
      // vig_d        : VIG_d[i] = ...*VIG_d[i-1] + ...*VIG_d[i-2] - eqn(B.22)
      d_temp = (1. / double(ii) / d_temp1 / d_temp2) *
               ((d_temp3)*Wigner[i - 1] -
                (d_temp4 * d_temp5 * d_temp6) * Wigner[i - 2]);
      if (i <= nMax)
        Wigner[i] = d_temp;
      // d_vig_d      : d_VIG_d[i-1] = ...*VIG_d[i-2] + ...*VIG_d[i-1] +
      // ...*VIG_d[i]         - eqn(B.26)
      dWigner[i - 1] = -d_temp4 * d_temp5 * d_temp6 /
                           (double(ii) * (2. * double(ii) + 1.)) *
                           Wigner[i - 2] -
                       0. // since l==0
                       +
                       double(ii) * d_temp1 * d_temp2 /
                           (d_temp4 * (2. * double(ii) + 1.)) * d_temp;
      dWigner[i - 1] /= sin(vig_the);
    } // end for loop
  }   // end singularity check

  // III.2 - else if (m!=0), no special case
  // ----------------------------------------
  else { // non-singular case : if (m!=0)
         // A - if (i<n_min), then set VIG_d to zero
         // -----------------------------------
    //              for(i=0; i<n_min; i++){                                 //
    //              for loop
    //              VIG_d[i]  =0.;                                          //
    //              == zero by default
    //      d_VIG_d[i]=0.;                                          // == zero
    //      by default
    //              } // end for loop
    // B - obtain the readily availble VIG_d[n_min] value
    // ------------------------
    // vig_d
    // ---------------------------------------------------------------------
    Wigner[n_min] = vig_d_n_min; // by definition
    // d_vig_d for the previous from current : d_VIG_d[i] = 0.*VIG_d[i-1] +
    // 0.*VIG_d[i] + ...*VIG_d[i+1]    - eqn(B.26)
    ii = n_min - 1;
    d_temp1 = std::sqrt((double(ii) + 1.) * (double(ii) + 1.) - double(l * l));
    d_temp2 =
        std::sqrt((double(ii) + 1.) * (double(ii) + 1.) - double(m_ * m_));
    d_temp3 = (2. * double(ii) + 1.);
    d_temp4 = (double(ii) + 1.);
    dWigner[ii] =
        double(ii) * d_temp1 * d_temp2 / (d_temp3 * d_temp4) * Wigner[n_min];
    dWigner[ii] /= sin(vig_the); // == zero by definition
    // C - obtain all other values in VIG_d[i] from recursive relationship
    // --------
    //   - Based on 'n_min' value and 'n' value; total recursive steps == n -
    //   n_min
    for (i = n_min + 1; i <= nMax + 1; i++) { // for loop
      ii = i - 1;
      d_temp1 =
          std::sqrt((double(ii) + 1.) * (double(ii) + 1.) - double(l * l));
      d_temp2 =
          std::sqrt((double(ii) + 1.) * (double(ii) + 1.) - double(m_ * m_));
      d_temp3 = (2. * double(ii) + 1.) * (ii * (ii + 1.) * vig_x - l * m_);
      d_temp4 = (double(ii) + 1.);
      d_temp5 = std::sqrt(double(ii * ii) - double(l * l));
      d_temp6 = std::sqrt(double(ii * ii) - double(m_ * m_));
      // evaluate from terms above
      d_temp = (1. / double(ii) / d_temp1 / d_temp2) *
               ((d_temp3)*Wigner[i - 1] -
                (d_temp4 * d_temp5 * d_temp6) * Wigner[i - 2]);
      // populate VIG_d[i] array
      if (i <= nMax)
        Wigner[i] = d_temp;
      // d_vig_d for the previous from current : d_VIG_d[i-1] = ...*VIG_d[i-2] +
      // ...*VIG_d[i-1] + ...*VIG_d[i]        - eqn(B.26)
      dWigner[ii] = -d_temp4 * d_temp5 * d_temp6 /
                        (double(ii) * (2. * double(ii) + 1.)) * Wigner[i - 2] -
                    0. // since l==0
                    +
                    double(ii) * d_temp1 * d_temp2 /
                        (d_temp4 * (2. * double(ii) + 1.)) * d_temp;
      dWigner[ii] /= sin(vig_the);
    } // end : for loop
  }   // end : non-singular case : if (m!=0)
  // -------------------------------------------------------------------------------

  // IV - if (m<0) : apply symmetry property eq(B.7) to eqs(B.22-B.24)
  // -------------
  if (check_m_negative == 1) {    // solve using symmetry relation eq(B.7)
    for (i = 0; i <= nMax; i++) { // for loop
      d_temp1 =
          1. / (pow(-1., double(i + l))); // this can be simplified, since l==0
      Wigner[i] *= d_temp1;               // obtain final VIG_d
      dWigner[i] *= -d_temp1;             // obtain final VIG_d
    }                                     // end : for loop
  }                                       // end : if (m<0) condition

  // ------------------------------------------------------------------------------

  return 0;
}

AuxCoefficients::AuxCoefficients(Spherical<double> R_,
                                 std::complex<double> waveK_, int regular_,
                                 int nMax_) {
  R = R_;
  waveK = waveK_;
  nMax = nMax_;
  besselType = regular_;

  if (regular_ == 0)
    besselType = 1;
  else
    besselType = 0;

  dataMp = new SphericalP<std::complex<double>>[Tools::iteratorMax(nMax)];
  dataNp = new SphericalP<std::complex<double>>[Tools::iteratorMax(nMax)];
  dataBp = new SphericalP<std::complex<double>>[Tools::iteratorMax(nMax)];
  dataCp = new SphericalP<std::complex<double>>[Tools::iteratorMax(nMax)];

  dn = new double[nMax + 1];

  compute_MpNp(R, waveK, besselType, nMax, dataMp, dataNp);
}
