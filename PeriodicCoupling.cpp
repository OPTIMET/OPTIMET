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
#include "TranslationAdditionCoefficients.h"

#include <cmath>
#include <cstdlib>
#include <iostream>

namespace optimet {

namespace {
t_real a_nm_p(t_int n, t_int m) {
  t_real const numerator = (n + m + 1) * (n - m + 1);
  t_real const denominator = (2 * n + 1) * (2 * n + 3);
  return -std::sqrt(numerator / denominator);
}

t_real a_nm_m(t_int n, t_int m) {
  t_real const numerator = (n + m) * (n - m);
  t_real const denominator = (2 * n + 1) * (2 * n - 1);
  return std::sqrt(numerator / denominator);
}

t_real b_nm_p(t_int n, t_int m) {
  t_real const numerator = (n + m + 2) * (n + m + 1);
  t_real const denominator = (2 * n + 1) * (2 * n + 3);
  return std::sqrt(numerator / denominator);
}

t_real b_nm_m(t_int n, t_int m) {
  t_real const numerator = (n - m) * (n - m - 1);
  t_real const denominator = (2 * n + 1) * (2 * n - 1);
  return std::sqrt(numerator / denominator);
}
} // anonymous namespace

void compute_AlBe_nmlk(Spherical<t_real> R, t_complex waveK, t_int BHreg, t_int n_max,
                       t_complex ****AlBe_nmlk) {

  // Translation coefficients calculation
  // -----------------------------------------
  // -------------------------------------------------------------------------------
  // temp variables
  t_complex c_temp(0., 0.);
  t_real d_temp(0.);
  // working variables
  t_int i(0), j(0);
  t_int ii(0), jj(0);
  t_int n(0), m(0), l(0), k(0);              // duale subscript variables
  t_real d_n(0.), d_m(0.), d_l(0.), d_k(0.); // duale subscript variables
  // auxiliary variables
  t_complex AlBe2(0., 0.);               // auxiliary coefficients
  t_complex AlBe3(0., 0.);               // auxiliary coefficients
  t_complex AlBe4(0., 0.);               // auxiliary coefficients
  t_complex AlBe2conj(0., 0.);           // auxiliary coefficients
  t_complex AlBe3conj(0., 0.);           // auxiliary coefficients
  t_complex AlBe4conj(0., 0.);           // auxiliary coefficients
  t_real A1(0.), A2(0.), A3(0.), A4(0.); // auxiliary coefficients
  t_real B1(0.), B2(0.), B3(0.);         // auxiliary coefficients
  t_int k_mirror(0), m_mirror(0);        // for calculating (negative)m
  CompoundIterator pl, ql;               // Create a compound iterator

  // Assign corresponding input values
  // ---------------------------------------------
  t_complex wkR = R.rrr * waveK;
  t_complex wkRconj = std::conj(wkR);
  t_complex exp_iwk_theji(0., 0.); // function of m & theji

  // prepare for computing transfer matrices
  // ---------------------------------------
  // Matrices sizes
  // ----------------------------------------------------------------
  // based on n_max
  t_int n_Matsize, m_Matsize; // dependent on (n_max)
  n_Matsize = (n_max) + 1;    // up to and including n_max    : indexed from 1
  // up to and including m_max    : indexed from 0 + 1 for m==0
  m_Matsize = 2 * (n_max) + 1;
  t_int e = 7;                        // required for extra values in 'AlBe_00lk'
  t_int n_Matsize1(0), m_Matsize1(0); // dependent on (n_max+e)
  // up to and including (n_max+e)  : indexed from 1
  n_Matsize1 = (n_max + e) + 1;
  // up to and including (n_max+e)  : indexed from 0 + 1 for m==0
  m_Matsize1 = 2 * (n_max + e) + 1;

  // prepare for calculating translation coefficients
  // -----------------------------
  auto const dataYp = optimet::compute_Yp(R, n_max + e);
  TranslationAdditionCoefficients ta(R, waveK, BHreg == 0);

  // 1.2 Prepare for AlBe_nmlk[0][n_max+e][ii][jj] evaluation
  // ---------------------
  std::vector<t_complex> HB_data, HB_ddata, HBconj_data, HBconj_ddata; // create Bessel object
  try {
    if(BHreg == 0) { // regular case
      // waveK
      std::tie(HB_data, HB_ddata) = optimet::bessel<optimet::Bessel>(wkR, n_Matsize1);
      // conj(waveK)
      std::tie(HBconj_data, HBconj_ddata) = optimet::bessel<optimet::Bessel>(wkRconj, n_Matsize1);
    } else if(BHreg == 1) { // irregular case
      // waveK
      std::tie(HB_data, HB_ddata) = optimet::bessel<optimet::Hankel1>(wkR, n_Matsize1);
      // conj(waveK)
      std::tie(HBconj_data, HBconj_ddata) =
          optimet::bessel<optimet::Hankel1>(consCm1 * wkRconj, n_Matsize1);
    }
  } catch(std::range_error &e) {
    std::cerr << e.what() << std::endl;
  }

  // Building blocks
  // --------------------------------------------------------------
  // I - fundamental building blocks - Ynm
  // ----------------------------------------
  auto Ynm = [&dataYp, &n_max, &e](t_int n, t_int m) {
    return n == 0 and m == n_max + e ? 1e0 / std::sqrt(4e0 * consPi) :
                                       dataYp[flatten_indices(n, m - n_max - e)];
  };

  // II - Global matrix - Alpha-Beta (n, m, l, k) - AlBe_nmlk
  // ---------------------
  t_complex ****AlBe_nmlkconj;
  AlBe_nmlkconj = Tools::Get_4D_c_double(n_Matsize1, m_Matsize1, n_Matsize1, m_Matsize1);
  for(i = 0; i < n_Matsize1; i++) {
    for(j = 0; j < m_Matsize1; j++) {
      for(ii = 0; ii < n_Matsize1; ii++) {
        for(jj = 0; jj < m_Matsize1; jj++) {
          AlBe_nmlk[i][j][ii][jj] = t_complex(0., 0.);
          AlBe_nmlkconj[i][j][ii][jj] = t_complex(0., 0.);
        }
      }
    }
  }

  // A - calculate AlBe_00lk - fundamental building block
  // -------------------------
  n = 0;
  m = 0; // n==m==0
  d_n = t_real(n);
  d_m = t_real(m); // n==m==0
  for(ii = 0, l = 0; ii < n_Matsize1; ii++, l++) {
    for(jj = 0, k = -(n_max + e); jj < m_Matsize1; jj++, k++) {

      d_l = t_real(l);
      d_k = t_real(k);

      // get mirror image in k
      k_mirror = n_max + e + (n_max + e - jj);

      if(abs(k) <= l) {

        // AlBe_00lk[ii][jj] : Bessel/Hankel - Legendre - exp(imthe)
        d_temp = std::sqrt(4. * consPi) * pow(-1., d_l + d_k);
        // wavek
        AlBe_nmlk[0][n_max + e][ii][jj] = d_temp * Ynm(ii, k_mirror) * HB_data[l]; // eqn (C3)
        // waveconj
        AlBe_nmlkconj[0][n_max + e][ii][jj] =
            d_temp * Ynm(ii, k_mirror) * HBconj_data[l]; // eqn (C3)
        if(std::abs(ta(0, 0, l, k) - AlBe_nmlk[0][n_max + e][ii][jj]) > 1e-12) {
          std::cout << l << " " << k << "\n";
          std::cout << ta(0, 0, l, k) << " " << AlBe_nmlk[0][n_max + e][ii][jj] << "\n";
          std::cout << "HB: " << HB_data[l] << " " << HBconj_data[l] << "\n";
          std::cout << optimet::Ynm(R, l, -k) << " " << Ynm(ii, k_mirror) << "\n";
          for(auto const d : dataYp)
            std::cout << "Ynm " << d << "\n";
          std::cout << optimet::Ynm(R, 7, 7) << "\n";
          std::cout << optimet::Ynm(R, 7, 6) << "\n";
          std::cout << optimet::Ynm(R, 7, 5) << "\n";
          std::cout << optimet::Ynm(R, 7, 4) << "\n";
          std::cout << optimet::Ynm(R, 7, 3) << "\n";
          std::cout << optimet::Ynm(R, 7, 2) << "\n";
          std::cout << optimet::Ynm(R, 7, 1) << "\n";
          std::cout << optimet::Ynm(R, 7, 0) << "\n";
          std::cout << optimet::Ynm(R, 7, -1) << "\n";
          std::cout << optimet::Ynm(R, 7, -2) << "\n";
          std::cout << optimet::Ynm(R, 7, -3) << "\n";
          std::cout << optimet::Ynm(R, 7, -4) << "\n";
          std::cout << optimet::Ynm(R, 7, -5) << "\n";
          std::cout << optimet::Ynm(R, 7, -6) << "\n";
          std::cout << optimet::Ynm(R, 7, -7) << "\n";
          throw std::runtime_error("uhoh");
        }

      } // if(abs(k)<=l)

    } // jj
  }   // ii

  // B.1 - calculate AlBe_nmlk - positive(m)==n
  // ----------------------------------
  for(i = 1, n = 1; i < n_Matsize1 - e; i++, n++) { // start at n==1 up to and including n_max
    for(j = e, m = -n_max; j < m_Matsize1 - e; j++, m++) { // increment by the padded e value

      d_n = t_real(n);
      d_m = t_real(m);

      // if n==positive(m) ---------------------------------------------------
      if(n == m) {

        for(ii = 0, l = 0; ii < n_Matsize1 - 1;
            ii++, l++) { // start at n==1 up to and including n_max+e
          for(jj = 0, k = -(n_max + e - 0); jj < m_Matsize1 - 0;
              jj++, k++) { // increment by the padded e value

            d_l = t_real(l);
            d_k = t_real(k);

            if(abs(k) <= l) {
              // obtain three coefficients ----------------------------
              B1 = b_nm_p(n - 1, n - 1);
              B2 = b_nm_p(l - 1, k - 1);
              B3 = b_nm_m(l + 1, k - 1);

              if(((ii - 1) < 0) || ((jj - 1) < 0)) {
                AlBe2 = t_complex(0., 0.);
                AlBe2conj = t_complex(0., 0.);
              } else {
                AlBe2 = AlBe_nmlk[i - 1][j - 1][ii - 1][jj - 1];
                AlBe2conj = AlBe_nmlkconj[i - 1][j - 1][ii - 1][jj - 1];
              }
              //
              if(((ii + 1) > n_Matsize1 - 1) || ((jj - 1) < 0)) {
                AlBe3 = t_complex(0., 0.);
                AlBe3conj = t_complex(0., 0.);
              } else {
                AlBe3 = AlBe_nmlk[i - 1][j - 1][ii + 1][jj - 1];
                AlBe3conj = AlBe_nmlkconj[i - 1][j - 1][ii + 1][jj - 1];
              }

              // evaluate current term --------------------------------
              // wavek
              AlBe_nmlk[i][j][ii][jj] = B2 * AlBe2 + B3 * AlBe3;
              AlBe_nmlk[i][j][ii][jj] /= B1; // eqn (C1-B)
              // wavekconj
              AlBe_nmlkconj[i][j][ii][jj] = B2 * AlBe2conj + B3 * AlBe3conj;
              AlBe_nmlkconj[i][j][ii][jj] /= B1; // eqn (C1-B)
            }
          } // jj
        }   // ii

      } // if(n>0 && n==m)

    } // j
  }   // i

  // B.2 - calculate AlBe_nmlk - for all other positive(m)
  // -----------------------
  for(i = 1, n = 1; i < n_Matsize1 - e; i++, n++) { // start at n==1 up to and including n_max
    for(j = e, m = -n_max; j < m_Matsize1 - e; j++, m++) { // increment by the padded e value

      d_n = t_real(n);
      d_m = t_real(m);

      if(abs(m) <= n) {

        // for all other positive m
        // ------------------------------------------------
        if(m >= 0 && m != n) {

          for(ii = 0, l = 0; ii < n_Matsize1 - 1;
              ii++, l++) { // start at n==1 up to and including n_max+e
            for(jj = 0, k = -(n_max + e - 0); jj < m_Matsize1 - 0;
                jj++, k++) { // increment by the padded e value

              d_l = t_real(l);
              d_k = t_real(k);

              if(abs(k) <= l) {
                // obtain three coefficients -------------------------------
                A1 = a_nm_p(n - 1, m);
                A2 = a_nm_m(n - 1, m);
                A3 = a_nm_p(l - 1, k);
                A4 = a_nm_m(l + 1, k);

                // evaluate current term -----------------------------------
                if(((i - 2) < 0)) {
                  AlBe2 = t_complex(0., 0.);
                  AlBe2conj = t_complex(0., 0.);
                } else {
                  AlBe2 = AlBe_nmlk[i - 2][j][ii][jj];
                  AlBe2conj = AlBe_nmlkconj[i - 2][j][ii][jj];
                }
                //
                if(((ii - 1) < 0)) {
                  AlBe3 = 0.;
                  AlBe3conj = 0.;
                } else {
                  AlBe3 = AlBe_nmlk[i - 1][j][ii - 1][jj];
                  AlBe3conj = AlBe_nmlkconj[i - 1][j][ii - 1][jj];
                }
                //
                if(((ii + 1) > n_Matsize1 - 1)) {
                  AlBe4 = 0.;
                  AlBe4conj = 0.;
                } else {
                  AlBe4 = AlBe_nmlk[i - 1][j][ii + 1][jj];
                  AlBe4conj = AlBe_nmlkconj[i - 1][j][ii + 1][jj];
                }

                // wavek
                AlBe_nmlk[i][j][ii][jj] = -A2 * AlBe2 + A3 * AlBe3 + A4 * AlBe4;
                AlBe_nmlk[i][j][ii][jj] /= A1; // eqn (C1-A)
                // wavek
                AlBe_nmlkconj[i][j][ii][jj] = -A2 * AlBe2conj + A3 * AlBe3conj + A4 * AlBe4conj;
                AlBe_nmlkconj[i][j][ii][jj] /= A1; // eqn (C1-A)
              }
            } // jj
          }   // ii

        } // if(m>=0 && m!=n)

      } // if(m<=n)

    } // j
  }   // i

  // C - calculate AlBe_nmlk - for all other negative(m)
  // --------------------------------
  for(i = 1, n = 1; i < n_Matsize1 - e; i++, n++) { // start at n==1 up to and including n_max
    for(j = e, m = -n_max; j < m_Matsize1 - e; j++, m++) { // increment by the padded e value

      d_n = t_real(n);
      d_m = t_real(m);

      if(abs(m) <= n) {

        // for all other negative m, including n==abs(m)
        // ---------------------------
        if(m < 0) {

          for(ii = 0, l = 0; ii < n_Matsize1 - 1;
              ii++, l++) { // start at n==1 up to and including n_max+e
            for(jj = 0, k = -(n_max + e - 0); jj < m_Matsize1 - 0;
                jj++, k++) { // increment by the padded e value

              d_l = t_real(l);
              d_k = t_real(k);

              if(abs(k) <= l) {

                // get mirror image in m -----------------------------------
                m_mirror = n_max + e + (n_max + e - j);
                // get mirror image in k
                k_mirror = n_max + e + (n_max + e - jj);

                if(BHreg == 0) { // regular    - Bessel
                  // obtain mirror coefficient -------------------------------
                  d_temp = pow(-1., d_k + d_m); // eqn (C4-A)
                  AlBe_nmlk[i][j][ii][jj] = d_temp * conj(AlBe_nmlkconj[i][m_mirror][ii][k_mirror]);
                } else if(BHreg == 1) { // iregular   - Hankel
                  // obtain mirror coefficient -------------------------------
                  d_temp = pow(-1., d_n + d_l + d_k + d_m); // eqn (C4-B)
                  AlBe_nmlk[i][j][ii][jj] =
                      1. * d_temp * conj(AlBe_nmlkconj[i][m_mirror][ii][k_mirror]);
                }
              }
            } // jj
          }   // ii

        } // if(m<0)

      } // if(m<=n)

    } // j
  }   // i

  // C - delete all auxiliary matrices
  // ---------------------------------------------------
  for(i = 0; i < n_Matsize1; i++) {
    for(j = 0; j < m_Matsize1; j++) {
      for(ii = 0; ii < n_Matsize1; ii++) {
        delete[] AlBe_nmlkconj[i][j][ii];
      }
    }
  }

  for(i = 0; i < n_Matsize1; i++) {
    for(j = 0; j < m_Matsize1; j++) {
      delete[] AlBe_nmlkconj[i][j];
    }
  }

  for(i = 0; i < n_Matsize1; i++) {
    delete[] AlBe_nmlkconj[i];
  }

  delete[] AlBe_nmlkconj;
}
}
