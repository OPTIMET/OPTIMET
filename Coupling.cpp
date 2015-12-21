#include "Coupling.h"

#include "PeriodicCoupling.h"
#include "CompoundIterator.h"
#include "constants.h"
#include "TranslationAdditionCoefficients.h"

#include <cmath>
#include <iostream>

namespace optimet {

namespace {
std::tuple<t_real, t_real, t_real> coefficients_A(t_int l, t_int n, t_int m, t_int k) {
  return {0.5 / std::sqrt(static_cast<t_real>(l * (l + 1) * n * (n + 1))),
          std::sqrt(static_cast<t_real>((n - m) * (n + m + 1) * (l - k) * (l + k + 1))),
          std::sqrt(static_cast<t_real>((n + m) * (n - m + 1) * (l + k) * (l - k + 1)))};
}
t_complex coefficients_A(t_int l, t_int n, t_int m, t_int k, t_int n_max, t_complex ****AlBe_nmlk) {
  if(std::abs(k) > l)
    return 0e0;
  auto const offset = 7; // something optimet
  auto const coeffs = coefficients_A(l, n, m, k);
  auto const zero = AlBe_nmlk[n][m + n_max + offset][l][k + n_max + offset];
  auto const one = AlBe_nmlk[n][m + n_max + offset + 1][l][k + n_max + offset + 1];
  auto const two = AlBe_nmlk[n][m + n_max + offset - 1][l][k + n_max + offset - 1];
  if(n == 1 and l == 1 and m == 1 and k == 1)
    std::cout << "old: " << zero << " " << one << " " << two << "\n";
  return std::get<0>(coeffs) * (static_cast<t_real>(2 * k * m) * zero + std::get<1>(coeffs) * one +
                                std::get<2>(coeffs) * two);
}
t_complex coefficients_A(t_int n, t_int m, t_int l, t_int k, TranslationAdditionCoefficients &ta) {
  if(std::abs(k) > l)
    return 0e0;
  auto const factor = 0.5 / std::sqrt(static_cast<t_real>(l * (l + 1) * n * (n + 1)));
  auto const c0 = static_cast<t_real>(2 * k * m);
  auto const c1 = std::sqrt(static_cast<t_real>((n - m) * (n + m + 1) * (l - k) * (l + k + 1)));
  auto const c2 = std::sqrt(static_cast<t_real>((n + m) * (n - m + 1) * (l + k) * (l - k + 1)));
  if(n == 1 and l == 1 and m == 1 and k == 1)
    std::cout << "new: " << ta(n, m, l, k) << " " << ta(n, m + 1, l, k + 1) << " "
              << ta(n, m - 1, l, k - 1) << "\n";
  return factor * (c0 * ta(n, m, l, k) + c1 * ta(n, m + 1, l, k + 1) + c2 * ta(n, m - 1, l, k - 1));
}

std::tuple<t_real, t_real, t_real> coefficients_B(t_int l, t_int n, t_int m, t_int k) {
  auto const a0 = 2 * l + 1;
  auto const a1 = (2 * l - 1) * l * (l + 1) * n * (n + 1);
  return {0.5 * std::sqrt(static_cast<t_real>(a0) / static_cast<t_real>(a1)),
          std::sqrt(static_cast<t_real>((n - m) * (n + m + 1) * (l - k) * (l - k - 1))),
          std::sqrt(static_cast<t_real>((n + m) * (n - m + 1) * (l + k) * (l + k - 1)))};
}
t_complex coefficients_B(t_int l, t_int n, t_int m, t_int k, t_int n_max, t_complex ****AlBe_nmlk) {
  if(std::abs(k) > l)
    return 0e0;
  auto const offset = 7; // something optimet
  auto const coeffs = coefficients_B(l, n, m, k);
  auto const zero = AlBe_nmlk[n][m + n_max + offset][l - 1][k + n_max + offset];
  auto const one = AlBe_nmlk[n][m + n_max + offset + 1][l - 1][k + n_max + offset + 1];
  auto const two = AlBe_nmlk[n][m + n_max + offset - 1][l - 1][k + n_max + offset - 1];
  return t_complex(0, -std::get<0>(coeffs)) *
         (static_cast<t_real>(2 * m) * std::sqrt(static_cast<t_real>(l * l - k * k)) * zero +
          std::get<1>(coeffs) * one - std::get<2>(coeffs) * two);
}
t_complex coefficients_B(t_int n, t_int m, t_int l, t_int k, TranslationAdditionCoefficients &ta) {
  if(std::abs(k) > l)
    return 0e0;
  t_real const a0 = 2 * l + 1;
  t_real const a1 = (2 * l - 1) * l * (l + 1) * n * (n + 1);
  auto const factor = -constant::i * 0.5 / std::sqrt(a0 / a1);
  auto const c0 = static_cast<t_real>(2 * m) * std::sqrt(static_cast<t_real>((l - k) * (l + k)));
  auto const c1 = std::sqrt(static_cast<t_real>((n - m) * (n + m + 1) * (l - k) * (l - k - 1)));
  auto const c2 = std::sqrt(static_cast<t_real>((n + m) * (n - m + 1) * (l + k) * (l + k - 1)));
  return factor * (c0 * ta(n, m, l - 1, k) + c1 * ta(n, m + 1, l - 1, k + 1) +
                   c2 * ta(n, m - 1, l - 1, k - 1));
}

std::tuple<Matrix<t_complex>, Matrix<t_complex>>
transfer_coefficients(Spherical<double> R, std::complex<double> waveK, bool regular, int n_max) {
  auto const N = Tools::iteratorMax(n_max);
  Matrix<t_complex> diagonal = Matrix<t_complex>::Zero(N, N);
  Matrix<t_complex> offdiagonal = Matrix<t_complex>::Zero(N, N);

  TranslationAdditionCoefficients ta(R, waveK, regular);

  // start at harmonic n = 1. (because n=0 spherical and hence symmetrically incompatible with
  // propagating wave?)
  for(t_int n(1); n <= n_max; ++n)
    for(t_int m(-n); m <= n; ++m) {

      auto const p = flatten_indices(n, m);
      for(t_int l(1); l <= n_max; ++l)
        for(t_int k(-l); k <= l; ++k) {

          auto const q = flatten_indices(l, k);
          diagonal(p, q) = coefficients_A(n, m, l, k, ta);
          offdiagonal(p, q) = coefficients_B(n, m, l, k, ta);
        }
    }
  return {diagonal, offdiagonal};
}

/**
 * Calculates the coupling coefficients for a relative vector R.
 * @param R the relative Spherical vector.
 * @param waveK the complex wave vector.
 * @param BHreg the type of Bessel function.
 * @param n_max the maximum value of the n iterator.
 * @param dataApq the storage vector for dataApq.
 * @param dataBpq the storage vector for dataBpq.
 */
void TransferCoefficients(Spherical<double> R, std::complex<double> waveK, int BHreg, int n_max,
                          Matrix<t_complex> &dataApq, Matrix<t_complex> &dataBpq) {
  // prepare for computing transfer matrices
  // ---------------------------------------
  // Matrices sizes
  // ----------------------------------------------------------------
  auto const e = 7;                        // required for extra values in 'AlBe_00lk'
  auto const n_Matsize1 = (n_max + e) + 1; // up to and including (n_max+e)  : indexed from 1
  auto const m_Matsize1 =
      2 * (n_max + e) + 1; // up to and including (n_max+e)  : indexed from 0 + 1 for m==0
  // scalar translation-addition theorem
  // ------------------------------------------
  auto AlBe_nmlk = Tools::Get_4D_c_double(n_Matsize1, m_Matsize1, n_Matsize1, m_Matsize1);
  compute_AlBe_nmlk(R, waveK, BHreg, n_max, AlBe_nmlk);

  for(int n = 1; n < n_Matsize1 - e; n++) { // start at n==1 up to and including n_max
    for(int j = e, m = -n_max; j < m_Matsize1 - e; j++, m++) { // increment by the padded e value

      if(std::abs(m) <= n) {
        auto const p = flatten_indices(n, m);

        for(int l = 1; l < n_Matsize1 - e; l++) {
          for(int jj = e, k = -n_max; jj < m_Matsize1 - e; jj++, k++) {

            if(std::abs(k) <= l) {
              auto const q = flatten_indices(l, k);
              dataApq(p, q) = coefficients_A(l, n, j - n_max - e, jj - n_max - e, n_max, AlBe_nmlk);
              dataBpq(p, q) = coefficients_B(l, n, j - n_max - e, jj - n_max - e, n_max, AlBe_nmlk);
            }

          } // jj
        }   // l
      }     // if(abs(m)<=n)
    }       // j
  }         // n

  // C - delete all other matrices
  // --------------------------------------------------------
  for(int i = 0; i < n_Matsize1; i++)
    for(int j = 0; j < m_Matsize1; j++)
      for(int ii = 0; ii < n_Matsize1; ii++)
        delete[] AlBe_nmlk[i][j][ii];

  for(int i = 0; i < n_Matsize1; i++)
    for(int j = 0; j < m_Matsize1; j++)
      delete[] AlBe_nmlk[i][j];

  for(int i = 0; i < n_Matsize1; i++)
    delete[] AlBe_nmlk[i];

  delete[] AlBe_nmlk;
}
} // anonymous namespace

Coupling::Coupling(Spherical<t_real> relR, t_complex waveK, t_uint nMax, bool regular) {
  auto const n = Tools::iteratorMax(nMax);
  offdiagonal = Matrix<t_complex>::Zero(n, n);
  if(std::abs(relR.rrr) < errEpsilon) { // Check for NO translation case
    diagonal = Matrix<t_complex>::Identity(n, n);
  } else {
    diagonal = Matrix<t_complex>::Zero(n, n);
    TransferCoefficients(relR, waveK, regular ? 1 : 0, nMax, diagonal, offdiagonal);

    // auto const r = transfer_coefficients(relR, waveK, not regular, nMax);
    // std::cout << std::get<0>(r).topLeftCorner(3, 3) << "\n\n" << diagonal.topLeftCorner(3, 3)
    //           << "\n";
    // if(not std::get<0>(r).isApprox(diagonal, 1e-12))
    //   throw std::runtime_error("Goodness gracious");
    // if(not std::get<1>(r).isApprox(offdiagonal, 1e-12))
    //   throw std::runtime_error("Goodness gracious");
  };
}
} // namespace optimet
