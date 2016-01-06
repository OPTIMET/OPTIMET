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
  return std::make_tuple(
      0.5 / std::sqrt(static_cast<t_real>(l * (l + 1) * n * (n + 1))),
      std::sqrt(static_cast<t_real>((n - m) * (n + m + 1) * (l - k) * (l + k + 1))),
      std::sqrt(static_cast<t_real>((n + m) * (n - m + 1) * (l + k) * (l - k + 1))));
}
t_complex coefficients_A(t_int l, t_int n, t_int m, t_int k, t_int n_max, t_complex ****AlBe_nmlk) {
  if(std::abs(k) > l)
    return 0e0;
  auto const offset = 7; // something optimet
  auto const coeffs = coefficients_A(l, n, m, k);
  auto const zero = AlBe_nmlk[n][m + n_max + offset][l][k + n_max + offset];
  auto const one = AlBe_nmlk[n][m + n_max + offset + 1][l][k + n_max + offset + 1];
  auto const two = AlBe_nmlk[n][m + n_max + offset - 1][l][k + n_max + offset - 1];
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
  return factor * (c0 * ta(n, m, l, k) + c1 * ta(n, m + 1, l, k + 1) + c2 * ta(n, m - 1, l, k - 1));
}

std::tuple<t_real, t_real, t_real> coefficients_B(t_int l, t_int n, t_int m, t_int k) {
  t_real const a0 = 2 * l + 1;
  t_real const a1 = (2 * l - 1) * l * (l + 1) * n * (n + 1);
  return std::make_tuple(0.5 * std::sqrt(a0 / a1),
                         std::sqrt((n - m) * (n + m + 1) * (l - k) * (l - k - 1)),
                         std::sqrt((n + m) * (n - m + 1) * (l + k) * (l + k - 1)));
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
         (static_cast<t_real>(2 * m) * std::sqrt(l * l - k * k) * zero + std::get<1>(coeffs) * one -
          std::get<2>(coeffs) * two);
}
t_complex coefficients_B(t_int n, t_int m, t_int l, t_int k, TranslationAdditionCoefficients &ta) {
  if(std::abs(k) > l)
    return 0e0;
  t_real const a0 = 2 * l + 1;
  t_real const a1 = (2 * l - 1) * l * (l + 1) * n * (n + 1);
  t_complex const factor(0, -0.5 * std::sqrt(a0 / a1));
  auto const c0 = static_cast<t_real>(2 * m) * std::sqrt((l - k) * (l + k));
  auto const c1 = std::sqrt((n - m) * (n + m + 1) * (l - k) * (l - k - 1));
  auto const c2 = std::sqrt((n + m) * (n - m + 1) * (l + k) * (l + k - 1));
  return factor * (c0 * ta(n, m, l - 1, k) + c1 * ta(n, m + 1, l - 1, k + 1) -
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
  return std::make_tuple(diagonal, offdiagonal);
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
  TranslationAdditionCoefficients ta(R, waveK, BHreg == 0);

  for(int n = 1; n <= n_max; n++)  // start at n==1 up to and including n_max
    for(int m = -n; m <= n; m++) { // increment by the padded e value

      auto const p = flatten_indices(n, m);

      for(int l = 1; l <= n_max; l++)
        for(int k = -l; k <= l; k++) {
          auto const q = flatten_indices(l, k);
          dataApq(p, q) = coefficients_A(l, n, m, k, n_max, AlBe_nmlk);
          dataBpq(p, q) = coefficients_B(l, n, m, k, n_max, AlBe_nmlk);
        }
    }

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
  if(std::abs(relR.rrr) < errEpsilon) { // Check for NO translation case
    offdiagonal = Matrix<t_complex>::Zero(n, n);
    diagonal = Matrix<t_complex>::Identity(n, n);
  } else
    std::tie(diagonal, offdiagonal) = transfer_coefficients(relR, waveK, not regular, nMax);
}
} // namespace optimet
