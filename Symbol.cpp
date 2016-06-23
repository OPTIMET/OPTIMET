#include "Symbol.h"

#include <cmath>
#include "constants.h"
#include "Bessel.h"
#include "CompoundIterator.h"
#include "gsl/gsl_sf_coupling.h"

namespace optimet {
namespace symbol {

double Wigner3j(int j1, int j2, int j3, int m1, int m2, int m3) {
  return gsl_sf_coupling_3j(2 * j1, 2 * j2, 2 * j3, 2 * m1, 2 * m2, 2 * m3);
}

namespace {

double Wigner6j(int j1, int j2, int j3, int j4, int j5, int j6) {
  return gsl_sf_coupling_6j(2 * j1, 2 * j2, 2 * j3, 2 * j4, 2 * j5, 2 * j6);
}

double Wigner9j(int j11, int j12, int j13, int j21, int j22, int j23, int j31,
                int j32, int j33) {
  return gsl_sf_coupling_9j(2 * j11, 2 * j12, 2 * j13, 2 * j21, 2 * j22,
                            2 * j23, 2 * j31, 2 * j32, 2 * j33);
}

double CleGor(int j, int m, int j1, int m1, int j2, int m2) {
  return std::pow(-1.0, m + j1 - j2) * std::sqrt(2.0 * j + 1.0) *
         Wigner3j(j1, j2, j, m1, m2, -m);
}

// C numbers
double C_01m1(int J1, int M1, int J2, int M2, int J, int M) {
  return std::sqrt(3.0 / 2.0 / consPi) * CleGor(J, M, J1, M1, J2, M2) *
         (std::sqrt((J1 + 1.0) * J2 * (2.0 * J1 - 1.0) * (2.0 * J2 - 1.0)) *
              Wigner9j(J1, J1 - 1, 1, J2, J2 - 1, 1, J, J, 1) *
              CleGor(J, 0, J1 - 1, 0, J2 - 1, 0) -
          std::sqrt((J1 + 1.0) * (J2 + 1.0) * (2.0 * J1 - 1.0) *
                    (2.0 * J2 + 3.0)) *
              Wigner9j(J1, J1 - 1, 1, J2, J2 + 1, 1, J, J, 1) *
              CleGor(J, 0, J1 - 1, 0, J2 + 1, 0) +
          std::sqrt(J1 * J2 * (2.0 * J1 + 3.0) * (2.0 * J2 - 1.0)) *
              Wigner9j(J1, J1 + 1, 1, J2, J2 - 1, 1, J, J, 1) *
              CleGor(J, 0, J1 + 1, 0, J2 - 1, 0) -
          std::sqrt(J1 * (J2 + 1) * (2.0 * J1 + 3.0) * (2.0 * J2 + 1.0)) *
              Wigner9j(J1, J1 + 1, 1, J2, J2 + 1, 1, J, J, 1) *
              CleGor(J, 0, J1 + 1, 0, J2 + 1, 0));
}

double C_10m1(int J1, int M1, int J2, int M2, int J, int M) {
  return std::sqrt(3.0 / 2.0 / consPi) * (2.0 * J1 + 1.0) *
         CleGor(J, M, J1, M1, J2, M2) *
         (std::sqrt(J2 * (2.0 * J1 - 1.0)) *
              Wigner9j(J1, J1, 1, J2, J2 - 1, 1, J, J + 1, 1) *
              CleGor(J + 1, 0, J1, 0, J2 - 1, 0) *
              std::sqrt(J / (2.0 * J + 1.0)) -
          std::sqrt((J2 + 1.0) * (2.0 * J2 + 3.0)) *
              Wigner9j(J1, J1, 1, J2, J2 + 1, 1, J, J + 1, 1) *
              CleGor(J + 1, 0, J1, 0, J2 + 1, 0) *
              std::sqrt(J / (2.0 * J + 1.0)) +
          std::sqrt(J2 * (2.0 * J2 - 1.0)) *
              Wigner9j(J1, J1, 1, J2, J2 - 1, 1, J, J - 1, 1) *
              CleGor(J - 1, 0, J1, 0, J2 - 1, 0) *
              std::sqrt((J + 1.0) / (2.0 * J + 1.0)) -
          std::sqrt((J2 + 1) * (2.0 * J1 + 3.0)) *
              Wigner9j(J1, J1, 1, J2, J2 + 1, 1, J, J - 1, 1) *
              CleGor(J - 1, 0, J1, 0, J2 + 1, 0) *
              std::sqrt((J + 1) / (2.0 * J + 1.0)));
}

double C_00m1(int J1, int M1, int J2, int M2, int J, int M) {
  return std::sqrt(3.0 / 2.0 / consPi) * (2.0 * J1 + 1.0) *
         CleGor(J, M, J1, M1, J2, M2) *
         (std::sqrt(J2 * (2.0 * J1 - 1.0)) *
              Wigner9j(J1, J1, 1, J2, J2 - 1, 1, J, J, 1) *
              CleGor(J, 0, J1, 0, J2 - 1, 0) -
          std::sqrt((J2 + 1.0) * (2.0 * J2 + 3.0)) *
              Wigner9j(J1, J1, 1, J2, J2 + 1, 1, J, J, 1) *
              CleGor(J + 1, 0, J1, 0, J2 + 1, 0));
}

double C_11m1(int J1, int M1, int J2, int M2, int J, int M) {
  return std::sqrt(3.0 / 2.0 / consPi) * CleGor(J, M, J1, M1, J2, M2) *
         (std::sqrt((J1 + 1.0) * J2 * (2.0 * J1 - 1.0) * (2.0 * J2 - 1.0)) *
              Wigner9j(J1, J1 - 1, 1, J2, J2 - 1, 1, J, J + 1, 1) *
              CleGor(J + 1, 0, J1 - 1, 0, J2 - 1, 0) *
              std::sqrt(J / (2.0 * J + 1.0)) -
          std::sqrt((J1 + 1.0) * (J2 + 1) * (2.0 * J1 - 1.0) *
                    (2.0 * J2 + 3.0)) *
              Wigner9j(J1, J1 - 1, 1, J2, J2 + 1, 1, J, J + 1, 1) *
              CleGor(J + 1, 0, J1 - 1, 0, J2 + 1, 0) *
              std::sqrt(J / (2.0 * J + 1.0)) +
          sqrt(J1 * J2 * (2.0 * J1 + 3.0) * (2.0 * J2 - 1.0)) *
              Wigner9j(J1, J1 + 1, 1, J2, J2 - 1, 1, J, J + 1, 1) *
              CleGor(J + 1, 0, J1 + 1, 0, J2 - 1, 0) *
              std::sqrt(J / (2.0 * J + 1.0)) -
          std::sqrt(J1 * (J2 + 1) * (2.0 * J1 + 3.0) * (2.0 * J2 + 3.0)) *
              Wigner9j(J1, J1 + 1, 1, J2, J2 + 1, 1, J, J + 1, 1) *
              CleGor(J + 1, 0, J1 + 1, 0, J2 + 1, 0) *
              std::sqrt(J / (2.0 * J + 1.0)) +
          std::sqrt((J1 + 1.0) * J2 * (2.0 * J1 - 1.0) * (2.0 * J2 - 1.0)) *
              Wigner9j(J1, J1 - 1, 1, J2, J2 - 1, 1, J, J - 1, 1) *
              CleGor(J - 1, 0, J1 - 1, 0, J2 - 1, 0) *
              std::sqrt((J + 1) / (2.0 * J + 1.0)) -
          std::sqrt((J1 + 1.0) * (J2 + 1) * (2.0 * J1 - 1.0) *
                    (2.0 * J2 + 3.0)) *
              Wigner9j(J1, J1 - 1, 1, J2, J2 + 1, 1, J, J - 1, 1) *
              CleGor(J - 1, 0, J1 - 1, 0, J2 + 1, 0) *
              std::sqrt((J + 1) / (2.0 * J + 1.0)) +
          std::sqrt(J1 * J2 * (2.0 * J1 + 3.0) * (2.0 * J2 - 1.0)) *
              Wigner9j(J1, J1 + 1, 1, J2, J2 - 1, 1, J, J - 1, 1) *
              CleGor(J - 1, 0, J1 + 1, 0, J2 - 1, 0) *
              std::sqrt((J + 1) / (2.0 * J + 1.0)) -
          std::sqrt(J1 * (J2 + 1) * (2.0 * J1 + 3.0) * (2.0 * J2 + 3.0)) *
              Wigner9j(J1, J1 + 1, 1, J2, J2 + 1, 1, J, J - 1, 1) *
              CleGor(J - 1, 0, J1 + 1, 0, J2 + 1, 0) *
              std::sqrt((J + 1) / (2.0 * J + 1.0)));
}

// A numbers
std::complex<double> A_0(int n, double R, const std::complex<double> &waveK_i,
                         const std::complex<double> &cmn_1) {
  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = bessel<Hankel1, false>(R * waveK_i, n);

  return data[n] * cmn_1;
}

std::complex<double> A_1(int n, double R, const std::complex<double> &waveK_i,
                         const std::complex<double> &dmn_1) {
  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = bessel<Hankel1, false>(R * waveK_i, n);
  // to be double checked - d/dr(J(kr))
  return std::complex<double>(0.0, 1.0) * (1.0 / waveK_i) *
         (waveK_i * ddata[n] + data[n] / R) * dmn_1;
}

std::complex<double> A_m1(int n, double R, const std::complex<double> &waveK_i,
                          const std::complex<double> &dmn_1) {
  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = bessel<Hankel1, false>(R * waveK_i, n);

  return std::complex<double>(0.0, 1.0) * std::sqrt(n * (n + 1.0)) *
         (1.0 / waveK_i / R) * data[n] * dmn_1;
}

double W(int L1, int J1, int M1, int L2, int J2, int M2, int L, int M) {
  return std::pow(-1.0, J2 + L1 + L) *
         std::sqrt((2.0 * J1 + 1.0) * (2.0 * J2 + 1.0) * (2.0 * L1 + 1.0) *
                   (2.0 * L2 + 1.0) / (4 / consPi / (2.0 * L + 1))) *
         Wigner6j(L1, L2, L, J2, J1, 1) * CleGor(L, 0, L1, 0, L2, 0) *
         CleGor(L, M, J1, M1, J2, M2);
}

} // namespace

// In what follows, need to devise CompoundIterator to start from n=0.
// otherwise, use double for loops

std::complex<double> up_mn(int m, int n, int nMax,
                           const std::complex<double> &cmn_1,
                           const std::complex<double> &dmn_1, double omega,
                           const Scatterer &object,
                           const ElectroMagnetic &bground) {
  // Basic relations
  const double R = object.radius;
  const std::complex<double> mu_0 = consMu0;
  const std::complex<double> mu_b = bground.mu;
  const std::complex<double> eps_b = bground.epsilon;
  // SH Basic relations
  const std::complex<double> b_non = object.elmag.b_SH;

  // Auxiliary variables
  const std::complex<double> waveK_j1 =
      (omega / 2.0) * std::sqrt(object.elmag.epsilon * object.elmag.mu);
  std::complex<double> sum(0.0, 0.0);

  for (CompoundIterator p(nMax, nMax); p < p.max(nMax); ++p) {
    for (CompoundIterator q(nMax, nMax); q < q.max(nMax); ++q)
      sum +=
          A_1(p.first, R, waveK_j1, dmn_1) * A_m1(q.first, R, waveK_j1, dmn_1) *
              C_10m1(p.first, p.second, q.first, q.second, n, m) +
          A_0(p.first, R, waveK_j1, cmn_1) * A_m1(q.first, R, waveK_j1, dmn_1) *
              C_11m1(p.first, p.second, q.first, q.second, n, m);
  }

  return sum * std::complex<double>(0.0, 1.0) * (-b_non / 2.0) *
         std::sqrt(mu_b / eps_b) / std::sqrt(mu_0 / eps_b);
}

std::complex<double> vp_mn(int m, int n, int nMax,
                           const std::complex<double> &cmn_1,
                           const std::complex<double> &dmn_1, double omega,
                           const Scatterer &object,
                           const ElectroMagnetic &bground) {
  // Basic relations
  const double R = object.radius;
  const std::complex<double> mu_0 = consMu0;
  const std::complex<double> eps_0 = consEpsilon0;
  const std::complex<double> mu_b = bground.mu;
  const std::complex<double> eps_b = bground.epsilon;
  // SH Basic relations
  const std::complex<double> b_non = object.elmag.b_SH;

  // Auxiliary variables
  const std::complex<double> waveK_j1 =
      (omega / 2.0) * std::sqrt(object.elmag.epsilon * object.elmag.mu);

  std::complex<double> sum(0.0, 0.0);
  for (CompoundIterator p(nMax, nMax); p < p.max(nMax); ++p) {
    for (CompoundIterator q(nMax, nMax); q < q.max(nMax); ++q)
      sum +=
          A_1(p.first, R, waveK_j1, dmn_1) * A_m1(q.first, R, waveK_j1, dmn_1) *
              C_00m1(p.first, p.second, q.first, q.second, n, m) +
          A_0(p.first, R, waveK_j1, cmn_1) * A_m1(q.first, R, waveK_j1, dmn_1) *
              C_01m1(p.first, p.second, q.first, q.second, n, m);
  }
  return sum * std::complex<double>(-2.0, 0.0) * (-b_non / 2.0) *
         std::sqrt(mu_b / eps_b) / std::sqrt(mu_0 / eps_0);
}

std::complex<double> upp_mn(int m, int n, int nMax,
                            const std::complex<double> &cmn_1,
                            const std::complex<double> &dmn_1, double omega,
                            const Scatterer &object) {

  // Basic relations
  const double R = object.radius;
  const std::complex<double> eps_0 = consEpsilon0;
  const std::complex<double> mu_j = object.elmag.mu;
  const std::complex<double> eps_j = object.elmag.epsilon;
  // SH Basic relations
  const std::complex<double> a_non = object.elmag.a_SH;
  const std::complex<double> d_non = object.elmag.d_SH;
  const std::complex<double> eps_j2 = object.elmag.epsilon_SH;

  // Auxiliary variables
  const std::complex<double> waveK_01 = (omega / 2.0) * std::sqrt(eps_j * mu_j);
  const std::complex<double> waveK_j1 = (omega / 2.0) * std::sqrt(eps_j * mu_j);

  std::complex<double> gmn(0.0, 0.0);
  std::complex<double> fmn(0.0, 0.0);

  for (CompoundIterator p(nMax, nMax); p < p.max(nMax); ++p) {
    for (CompoundIterator q(nMax, nMax); q < q.max(nMax); ++q) {
      gmn += A_m1(p.first, R, waveK_j1, dmn_1) *
             A_m1(q.first, R, waveK_j1, dmn_1) *
             (std::sqrt(p.first / (2.0 * p.first + 1.0)) *
                  std::sqrt(q.first / (2.0 * q.first + 1.0)) *
                  W(p.first - 1, p.first, p.second, q.first - 1, q.first,
                    q.second, n, m) +
              std::sqrt((p.first + 1) / (2.0 * p.first + 1.0)) *
                  std::sqrt((q.first + 1) / (2.0 * q.first + 1.0)) *
                  W(p.first + 1, p.first, p.second, q.first + 1, q.first,
                    q.second, n, m) -
              std::sqrt(p.first / (2.0 * p.first + 1.0)) *
                  std::sqrt((q.first + 1) / (2.0 * q.first + 1.0)) *
                  W(p.first - 1, p.first, p.second, q.first + 1, q.first,
                    q.second, n, m) -
              std::sqrt((p.first + 1) / (2.0 * p.first + 1.0)) *
                  std::sqrt(q.first / (2.0 * q.first + 1.0)) *
                  W(p.first + 1, p.first, p.second, q.first - 1, q.first,
                    q.second, n, m));

      fmn +=
          A_1(p.first, R, waveK_j1, dmn_1) * A_1(q.first, R, waveK_j1, dmn_1) *
          (std::sqrt((p.first + 1) / (2.0 * p.first + 1.0)) *
               std::sqrt((q.first + 1) / (2.0 * q.first + 1.0)) *
               W(p.first - 1, p.first, p.second, q.first - 1, q.first, q.second,
                 n, m) +
           std::sqrt(p.first / (2.0 * p.first + 1.0)) *
               std::sqrt(q.first / (2.0 * q.first + 1.0)) *
               W(p.first + 1, p.first, p.second, q.first + 1, q.first, q.second,
                 n, m) -
           std::sqrt((p.first + 1) / (2.0 * p.first + 1.0)) *
               std::sqrt(q.first / (2.0 * q.first + 1.0)) *
               W(p.first - 1, p.first, p.second, q.first + 1, q.first, q.second,
                 n, m) -
           std::sqrt(p.first / (2.0 * p.first + 1.0)) *
               std::sqrt((q.first + 1) / (2.0 * q.first + 1.0)) *
               W(p.first + 1, p.first, p.second, q.first - 1, q.first, q.second,
                 n, m) +
           A_0(p.first, R, waveK_j1, cmn_1) * A_0(q.first, R, waveK_j1, cmn_1) *
               W(p.first, p.first, p.second, q.first, q.first, q.second, n, m) +
           A_1(p.first, R, waveK_j1, dmn_1) * A_0(q.first, R, waveK_j1, cmn_1) *
               (std::sqrt((p.first + 1) / (2.0 * p.first + 1.0)) *
                    W(p.first - 1, p.first, p.second, q.first, q.first,
                      q.second, n, m) +
                std::sqrt(p.first / (2.0 * p.first + 1.0)) *
                    W(p.first + 1, p.first, p.second, q.first, q.first,
                      q.second, n, m)) +
           A_0(p.first, R, waveK_j1, cmn_1) * A_1(q.first, R, waveK_j1, dmn_1) *
               (std::sqrt((p.first + 1) / (2.0 * p.first + 1.0)) *
                    W(p.first, p.first, p.second, q.first - 1, q.first,
                      q.second, n, m) +
                std::sqrt(q.first / (2.0 * p.first + 1.0)) *
                    W(p.first, p.first, p.second, q.first + 1, q.first,
                      q.second, n, m)));
    }
  }

  return std::complex<double>(0.0, 1.0) * (-a_non / 4.0) *
             std::sqrt(n * (n + 1)) * gmn / waveK_01 / R +
         std::complex<double>(0.0, 1.0) * (-d_non / 8.0) * (eps_0 / eps_j2) *
             std::sqrt(n * (n + 1)) * (gmn + fmn) / waveK_01 / R;
}

} // namespace symbol
} // namespace optimet
