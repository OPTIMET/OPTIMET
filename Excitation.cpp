#include "Excitation.h"

#include "Algebra.h"
#include "Tools.h"
#include "constants.h"
#include "AuxCoefficients.h"
#include "CompoundIterator.h"
#include "Coupling.h"

#include <iostream>
#include <cmath>

Excitation::Excitation() { initDone = false; }

Excitation::Excitation(unsigned long type_,
                       SphericalP<std::complex<double>> Einc_,
                       Spherical<double> waveKInc_, int nMax_) {
  init(type_, Einc_, waveKInc_, nMax_);
}

Excitation::~Excitation() {
  if (initDone) {
    delete[] dataIncAp;
    delete[] dataIncBp;
  }
}

void Excitation::init(unsigned long type_,
                      SphericalP<std::complex<double>> Einc_,
                      Spherical<double> vKInc_, int nMax_) {
  type = type_;
  Einc = Einc_;
  vKInc = vKInc_;
  nMax = nMax_;

  waveK = vKInc.rrr;
  lambda = 2 * consPi / waveK.real();
  omega = consC * waveK.real();

  dataIncAp = new std::complex<double>[Tools::iteratorMax(nMax)];
  dataIncBp = new std::complex<double>[Tools::iteratorMax(nMax)];

  initDone = true;
}

void Excitation::update(unsigned long type_,
                        SphericalP<std::complex<double>> Einc_,
                        Spherical<double> vKInc_, int nMax_) {
  type = type_;
  Einc = Einc_;
  vKInc = vKInc_;
  nMax = nMax_;

  waveK = vKInc.rrr;
  lambda = 2 * consPi / waveK.real();
  omega = consC * waveK.real();

  populate();
}

int Excitation::populate() {
  if (!initDone) {
    std::cerr << "Excitation was not initialized!";
    return 1;
  }

  AuxCoefficients coef(Spherical<double>(0.0, vKInc.the, vKInc.phi), waveK, 1,
                       nMax);

  CompoundIterator p;

  for (p = 0; p < p.max(nMax); p++) {
    SphericalP<std::complex<double>> C_local = coef.dataCp[p];
    SphericalP<std::complex<double>> B_local = coef.dataBp[p];

    SphericalP<std::complex<double>> conjAux(
        std::conj(C_local.rrr), std::conj(C_local.the),
        std::conj(C_local.phi)); // std::complex conjugate of C
    dataIncAp[p] = 4 * consPi * std::pow(-1.0, p.second) *
                   std::pow(consCi, p.first) * coef.dn[p.first] *
                   (conjAux * Einc) *
                   std::exp(consCmi * (double)p.second * vKInc.phi);

    conjAux = SphericalP<std::complex<double>>(
        std::conj(B_local.rrr), std::conj(B_local.the),
        std::conj(B_local.phi)); // std::complex conjugate of B
    dataIncBp[p] = 4 * consPi * std::pow(-1.0, p.second) *
                   std::pow(consCi, p.first - 1) * coef.dn[p.first] *
                   (conjAux * Einc) *
                   std::exp(consCmi * (double)p.second * vKInc.phi);
  }

  return 0;
}

int Excitation::getIncLocal(Spherical<double> point_,
                            std::complex<double> *Inc_local_, int nMax_) {
  if (!initDone) {
    std::cerr << "Excitation was not initialized!";
    return 1;
  }

  Spherical<double> Rrel = point_ - Spherical<double>(0.0, 0.0, 0.0);
  optimet::Coupling coupling(Rrel, waveK, 1, nMax_);
  coupling.populate();

  CompoundIterator p, q;

  int pMax = p.max(nMax_);
  int qMax = q.max(nMax_);

  std::complex<double> *Inc_direct = new std::complex<double>[2 * pMax];
  std::complex<double> **T_AB = new std::complex<double> *[2 * (p.max(nMax))];
  for (p = 0; p < (int)(2 * p.max(nMax)); p++) {
    T_AB[p] = new std::complex<double>[2 * p.max(nMax)];
  }

  for (p = 0; p < pMax; p++) {
    Inc_direct[p] = dataIncAp[p];
    Inc_direct[p + pMax] = dataIncBp[p];
  }

  for (p = 0; p < pMax; p++) {
    for (q = 0; q < qMax; q++) {
      T_AB[p][q] = coupling.dataApq(q, p);
      T_AB[p + pMax][q + qMax] = coupling.dataApq(q, p);
      T_AB[p + pMax][q] = coupling.dataBpq(q, p);
      T_AB[p][q + qMax] = coupling.dataBpq(q, p);
    }
  }

  Algebra::multiplyVectorMatrix(T_AB, 2 * pMax, 2 * pMax, Inc_direct,
                                Inc_local_, consC1, consC0);

  delete[] Inc_direct;

  for (p = 0; p < 2 * pMax; p++) {
    delete[] T_AB[p];
  }

  delete[] T_AB;

  return 0;
}

void Excitation::updateWavelength(double lambda_) {
  Spherical<double> vKInc_local = vKInc;
  vKInc_local.rrr = 2 * consPi / lambda_;

  update(type, Einc, vKInc_local, nMax);
}
