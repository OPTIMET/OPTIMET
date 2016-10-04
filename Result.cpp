#include "Result.h"
#include "Coupling.h"
#include "Algebra.h"
#include "constants.h"
#include "Tools.h"
#include "AuxCoefficients.h"

#include <iostream>
#include <fstream>
#include <complex>
#include <cstdlib>

namespace optimet {
Result::Result(Geometry *geometry_, std::shared_ptr<Excitation> excitation_, int nMax_)
  : flagSH(false), result_FF(nullptr) {
  init(geometry_, excitation_, nMax_);
}

Result::Result(Geometry *geometry_, std::shared_ptr<Excitation> excitation_, Result *result_FF_,
               int nMax_) {
  init(geometry_, excitation_, result_FF_, nMax_);
}

void Result::init(Geometry *geometry_, std::shared_ptr<Excitation> excitation_, int nMax_) {
  geometry = geometry_;
  nMax = nMax_;
  excitation = excitation_;
  waveK = excitation->waveK;
  flagSH = false;
  result_FF = NULL;

  scatter_coef.resize(2 * Tools::iteratorMax(nMax) * geometry->objects.size());
  internal_coef.resize(2 * Tools::iteratorMax(nMax) * geometry->objects.size());
  c_scatter_coef.resize(2 * Tools::iteratorMax(nMax));
}

void Result::update(Geometry *geometry_, std::shared_ptr<Excitation> excitation_, int nMax_) {
  geometry = geometry_;
  nMax = nMax_;
  excitation = excitation_;
  waveK = excitation->waveK;
}

void Result::init(Geometry *geometry_, std::shared_ptr<Excitation> excitation_,
                  Result *result_FF_, int nMax_) {
  geometry = geometry_;
  nMax = nMax_;
  excitation = excitation_;
  waveK = excitation->waveK;
  flagSH = true;
  result_FF = result_FF_;

  scatter_coef.resize(2 * Tools::iteratorMax(nMax) * geometry->objects.size());
  internal_coef.resize(2 * Tools::iteratorMax(nMax) * geometry->objects.size());
  c_scatter_coef.resize(2 * Tools::iteratorMax(nMax));
}

void Result::getEHFieldsModal(Spherical<double> R_,
                              SphericalP<std::complex<double>> &EField_,
                              SphericalP<std::complex<double>> &HField_,
                              int projection_, CompoundIterator p,
                              int singleComponent_) {
  SphericalP<std::complex<double>> Efield = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Einc = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Hfield = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Hinc = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));

  Spherical<double> Rrel;

  std::complex<double> iZ =
      (consCmi / sqrt(geometry->bground.mu / geometry->bground.epsilon));

  int intInd = geometry->checkInner(R_);

  if (intInd < 0) // Outside a sphere
  {
    if (!flagSH) // this a fundamental frequency result - calculate the incoming
                 // field
    {
      // Incoming field
      optimet::AuxCoefficients aCoefInc(R_, waveK, 1, nMax);

      if (singleComponent_ == 1) // Only TE part
      {
        Einc =
            Einc + aCoefInc.M(static_cast<long>(p)) * excitation->dataIncAp[p];
        Hinc = Hinc +
               aCoefInc.M(static_cast<long>(p)) * excitation->dataIncBp[p] * iZ;
      }

      if (singleComponent_ == 2) // Only TM part
      {
        Einc =
            Einc + aCoefInc.N(static_cast<long>(p)) * excitation->dataIncBp[p];
        Hinc = Hinc +
               aCoefInc.N(static_cast<long>(p)) * excitation->dataIncAp[p] * iZ;
      }

      if (!singleComponent_) // Both TE and TM part
      {
        Einc = Einc +
               (aCoefInc.M(static_cast<long>(p)) * excitation->dataIncAp[p] +
                aCoefInc.N(static_cast<long>(p)) * excitation->dataIncBp[p]);
        Hinc = Hinc +
               (aCoefInc.N(static_cast<long>(p)) * excitation->dataIncAp[p] +
                aCoefInc.M(static_cast<long>(p)) * excitation->dataIncBp[p]) *
                   iZ;
      }
    } else // this a second harmonic frequency result - calculate the source
           // fields (save it in Einc for convenience)
    {
      // Source fields
      for (size_t j = 0; j < geometry->objects.size(); j++) {
        Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
        optimet::AuxCoefficients aCoef(Rrel, waveK, 0, nMax);

        Einc = Einc +
               aCoef.M(static_cast<long>(p)) *
                   geometry->objects[j].sourceCoef[static_cast<int>(p)] +
               aCoef.N(static_cast<long>(p)) *
                   geometry->objects[j]
                       .sourceCoef[static_cast<int>(p) + p.max(nMax)];
      }
    }

    // Scattered field
    for (size_t j = 0; j < geometry->objects.size(); j++) {
      SphericalP<std::complex<double>> Efield_local =
          SphericalP<std::complex<double>>(std::complex<double>(0.0, 0.0),
                                           std::complex<double>(0.0, 0.0),
                                           std::complex<double>(0.0, 0.0));

      Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
      optimet::AuxCoefficients aCoef(Rrel, waveK, 0, nMax);

      if (singleComponent_ == 1) // TE Part only
      {
        Efield = Efield +
                 aCoef.M(static_cast<long>(p)) *
                     scatter_coef[j * 2 * p.max(nMax) + p.compound];
        Hfield =
            Hfield +
            aCoef.M(static_cast<long>(p)) *
                scatter_coef[p.max(nMax) + j * 2 * p.max(nMax) + p.compound] *
                iZ;
      }

      if (singleComponent_ == 2) // TM Part only
      {
        Efield =
            Efield +
            aCoef.N(static_cast<long>(p)) *
                scatter_coef[p.max(nMax) + j * 2 * p.max(nMax) + p.compound];
        Hfield = Hfield +
                 aCoef.N(static_cast<long>(p)) *
                     scatter_coef[j * 2 * p.max(nMax) + p.compound] * iZ;
      }

      if (!singleComponent_) {
        Efield =
            Efield +
            aCoef.M(static_cast<long>(p)) *
                scatter_coef[j * 2 * p.max(nMax) + p.compound] +
            aCoef.N(static_cast<long>(p)) *
                scatter_coef[p.max(nMax) + j * 2 * p.max(nMax) + p.compound];
        Hfield =
            Hfield +
            (aCoef.N(static_cast<long>(p)) *
                 scatter_coef[j * 2 * p.max(nMax) + p.compound] +
             aCoef.M(static_cast<long>(p)) *
                 scatter_coef[p.max(nMax) + j * 2 * p.max(nMax) + p.compound]) *
                iZ;
      }
    }
  } else // Inside a sphere
  {
    Rrel = Tools::toPoint(R_, geometry->objects[intInd].vR);
    optimet::AuxCoefficients aCoef(
        Rrel, waveK * sqrt(geometry->objects[intInd].elmag.epsilon_r *
                           geometry->objects[intInd].elmag.mu_r),
        1, nMax);

    std::complex<double> iZ_object =
        (consCmi / sqrt(geometry->objects[intInd].elmag.mu /
                        geometry->objects[intInd].elmag.epsilon));

    if (singleComponent_ == 1) // TE Part Only
    {
      Efield = Efield +
               aCoef.M(static_cast<long>(p)) *
                   internal_coef[intInd * 2 * p.max(nMax) + p.compound];
      Hfield = Hfield +
               aCoef.M(static_cast<long>(p)) *
                   internal_coef[p.max(nMax) + intInd * 2 * p.max(nMax) +
                                 p.compound] *
                   iZ_object;
    }

    if (singleComponent_ == 2) // TM Part Only
    {
      Efield = Efield +
               aCoef.N(static_cast<long>(p)) *
                   internal_coef[p.max(nMax) + intInd * 2 * p.max(nMax) +
                                 p.compound];
      Hfield = Hfield +
               aCoef.N(static_cast<long>(p)) *
                   internal_coef[intInd * 2 * p.max(nMax) + p.compound] *
                   iZ_object;
    }

    if (!singleComponent_) // Both TE and TM
    {
      Efield =
          Efield + (aCoef.M(static_cast<long>(p)) *
                        internal_coef[intInd * 2 * p.max(nMax) + p.compound] +
                    aCoef.N(static_cast<long>(p)) *
                        internal_coef[p.max(nMax) + intInd * 2 * p.max(nMax) +
                                      p.compound]);
      Hfield = Hfield +
               (aCoef.N(static_cast<long>(p)) *
                    internal_coef[intInd * 2 * p.max(nMax) + p.compound] +
                aCoef.M(static_cast<long>(p)) *
                    internal_coef[p.max(nMax) + intInd * 2 * p.max(nMax) +
                                  p.compound]) *
                   iZ_object;
    }
  }

  if (projection_) {
    SphericalP<std::complex<double>> SphEField;
    SphericalP<std::complex<double>> SphHField;
    Rrel = Tools::toPoint(R_, geometry->objects[0].vR);

    SphEField = Tools::fromProjection(Rrel, Einc + Efield);
    SphHField = Tools::fromProjection(Rrel, Hinc + Hfield);

    EField_ = SphEField;
    HField_ = SphHField;
  } else {
    EField_ = Einc + Efield;
    HField_ = Hinc + Hfield;
  }
}

void Result::getEHFields(Spherical<double> R_,
                         SphericalP<std::complex<double>> &EField_,
                         SphericalP<std::complex<double>> &HField_,
                         int projection_) {
  SphericalP<std::complex<double>> Efield = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Einc = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Hfield = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Hinc = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));

  Spherical<double> Rrel;

  std::complex<double> iZ =
      (consCmi / sqrt(geometry->bground.mu / geometry->bground.epsilon));

  int pMax = Tools::iteratorMax(nMax);

  CompoundIterator p;

  // Check for inner point and set to 0
  int intInd = geometry->checkInner(R_);

  if (intInd < 0) // Outside a sphere
  {
    if (!flagSH) // this a fundamental frequency result - calculate the incoming
                 // field
    {
      // Incoming field
      optimet::AuxCoefficients aCoefInc(R_, waveK, 1, nMax);
      for (p = 0; p < p.max(nMax); p++) {
        Einc = Einc +
               (aCoefInc.M(static_cast<long>(p)) * excitation->dataIncAp[p] +
                aCoefInc.N(static_cast<long>(p)) * excitation->dataIncBp[p]);
        Hinc = Hinc +
               (aCoefInc.N(static_cast<long>(p)) * excitation->dataIncAp[p] +
                aCoefInc.M(static_cast<long>(p)) * excitation->dataIncBp[p]) *
                   iZ;
      }
    } else // this a second harmonic frequency result - calculate the source
           // fields (save it in Einc for convenience)
    {
      // Source fields
      for (size_t j = 0; j < geometry->objects.size(); j++) {
        Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
        optimet::AuxCoefficients aCoef(Rrel, waveK, 0, nMax);

        for (p = 0; p < pMax; p++) {
          Einc =
              Einc +
              aCoef.M(static_cast<long>(p)) *
                  geometry->objects[j].sourceCoef[static_cast<int>(p)] +
              aCoef.N(static_cast<long>(p)) *
                  geometry->objects[j].sourceCoef[static_cast<int>(p) + pMax];
        }
      }
    }

    // Scattered field
    for (size_t j = 0; j < geometry->objects.size(); j++) {
      SphericalP<std::complex<double>> Efield_local =
          SphericalP<std::complex<double>>(std::complex<double>(0.0, 0.0),
                                           std::complex<double>(0.0, 0.0),
                                           std::complex<double>(0.0, 0.0));

      Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
      optimet::AuxCoefficients aCoef(Rrel, waveK, 0, nMax);

      for (p = 0; p < p.max(nMax); p++) {
        Efield = Efield +
                 aCoef.M(static_cast<long>(p)) *
                     scatter_coef[j * 2 * pMax + p.compound] +
                 aCoef.N(static_cast<long>(p)) *
                     scatter_coef[pMax + j * 2 * pMax + p.compound];
        Hfield = Hfield +
                 (aCoef.N(static_cast<long>(p)) *
                      scatter_coef[j * 2 * pMax + p.compound] +
                  aCoef.M(static_cast<long>(p)) *
                      scatter_coef[pMax + j * 2 * pMax + p.compound]) *
                     iZ;
      }
    }
  } else // Inside a sphere
  {
    Rrel = Tools::toPoint(R_, geometry->objects[intInd].vR);
    optimet::AuxCoefficients aCoef(
        Rrel, waveK * sqrt(geometry->objects[intInd].elmag.epsilon_r *
                           geometry->objects[intInd].elmag.mu_r),
        1, nMax);

    std::complex<double> iZ_object =
        (consCmi / sqrt(geometry->objects[intInd].elmag.mu /
                        geometry->objects[intInd].elmag.epsilon));

    for (p = 0; p < p.max(nMax); p++) {
      Efield =
          Efield + (aCoef.M(static_cast<long>(p)) *
                        internal_coef[intInd * 2 * pMax + p.compound] +
                    aCoef.N(static_cast<long>(p)) *
                        internal_coef[pMax + intInd * 2 * pMax + p.compound]);
      Hfield = Hfield +
               (aCoef.N(static_cast<long>(p)) *
                    internal_coef[intInd * 2 * pMax + p.compound] +
                aCoef.M(static_cast<long>(p)) *
                    internal_coef[pMax + intInd * 2 * pMax + p.compound]) *
                   iZ_object;
    }
  }

  if (projection_) {
    SphericalP<std::complex<double>> SphEField;
    SphericalP<std::complex<double>> SphHField;
    Rrel = Tools::toPoint(R_, geometry->objects[0].vR);

    SphEField = Tools::fromProjection(Rrel, Einc + Efield);
    SphHField = Tools::fromProjection(Rrel, Hinc + Hfield);

    EField_ = SphEField;
    HField_ = SphHField;
  } else {
    EField_ = Einc + Efield;
    HField_ = Hinc + Hfield;
  }
}

SphericalP<std::complex<double>> Result::getEFieldC(Spherical<double> R_,
                                                    int projection) {
  /* TEST FUNCTION. NOT USED IN PRODUCTION CODE! */

  SphericalP<std::complex<double>> Efield = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Einc = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  Spherical<double> Rrel;

  int pMax = Tools::iteratorMax(nMax);

  CompoundIterator p;

  // Check for inner point and set to 0
  int intInd = geometry->checkInner(R_);

  if (intInd < 0) // Outside a sphere
  {
    if (!flagSH) // this a fundamental frequency result - calculate the incoming
                 // field
    {
      // Incoming field
      optimet::AuxCoefficients aCoefInc(R_, waveK, 1, nMax);
      for (p = 0; p < p.max(nMax); p++) {
        Einc = Einc +
               (aCoefInc.M(static_cast<long>(p)) * excitation->dataIncAp[p] +
                aCoefInc.N(static_cast<long>(p)) * excitation->dataIncBp[p]);
      }

      if (projection) // Spherical projection
      {
        Einc = Tools::fromProjection(R_, Einc);
      }
    } else // this a second harmonic frequency result - calculate the source
           // fields (save it in Einc for convenience)
    {
      // Source fields
      for (size_t j = 0; j < geometry->objects.size(); j++) {
        Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
        optimet::AuxCoefficients aCoef(Rrel, waveK, 0, nMax);

        for (p = 0; p < pMax; p++) {
          Einc =
              Einc +
              aCoef.M(static_cast<long>(p)) *
                  geometry->objects[j].sourceCoef[static_cast<int>(p)] +
              aCoef.N(static_cast<long>(p)) *
                  geometry->objects[j].sourceCoef[static_cast<int>(p) + pMax];
        }
      }

      if (projection) // Spherical projection
      {
        Einc = Tools::fromProjection(R_, Einc);
      }
    }

    optimet::AuxCoefficients aCoef(R_, waveK, 0, nMax);

    // Scattered field
    for (p = 0; p < p.max(nMax); p++) {
      Efield =
          Efield +
          (aCoef.M(static_cast<long>(p)) * c_scatter_coef[p.compound] +
           aCoef.N(static_cast<long>(p)) * c_scatter_coef[pMax + p.compound]);
    }

    if (projection) // Spherical projection
    {
      Efield = Tools::fromProjection(R_, Efield);
    }
  } else // Inside a sphere
  {
    Rrel = Tools::toPoint(R_, geometry->objects[intInd].vR);
    optimet::AuxCoefficients aCoef(Rrel, waveK, 1, nMax);

    for (p = 0; p < p.max(nMax); p++) {
      Efield =
          Efield + (aCoef.M(static_cast<long>(p)) *
                        internal_coef[intInd * 2 * pMax + p.compound] +
                    aCoef.N(static_cast<long>(p)) *
                        internal_coef[pMax + intInd * 2 * pMax + p.compound]);
    }

    if (projection) // Spherical projection
    {
      Efield = Tools::fromProjection(R_, Efield);
    }
  }

  return Einc + Efield;
}

double Result::getExtinctionCrossSection() {
  CompoundIterator p;
  int pMax = Tools::iteratorMax(nMax);

  double Cext(0.);
  std::complex<double> *Q_local = new std::complex<double>[2 * pMax];

  for (size_t j = 0; j < geometry->objects.size(); j++) {
    excitation->getIncLocal(geometry->objects[j].vR, Q_local, nMax);
    for (p = 0; p < pMax; p++) {
      Cext += std::real(std::conj(Q_local[p]) *
                            scatter_coef[j * 2 * pMax + p.compound] +
                        std::conj(Q_local[p.compound + pMax]) *
                            scatter_coef[pMax + j * 2 * pMax + p.compound]);
    }
  }

  delete[] Q_local;
  return (-1. / (std::real(waveK) * std::real(waveK))) * Cext;
}

double Result::getAbsorptionCrossSection() {
  CompoundIterator p;
  int pMax = Tools::iteratorMax(nMax);

  double Cabs(0.);
  double temp1(0.), temp2(0.);
  double *Cabs_aux = new double[2 * pMax];

  auto const omega = excitation->omega();
  for (size_t j = 0; j < geometry->objects.size(); j++) {

    geometry->getCabsAux(omega, j, nMax, Cabs_aux);

    for (p = 0; p < pMax; p++) {
      temp1 = abs(scatter_coef[j * 2 * pMax + p.compound]);
      temp1 *= temp1;
      temp2 = abs(scatter_coef[pMax + j * 2 * pMax + p.compound]);
      temp2 *= temp2;
      Cabs +=
          temp1 * Cabs_aux[p.compound] + temp2 * Cabs_aux[pMax + p.compound];
    }
  }
  delete[] Cabs_aux;
  return (1 / (std::real(waveK) * std::real(waveK))) * Cabs;
}

int Result::setFields(OutputGrid &oEGrid_, OutputGrid &oHGrid_, int projection_) {
  Spherical<double> Rloc;

  // centerScattering();

  // Calculate the fields
  while (!oEGrid_.gridDone) {
    Rloc = oEGrid_.getPoint();
    oHGrid_.getPoint();
    std::cout << "Calculating fields for point " << oEGrid_.iterator + 1
              << " out of " << oEGrid_.gridPoints << std::endl;

    SphericalP<std::complex<double>> EField;
    SphericalP<std::complex<double>> HField;

    getEHFields(Rloc, EField, HField, projection_);

    oHGrid_.pushDataNext(HField);
    oEGrid_.pushDataNext(EField);
  }

  return 0;
}

int Result::setFieldsModal(OutputGrid &oEGrid_, OutputGrid &oHGrid_, int projection_,
                           CompoundIterator p_, int singleComponent_) {
  Spherical<double> Rloc;

  // Calculate the fields
  while (!oEGrid_.gridDone) {
    Rloc = oEGrid_.getPoint();
    oHGrid_.getPoint();
    std::cout << "Calculating fields for point " << oEGrid_.iterator + 1
              << " out of " << oEGrid_.gridPoints << std::endl;

    SphericalP<std::complex<double>> EField;
    SphericalP<std::complex<double>> HField;

    getEHFieldsModal(Rloc, EField, HField, projection_, p_, singleComponent_);

    oHGrid_.pushDataNext(HField);
    oEGrid_.pushDataNext(EField);
  }
  return 0;
}

void Result::centerScattering() {
  CompoundIterator p, q;

  int pMax = p.max(nMax);
  int qMax = q.max(nMax);

  for (p = 0; p < 2 * p.max(nMax); p++) {
    c_scatter_coef[p] = std::complex<double>(0.0, 0.0);
  }

  std::complex<double> **T_AB = new std::complex<double> *[2 * (p.max(nMax))];
  std::complex<double> *scatter_aux = new std::complex<double>[2 * p.max(nMax)];
  std::complex<double> *scatter_fin = new std::complex<double>[2 * p.max(nMax)];

  for (p = 0; p < (int)(2 * p.max(nMax)); p++) {
    T_AB[p] = new std::complex<double>[2 * p.max(nMax)];
  }

  for (size_t i = 0; i < geometry->objects.size(); i++) {
    Spherical<double> Rrel = Tools::toPoint(Spherical<double>(0.0, 0.0, 0.0),
                                            geometry->objects[i].vR);

    optimet::Coupling const coupling(Rrel, excitation->waveK, nMax);

    for (p = 0; p < pMax; p++)
      for (q = 0; q < qMax; q++) {
        T_AB[p][q] = coupling.diagonal(p, q);
        T_AB[p + pMax][q + qMax] = coupling.diagonal(p, q);
        T_AB[p + pMax][q] = coupling.offdiagonal(p, q);
        T_AB[p][q + qMax] = coupling.offdiagonal(p, q);
      }

    for (p = 0; p < 2 * p.max(nMax); p++) {
      scatter_aux[p] +=
          scatter_coef[p.compound + i * 2 * pMax * geometry->objects.size()];
    }

    Algebra::multiplyVectorMatrix(T_AB, 2 * pMax, 2 * pMax, scatter_aux,
                                  scatter_fin, consC1, consC0);

    for (p = 0; p < 2 * pMax; p++) {
      c_scatter_coef[p] += scatter_fin[p];
    }
  }

  for (p = 0; p < (int)(2 * p.max(nMax)); p++) {
    delete[] T_AB[p];
  }

  delete[] T_AB;
  delete[] scatter_aux;
  delete[] scatter_fin;
}

CompoundIterator Result::getDominant() {
  CompoundIterator p, q;

  q = 0;

  std::complex<double> TEMax = scatter_coef[0];
  std::complex<double> TMMax = scatter_coef[p.max(nMax)];

  for (p = 0; p < p.max(nMax); p++) {
    if ((abs(scatter_coef[p]) > abs(TEMax)) ||
        (abs(scatter_coef[p.max(nMax) + p.compound]) > abs(TMMax))) {
      q = p;

      TEMax = scatter_coef[p];
      TMMax = scatter_coef[p.max(nMax) + p.compound];
    }
  }

  return q;
}

void Result::getEHFieldsContCheck(Spherical<double> R_,
                                  SphericalP<std::complex<double>> &EField_,
                                  SphericalP<std::complex<double>> &HField_,
                                  int projection_, int inside_) {
  SphericalP<std::complex<double>> Efield = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Einc = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Hfield = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Hinc = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));

  Spherical<double> Rrel;

  std::complex<double> iZ =
      (consCmi / sqrt(geometry->bground.mu / geometry->bground.epsilon));

  int pMax = Tools::iteratorMax(nMax);

  CompoundIterator p;

  // Check for inner point and set to 0
  int intInd = geometry->checkInner(R_);
  intInd = inside_;

  if (intInd < 0) // Outside a sphere
  {
    if (!flagSH) // this a fundamental frequency result - calculate the incoming
                 // field
    {
      // Incoming field
      optimet::AuxCoefficients aCoefInc(R_, waveK, 1, nMax);
      for (p = 0; p < p.max(nMax); p++) {
        Einc = Einc +
               (aCoefInc.M(static_cast<long>(p)) * excitation->dataIncAp[p] +
                aCoefInc.N(static_cast<long>(p)) * excitation->dataIncBp[p]);
        Hinc = Hinc +
               (aCoefInc.N(static_cast<long>(p)) * excitation->dataIncAp[p] +
                aCoefInc.M(static_cast<long>(p)) * excitation->dataIncBp[p]) *
                   iZ;
      }
    } else // this a second harmonic frequency result - calculate the source
           // fields (save it in Einc for convenience)
    {
      // Source fields
      for (size_t j = 0; j < geometry->objects.size(); j++) {
        Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
        optimet::AuxCoefficients aCoef(Rrel, waveK, 0, nMax);

        for (p = 0; p < pMax; p++) {
          Einc =
              Einc +
              aCoef.M(static_cast<long>(p)) *
                  geometry->objects[j].sourceCoef[static_cast<int>(p)] +
              aCoef.N(static_cast<long>(p)) *
                  geometry->objects[j].sourceCoef[static_cast<int>(p) + pMax];
        }
      }
    }

    // Scattered field
    for (size_t j = 0; j < geometry->objects.size(); j++) {
      SphericalP<std::complex<double>> Efield_local =
          SphericalP<std::complex<double>>(std::complex<double>(0.0, 0.0),
                                           std::complex<double>(0.0, 0.0),
                                           std::complex<double>(0.0, 0.0));

      Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
      optimet::AuxCoefficients aCoef(Rrel, waveK, 0, nMax);

      for (p = 0; p < p.max(nMax); p++) {
        Efield = Efield +
                 aCoef.M(static_cast<long>(p)) *
                     scatter_coef[j * 2 * pMax + p.compound] +
                 aCoef.N(static_cast<long>(p)) *
                     scatter_coef[pMax + j * 2 * pMax + p.compound];
        Hfield = Hfield +
                 (aCoef.N(static_cast<long>(p)) *
                      scatter_coef[j * 2 * pMax + p.compound] +
                  aCoef.M(static_cast<long>(p)) *
                      scatter_coef[pMax + j * 2 * pMax + p.compound]) *
                     iZ;
      }
    }
  } else // Inside a sphere
  {
    Rrel = Tools::toPoint(R_, geometry->objects[intInd].vR);
    optimet::AuxCoefficients aCoef(
        Rrel, waveK * sqrt(geometry->objects[intInd].elmag.epsilon_r *
                           geometry->objects[intInd].elmag.mu_r),
        1, nMax);

    std::complex<double> iZ_object =
        (consCmi / sqrt(geometry->objects[intInd].elmag.mu /
                        geometry->objects[intInd].elmag.epsilon));

    for (p = 0; p < p.max(nMax); p++) {
      Efield =
          Efield + (aCoef.M(static_cast<long>(p)) *
                        internal_coef[intInd * 2 * pMax + p.compound] +
                    aCoef.N(static_cast<long>(p)) *
                        internal_coef[pMax + intInd * 2 * pMax + p.compound]);
      Hfield = Hfield +
               (aCoef.N(static_cast<long>(p)) *
                    internal_coef[intInd * 2 * pMax + p.compound] +
                aCoef.M(static_cast<long>(p)) *
                    internal_coef[pMax + intInd * 2 * pMax + p.compound]) *
                   iZ_object;
    }
  }

  if (projection_) {
    SphericalP<std::complex<double>> SphEField;
    SphericalP<std::complex<double>> SphHField;
    Rrel = Tools::toPoint(R_, geometry->objects[0].vR);

    SphEField = Tools::fromProjection(Rrel, Einc + Efield);
    SphHField = Tools::fromProjection(Rrel, Hinc + Hfield);

    EField_ = SphEField;
    HField_ = SphHField;
  }

  else { // AJ - no spherical projection
    EField_ = Einc + Efield;
    HField_ = Hinc + Hfield;
  }
}

void Result::writeContinuityCheck(int objectIndex_) {
  SphericalP<std::complex<double>> AnEField_in, AnEField_out;
  SphericalP<std::complex<double>> AnHField_in, AnHField_out;
  Spherical<double> APoint(0.0, 0.0, 0.0);
  int projection = 1; // Spherical projection - True - projection is internally
                      // set to be evaluated w.r.t. object[0]
  int outside = -1;   // Forces result to be outside an object
  int inside = 0;     // Forces result to be inside an object
  double radius = geometry->objects[objectIndex_].radius;
  ;
  std::complex<double> eps_r = geometry->objects[objectIndex_].elmag.epsilon_r;
  std::complex<double> mu_r = geometry->objects[objectIndex_].elmag.mu_r;

  // the-phi - 2D plot
  // ----------------------------------------------------------------------------------------------------------------
  std::ofstream E1_err_mag("E1_err_mag");
  std::ofstream E2_err_mag("E2_err_mag");
  std::ofstream E3_err_mag("E3_err_mag");
  std::ofstream H1_err_mag("H1_err_mag");
  std::ofstream H2_err_mag("H2_err_mag");
  std::ofstream H3_err_mag("H3_err_mag");
  int max_ii =
      180; // theta observation range (increment by 1 degree) - [1, max_ii-1]
  int max_jj =
      180; // phi   observation range (increment by 1 degree) - [1, max_jj-1]
  for (int ii = 1; ii <= max_ii - 1; ii++) {
    for (int jj = 1; jj <= max_jj - 1; jj++) {

      APoint = Spherical<double>(radius, consPi * (double(ii) / 180.),
                                 consPi * (double(jj) / 180.));
      getEHFieldsContCheck(APoint, AnEField_out, AnHField_out, projection,
                           outside);
      getEHFieldsContCheck(APoint, AnEField_in, AnHField_in, projection,
                           inside);
      // op -------------------------------------------------------------------
      std::cout << "Continuity check : computed "
                << jj + ((ii - 1) * (max_ii - 1)) << " out of a total of "
                << (max_ii - 1) * (max_jj - 1) << std::endl;
      // EF
      E1_err_mag << (abs(AnEField_out.rrr) - abs(AnEField_in.rrr * eps_r)) /
                        abs(AnEField_in.rrr * eps_r)
                 << " ";
      E2_err_mag << (abs(AnEField_out.the) - abs(AnEField_in.the)) /
                        abs(AnEField_in.the)
                 << " ";
      E3_err_mag << (abs(AnEField_out.phi) - abs(AnEField_in.phi)) /
                        abs(AnEField_in.phi)
                 << " ";
      // EF
      H1_err_mag << (abs(AnHField_out.rrr) - abs(AnHField_in.rrr * mu_r)) /
                        abs(AnHField_in.rrr * mu_r)
                 << " ";
      H2_err_mag << (abs(AnHField_out.the) - abs(AnHField_in.the)) /
                        abs(AnHField_in.the)
                 << " ";
      H3_err_mag << (abs(AnHField_out.phi) - abs(AnHField_in.phi)) /
                        abs(AnHField_in.phi)
                 << " ";
    }
    E1_err_mag << std::endl;
    E2_err_mag << std::endl;
    E3_err_mag << std::endl;
    H1_err_mag << std::endl;
    H2_err_mag << std::endl;
    H3_err_mag << std::endl;
  }
  E1_err_mag.flush();
  E2_err_mag.flush();
  E3_err_mag.flush();
  H1_err_mag.flush();
  H2_err_mag.flush();
  H3_err_mag.flush();
  E1_err_mag.close();
  E2_err_mag.close();
  E3_err_mag.close();
  H1_err_mag.close();
  H2_err_mag.close();
  H3_err_mag.close();
}
}
