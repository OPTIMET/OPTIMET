// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.


#include "Algebra.h"
#include "AuxCoefficients.h"
#include "Coupling.h"
#include "Result.h"
#include "Tools.h"
#include "constants.h"

#include <complex>
#include <cstdlib>
#include <fstream>
#include <iostream>

namespace optimet {
Result::Result(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> excitation_)
    : flagSH(true), result_FF(nullptr) {
     
  init(geometry_, excitation_);
}

Result::Result(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> excitation_,
               Result *result_FF_) {
  init(geometry_, excitation_, result_FF_);
}

void Result::init(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> excitation_) {
  geometry = geometry_;
  excitation = excitation_;
  waveK = excitation->waveK;
  flagSH = true;
  result_FF = NULL;
  result_SH = NULL;
  nMax = geometry_->nMax();
  nMaxS = geometry_->nMaxS();

  CLGcoeff.resize(9); 
  scatter_coef.resize(2 * Tools::iteratorMax(nMax) * geometry->objects.size());
  internal_coef.resize(2 * Tools::iteratorMax(nMax) * geometry->objects.size());
  scatter_coef_SH.resize(4 * Tools::iteratorMax(nMaxS) * geometry->objects.size());
  internal_coef_SH.resize(4 * Tools::iteratorMax(nMaxS) * geometry->objects.size());
  c_scatter_coef.resize(2 * Tools::iteratorMax(nMax));
}

void Result::update(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> excitation_) {
  geometry = geometry_;
  excitation = excitation_;
  waveK = excitation->waveK;
}

void Result::init(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> excitation_,
                  Result *result_FF_) {
  init(geometry_, excitation_);
  flagSH = true;
  result_FF = result_FF_;
}

void Result::getEHFieldsModal(Spherical<double> R_, SphericalP<std::complex<double>> &EField_,
                              SphericalP<std::complex<double>> &HField_, bool projection_,
                              CompoundIterator p, int singleComponent_) {
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

  std::complex<double> iZ = (consCmi / sqrt(geometry->bground.mu / geometry->bground.epsilon));

  int intInd = geometry->checkInner(R_);

  if(intInd < 0) // Outside a sphere
  {
    if(!flagSH) // this a fundamental frequency result - calculate the incoming
                // field
    {
      // Incoming field
      optimet::AuxCoefficients aCoefInc(R_, waveK, 1, nMax);

      if(singleComponent_ == 1) // Only TE part
      {
        Einc = Einc + aCoefInc.M(static_cast<long>(p)) * excitation->dataIncAp[p];
        Hinc = Hinc + aCoefInc.M(static_cast<long>(p)) * excitation->dataIncBp[p] * iZ;
      }

      if(singleComponent_ == 2) // Only TM part
      {
        Einc = Einc + aCoefInc.N(static_cast<long>(p)) * excitation->dataIncBp[p];
        Hinc = Hinc + aCoefInc.N(static_cast<long>(p)) * excitation->dataIncAp[p] * iZ;
      }

      if(!singleComponent_) // Both TE and TM part
      { 
        Einc = Einc + (aCoefInc.M(static_cast<long>(p)) * excitation->dataIncAp[p] +
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
      for(size_t j = 0; j < geometry->objects.size(); j++) {
        Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
        optimet::AuxCoefficients aCoef(Rrel, waveK, 0, nMax);
		 
        Einc =
            Einc +
            aCoef.M(static_cast<long>(p)) * geometry->objects[j].sourceCoef[static_cast<int>(p)] +
            aCoef.N(static_cast<long>(p)) *
                geometry->objects[j].sourceCoef[static_cast<int>(p) + p.max(nMax)];
      }
    }

    // Scattered field
    for(size_t j = 0; j < geometry->objects.size(); j++) {
      SphericalP<std::complex<double>> Efield_local = SphericalP<std::complex<double>>(
          std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
          std::complex<double>(0.0, 0.0));
         

      Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
      optimet::AuxCoefficients aCoef(Rrel, waveK, 0, nMax);

      if(singleComponent_ == 1) // TE Part only
      {
        Efield =
            Efield + aCoef.M(static_cast<long>(p)) * scatter_coef[j * 2 * p.max(nMax) + p.compound];
        Hfield = Hfield +
                 aCoef.M(static_cast<long>(p)) *
                     scatter_coef[p.max(nMax) + j * 2 * p.max(nMax) + p.compound] * iZ;
      }

      if(singleComponent_ == 2) // TM Part only
      {
        Efield = Efield +
                 aCoef.N(static_cast<long>(p)) *
                     scatter_coef[p.max(nMax) + j * 2 * p.max(nMax) + p.compound];
        Hfield =
            Hfield +
            aCoef.N(static_cast<long>(p)) * scatter_coef[j * 2 * p.max(nMax) + p.compound] * iZ;
            
      }

      if(!singleComponent_) {
        Efield = Efield +
                 aCoef.M(static_cast<long>(p)) * scatter_coef[j * 2 * p.max(nMax) + p.compound] +
                 aCoef.N(static_cast<long>(p)) *
                     scatter_coef[p.max(nMax) + j * 2 * p.max(nMax) + p.compound];
        Hfield = Hfield +
                 (aCoef.N(static_cast<long>(p)) * scatter_coef[j * 2 * p.max(nMax) + p.compound] +
                  aCoef.M(static_cast<long>(p)) *
                      scatter_coef[p.max(nMax) + j * 2 * p.max(nMax) + p.compound]) *
                     iZ;
                                 
      }
    }
  } else // Inside a sphere
  {
    Rrel = Tools::toPoint(R_, geometry->objects[intInd].vR);
    optimet::AuxCoefficients aCoef(Rrel, waveK * sqrt(geometry->objects[intInd].elmag.epsilon_r *
                                                      geometry->objects[intInd].elmag.mu_r),
                                   1, nMax);

    std::complex<double> iZ_object = (consCmi / sqrt(geometry->objects[intInd].elmag.mu /
                                                     geometry->objects[intInd].elmag.epsilon));

    if(singleComponent_ == 1) // TE Part Only
    {
      Efield = Efield +
               aCoef.M(static_cast<long>(p)) * internal_coef[intInd * 2 * p.max(nMax) + p.compound];
      Hfield = Hfield +
               aCoef.M(static_cast<long>(p)) *
                   internal_coef[p.max(nMax) + intInd * 2 * p.max(nMax) + p.compound] * iZ_object;
    }

    if(singleComponent_ == 2) // TM Part Only
    {
      Efield = Efield +
               aCoef.N(static_cast<long>(p)) *
                   internal_coef[p.max(nMax) + intInd * 2 * p.max(nMax) + p.compound];
      Hfield = Hfield +
               aCoef.N(static_cast<long>(p)) *
                   internal_coef[intInd * 2 * p.max(nMax) + p.compound] * iZ_object;
    }

    if(!singleComponent_) // Both TE and TM
    {
      Efield = Efield + (aCoef.M(static_cast<long>(p)) *
                             internal_coef[intInd * 2 * p.max(nMax) + p.compound] +
                         aCoef.N(static_cast<long>(p)) *
                             internal_coef[p.max(nMax) + intInd * 2 * p.max(nMax) + p.compound]);
      Hfield =
          Hfield +
          (aCoef.N(static_cast<long>(p)) * internal_coef[intInd * 2 * p.max(nMax) + p.compound] +
           aCoef.M(static_cast<long>(p)) *
               internal_coef[p.max(nMax) + intInd * 2 * p.max(nMax) + p.compound]) *
              iZ_object;
    }
  }

  if(projection_) {
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

void Result::getEHFields(Spherical<double> R_, SphericalP<std::complex<double>> &EField_FF,
                         SphericalP<std::complex<double>> &HField_FF, SphericalP<std::complex<double>> &EField_SH,
                         SphericalP<std::complex<double>> &HField_SH, bool projection_, std::complex<double> *coeffXmn, std::complex<double> *coeffXpl) const {
                                                 
  SphericalP<std::complex<double>> Efield_FF = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Einc_FF = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Hfield_FF = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Hinc_FF = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
                                                       
   
      
  SphericalP<std::complex<double>> Efield_SH = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Egamma_SH = SphericalP<std::complex<double>>( // particular solution of diff equations
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
  SphericalP<std::complex<double>> Hfield_SH = SphericalP<std::complex<double>>(
      std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
 
  int NO = geometry->objects.size(); // how many particles there is
  
   std::complex<double> amnSH, bmnSH, cmnSH, dmnSH;
  
  auto const omega = excitation->omega();
  
  const std::complex<double> waveK_0 = (omega) * std::sqrt(consEpsilon0 * consMu0);
  
  Spherical<double> Rrel;

  std::complex<double> iZ = (consCmi / sqrt(geometry->bground.mu / geometry->bground.epsilon));
  
  int pMax = Tools::iteratorMax(nMax);
  int pMaxS = Tools::iteratorMax(nMaxS);
  
  CompoundIterator p;

  // Check for inner point and set to 0
  int intInd = geometry->checkInner(R_);
  
  
  if(intInd < 0) // Outside a sphere
  {
    // this a fundamental frequency result - calculate the incoming
                // field
    
      // Incoming field
      optimet::AuxCoefficients aCoefInc(R_, waveK, 1, nMax); //regular spherical functions
      
      for(p = 0; p < p.max(nMax); p++) {
      
        Einc_FF = Einc_FF + (aCoefInc.M(static_cast<long>(p)) * excitation->dataIncAp[p] +
                       aCoefInc.N(static_cast<long>(p)) * excitation->dataIncBp[p]); 
        Hinc_FF = Hinc_FF +
               (aCoefInc.N(static_cast<long>(p)) * excitation->dataIncAp[p] +
                aCoefInc.M(static_cast<long>(p)) * excitation->dataIncBp[p]) *
                   iZ;           
      }
    

    // Scattered field
       // this a fundamental frequency result 
     
    for(size_t j = 0; j < geometry->objects.size(); j++) {
      

      Rrel = Tools::toPoint(R_, geometry->objects[j].vR);

 
      optimet::AuxCoefficients aCoefFF(Rrel, waveK, 0, nMax); //radiative spherical vector functions

       
      for(p = 0; p < p.max(nMax); p++) {
     
        
        Efield_FF = Efield_FF + aCoefFF.M(static_cast<long>(p)) * scatter_coef[j * 2 * pMax + p.compound] +
                 aCoefFF.N(static_cast<long>(p)) * scatter_coef[pMax + j * 2 * pMax + p.compound];
                 
        Hfield_FF = Hfield_FF +
                 (aCoefFF.N(static_cast<long>(p)) * scatter_coef[j * 2 * pMax + p.compound] +
                  aCoefFF.M(static_cast<long>(p)) * scatter_coef[pMax + j * 2 * pMax + p.compound]) *
                     (iZ);
      }
            
    }
    
  
      // this a second harmonic frequency result 
  
  for(size_t j = 0; j < geometry->objects.size(); j++) {


      Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
      
      optimet::AuxCoefficients aCoefSH(Rrel, std::complex<double>(2.0, 0.0) * waveK, 0, nMaxS); //radiative spherical

      for(p = 0; p < p.max(nMaxS); p++) {
      
      bmnSH = scatter_coef_SH[j * 4 * pMaxS + p.compound] + scatter_coef_SH[2 * pMaxS + j * 4 * pMaxS + p.compound];
      
      amnSH = scatter_coef_SH[j * 4 * pMaxS + pMaxS + p.compound] + scatter_coef_SH[j * 4 * pMaxS + 3 * pMaxS +  p.compound];
     
        Efield_SH = Efield_SH + (aCoefSH.M(static_cast<long>(p)) * (bmnSH) +
                                 aCoefSH.N(static_cast<long>(p)) * (amnSH)) * waveK_0; 
        
        Hfield_SH = Hfield_SH +
                 (aCoefSH.N(static_cast<long>(p)) * (bmnSH) +
                  aCoefSH.M(static_cast<long>(p)) * (amnSH)) * (iZ * waveK_0);
                      
      }
      
    }    
    
   } // if
 
 
  else // Inside a sphere
 
    {
    
    Rrel = Tools::toPoint(R_, geometry->objects[intInd].vR);
    
    optimet::AuxCoefficients aCoefFF(Rrel, waveK_0 * sqrt(geometry->objects[intInd].elmag.epsilon_r *
                                                      geometry->objects[intInd].elmag.mu_r), 1, nMax); //regular vector spherical functions
                                                      
                                   
    std::complex<double> iZ_object = (consCmi / sqrt(geometry->objects[intInd].elmag.mu /
                                                     geometry->objects[intInd].elmag.epsilon));
                                                                                                     
  
    for(p = 0; p < p.max(nMax); p++) {
      Efield_FF =
          Efield_FF +
          (aCoefFF.M(static_cast<long>(p)) * internal_coef[intInd * 2 * pMax + p.compound] +
           aCoefFF.N(static_cast<long>(p)) * internal_coef[pMax + intInd * 2 * pMax + p.compound]);
           
                
      Hfield_FF =
          Hfield_FF +
          (aCoefFF.N(static_cast<long>(p)) * internal_coef[intInd * 2 * pMax + p.compound] +
           aCoefFF.M(static_cast<long>(p)) * internal_coef[pMax + intInd * 2 * pMax + p.compound]) *
              iZ_object; 
                       
              
         
    }
    
    
    // this a second harmonic frequency result 
    optimet::AuxCoefficients aCoefSH(Rrel, std::complex<double> (2.0 , 0.0) * waveK_0 * sqrt(geometry->objects[intInd].elmag.epsilon_r_SH * geometry->objects[intInd].elmag.mu_r_SH),
                                   1, nMaxS); // regular vector spherical functions
                                   
                                   
    std::complex<double> iZ_object_SH = (consCmi / sqrt(geometry->objects[intInd].elmag.mu_SH /
                                                     geometry->objects[intInd].elmag.epsilon_SH));
                                                                                                                                                       

    for(p = 0; p < p.max(nMaxS); p++) {
    
    cmnSH = internal_coef_SH[intInd * 4 * pMaxS + p.compound] + internal_coef_SH[2 * pMaxS + intInd * 4 * pMaxS + p.compound];

    dmnSH = internal_coef_SH[pMaxS + intInd * 4 * pMaxS + p.compound] + 
    internal_coef_SH[3 * pMaxS + intInd * 4 * pMaxS + p.compound];
    
     Efield_SH =
          Efield_SH +
          (aCoefSH.M(static_cast<long>(p)) * (cmnSH) +
           aCoefSH.N(static_cast<long>(p)) * (dmnSH)) * (waveK_0); //predznak postavljen kao u papiru
           
      Egamma_SH = Egamma_SH +
             (aCoefSH.Xm(static_cast<long>(p)) * coeffXmn[p] +
            aCoefSH.Xp(static_cast<long>(p)) * coeffXpl[p]);
        
      Hfield_SH =
          Hfield_SH +
          (aCoefSH.N(static_cast<long>(p)) * (cmnSH) +
           aCoefSH.M(static_cast<long>(p)) * (dmnSH)) * iZ_object_SH * waveK_0;  // predznak postavljen kao u papiru
              

   } 
     
}

  if(projection_) {
    SphericalP<std::complex<double>> SphEField;
    SphericalP<std::complex<double>> SphHField;
    Rrel = Tools::toPoint(R_, geometry->objects[0].vR);

    SphEField = Tools::fromProjection(Rrel, Einc_FF + Efield_FF);
    SphHField = Tools::fromProjection(Rrel, Hinc_FF + Hfield_FF);
    
    EField_FF = SphEField;
    HField_FF = SphHField;
      
    
  } else {
    EField_FF = Einc_FF + Efield_FF;
    HField_FF = Hinc_FF + Hfield_FF;
   
    EField_SH = Efield_SH + Egamma_SH;
    HField_SH = Hfield_SH;
 
 
  }
  
}



double Result::getExtinctionCrossSection(int gran1, int gran2) {
  CompoundIterator p;
  int pMax = Tools::iteratorMax(nMax);

  double Cext(0.);
  std::complex<double> *Q_local = new std::complex<double>[2 * pMax];

  for(int j = gran1; j < gran2; j++) {

    excitation->getIncLocal(geometry->objects[j].vR, Q_local, nMax);
    for(p = 0; p < pMax; p++) {
      Cext += std::real(std::conj(Q_local[p]) * scatter_coef[j * 2 * pMax + p.compound] +
                        std::conj(Q_local[p.compound + pMax]) *
                            scatter_coef[pMax + j * 2 * pMax + p.compound]);
    }
  }
  
  

  delete[] Q_local;
  return (-1. / (std::real(waveK) * std::real(waveK))) * Cext;
}



double Result::getAbsorptionCrossSection(int gran1, int gran2) {
 
  CompoundIterator p;
  int pMax = Tools::iteratorMax(nMax);

  double Cabs(0.);
  double temp1(0.), temp2(0.);
  double *Cabs_aux = new double[2 * pMax];
  auto const omega = excitation->omega();
  for(int j = gran1; j < gran2; j++) {

    geometry->getCabsAux(omega, j, nMax, Cabs_aux);

    for(p = 0; p < pMax; p++) {
      temp1 = abs(scatter_coef[j * 2 * pMax + p.compound]);
      temp1 *= temp1;
      temp2 = abs(scatter_coef[pMax + j * 2 * pMax + p.compound]);
      temp2 *= temp2;
      Cabs += temp1 * Cabs_aux[p.compound] + temp2 * Cabs_aux[pMax + p.compound];
      
    }
  }
 
  delete[] Cabs_aux;
  
  return (1 / (std::real(waveK) * std::real(waveK))) * Cabs;
}


double Result::getScatteringCrossSection_SH(int gran1, int gran2) {
 
  CompoundIterator p, q;
  int pMax = Tools::iteratorMax(nMaxS);
  int qMax = q.max(nMaxS);

  double temp1(0.0);
  
  std::complex<double> **T_AB = new std::complex<double> *[2 * (p.max(nMaxS))];
  
  std::complex<double> *SCcoefpom_SH = new std::complex<double>[2 * pMax];
  
  for(p = 0; p < (int)(2 * p.max(nMaxS)); p++) {
    T_AB[p] = new std::complex<double>[2 * p.max(nMaxS)];
  }

  double mu_b_r = std::real(geometry->bground.mu_r);

  double eps_b_r = std::real(geometry->bground.epsilon_r);

  auto const omega = excitation->omega();
  
  for(int j = gran1; j < gran2; j++) {

   Spherical<double> point_ = geometry->objects[j].vR;
   
   Spherical<double> Rrel = point_ - Spherical<double>(0.0, 0.0, 0.0);
  
   optimet::Coupling const coupling(Rrel,  2.0 * waveK, nMaxS, false);  // paziti na waveK
   
    for(p = 0; p < pMax; p++)   {
    for(q = 0; q < qMax; q++) {
    
      T_AB[p][q] = coupling.diagonal(q, p);
      T_AB[p + pMax][q + qMax] = coupling.diagonal(q, p);  // This is the transfer matrix for the incident wave
      T_AB[p + pMax][q] = coupling.offdiagonal(q, p);
      T_AB[p][q + qMax] = coupling.offdiagonal(q, p);
      
    }
  }
  
  for(p = 0; p < 2*pMax; p++) {

  SCcoefpom_SH[ p ] = scatter_coef_SH[j * 4 * pMax + p.compound] + scatter_coef_SH[2 * pMax + j * 4 * pMax + p.compound] ;
 
  }
  
   for(p = 0; p < 2 * pMax; p++) {
   
    for(q = 0; q < 2 * qMax; q++) {


      temp1 += std::real( (T_AB[p][q]) * std::conj(SCcoefpom_SH[q.compound]) * 
      
        std::conj(T_AB[p][q]) * SCcoefpom_SH[q.compound] );
      

    }
  }  
}

   delete[] SCcoefpom_SH;

  for(p = 0; p < 2 * pMax; p++) {
    delete[] T_AB[p];
  }

  delete[] T_AB;
  
  return  (1.0 / (4.0 * mu_b_r*eps_b_r))* temp1 ; // frequency doubled because of SH
}


 double Result::getAbsorptionCrossSection_SH(std::vector<double *> CLGcoeff, int gran1, int gran2) {


  auto const omega = excitation->omega();
      
  int NO = geometry->objects.size(); // how many particles there is
  
  double absCS (0.0);
  
  CompoundIterator p;
  
  std::complex<double> mu_b = geometry->bground.mu;
  
  std::complex<double> eps_b = geometry->bground.epsilon;
  
  std::complex<double> eta = std::sqrt (mu_b / eps_b);
  
  int pMaxS = p.max(nMaxS);
 
  std::complex<double>  sigma;
 
  int sizeCFsh = gran2-gran1; 

  Vector<t_complex> coeffSH(sizeCFsh);

  geometry->AbsCSSHcoeff(CLGcoeff, gran1, gran2, excitation, internal_coef, internal_coef_SH, nMaxS, coeffSH.data());
  
  int brojac(0);

  int objIndex;
 
  for(int ii = gran1; ii < gran2; ii++) {

   objIndex = ii / pMaxS;

  sigma = - std::complex<double>(0.0, 1.0) * consEpsilon0 * 2.0 * omega * (geometry->objects[objIndex].elmag.epsilon_r_SH - 1.0); 
   
    absCS = absCS +  std::real( (2.0 * eta) * 0.5 * sigma * coeffSH(brojac));

   brojac++;  
  } // for  ii
   
  return  absCS ;
   
  }



     int Result::setFields(std::vector<double> &Rr, std::vector<double> 
                           &Rthe, std::vector<double> &Rphi, bool projection_, std::vector<double *> CLGcoeff) {
  
  Spherical<double> Rloc, Rrel;
  
  int sizeField;
  
  CompoundIterator p;

  int pMaxS = p.max(nMaxS);
 
  std::complex<double> coeffXpl[pMaxS], coeffXmn[pMaxS];
    
  // Calculate the fields Fundamental Frequency and SH Frequency

    SphericalP<std::complex<double>> EField_FF;
    SphericalP<std::complex<double>> HField_FF;
    SphericalP<std::complex<double>> EField_SH;
    SphericalP<std::complex<double>> HField_SH;
    
    std::vector<std::complex<double>> EField_FF_x;
    std::vector<std::complex<double>> EField_FF_y;
    std::vector<std::complex<double>> EField_FF_z;
    std::vector<std::complex<double>> HField_FF_x;
    std::vector<std::complex<double>> HField_FF_y;
    std::vector<std::complex<double>> HField_FF_z;
    
    std::vector<std::complex<double>> EField_SH_x;
    std::vector<std::complex<double>> EField_SH_y;
    std::vector<std::complex<double>> EField_SH_z;
    std::vector<std::complex<double>> HField_SH_x;
    std::vector<std::complex<double>> HField_SH_y;
    std::vector<std::complex<double>> HField_SH_z;
   
    
    for (int ii=0; ii<Rr.size(); ii++){

    Rloc.rrr = Rr[ii];
    Rloc.the = Rthe[ii];
    Rloc.phi = Rphi[ii];

    int intInd = geometry->checkInner(Rloc);
    
    if(intInd >= 0){
    
    Rrel = Tools::toPoint(Rloc, geometry->objects[intInd].vR);
    geometry->COEFFpartSH(intInd, excitation, internal_coef, Rrel.rrr, nMaxS, coeffXmn, coeffXpl, CLGcoeff);
    
    }
   
    getEHFields(Rloc, EField_FF, HField_FF, EField_SH, HField_SH, projection_, coeffXmn, coeffXpl);
    
    EField_FF_x.push_back (EField_FF.rrr);
    EField_FF_y.push_back (EField_FF.the);
    EField_FF_z.push_back (EField_FF.phi);
    HField_FF_x.push_back (HField_FF.rrr);
    HField_FF_y.push_back (HField_FF.the);
    HField_FF_z.push_back (HField_FF.phi);

    EField_SH_x.push_back (EField_SH.rrr);
    EField_SH_y.push_back (EField_SH.the);
    EField_SH_z.push_back (EField_SH.phi);
    HField_SH_x.push_back (HField_SH.rrr);
    HField_SH_y.push_back (HField_SH.the);
    HField_SH_z.push_back (HField_SH.phi);

   
}//for


   sizeField = EField_FF_x.size();
   
   MPI_Send(&EField_FF_x[0], EField_FF_x.size(), MPI_DOUBLE_COMPLEX, 0, 1, MPI_COMM_WORLD);
   MPI_Send(&EField_FF_y[0], EField_FF_y.size(), MPI_DOUBLE_COMPLEX, 0, 2, MPI_COMM_WORLD);
   MPI_Send(&EField_FF_z[0], EField_FF_z.size(), MPI_DOUBLE_COMPLEX, 0, 3, MPI_COMM_WORLD);   
   MPI_Send(&HField_FF_x[0], HField_FF_x.size(), MPI_DOUBLE_COMPLEX, 0, 4, MPI_COMM_WORLD);
   MPI_Send(&HField_FF_y[0], HField_FF_y.size(), MPI_DOUBLE_COMPLEX, 0, 5, MPI_COMM_WORLD);
   MPI_Send(&HField_FF_z[0], HField_FF_z.size(), MPI_DOUBLE_COMPLEX, 0, 6, MPI_COMM_WORLD);

   MPI_Send(&EField_SH_x[0], EField_SH_x.size(), MPI_DOUBLE_COMPLEX, 0, 7, MPI_COMM_WORLD);
   MPI_Send(&EField_SH_y[0], EField_SH_y.size(), MPI_DOUBLE_COMPLEX, 0, 8, MPI_COMM_WORLD);
   MPI_Send(&EField_SH_z[0], EField_SH_z.size(), MPI_DOUBLE_COMPLEX, 0, 9, MPI_COMM_WORLD);
   MPI_Send(&HField_SH_x[0], HField_SH_x.size(), MPI_DOUBLE_COMPLEX, 0, 10, MPI_COMM_WORLD);
   MPI_Send(&HField_SH_y[0], HField_SH_y.size(), MPI_DOUBLE_COMPLEX, 0, 11, MPI_COMM_WORLD);
   MPI_Send(&HField_SH_z[0], HField_SH_z.size(), MPI_DOUBLE_COMPLEX, 0, 12, MPI_COMM_WORLD);
   
   MPI_Send(&sizeField, 1, MPI_INT, 0, 13, MPI_COMM_WORLD);

  return 0;

}

int Result::setFieldsModal(OutputGrid &oEGrid_, OutputGrid &oHGrid_, bool projection_,
                           CompoundIterator p_, int singleComponent_) {
  Spherical<double> Rloc;

  // Calculate the fields
  while(!oEGrid_.gridDone) {
    Rloc = oEGrid_.getPoint();
    oHGrid_.getPoint();
    //std::cout << "Calculating fields for point " << oEGrid_.iterator + 1 << " out of "
              //<< oEGrid_.gridPoints << std::endl;

    SphericalP<std::complex<double>> EField;
    SphericalP<std::complex<double>> HField;

    getEHFieldsModal(Rloc, EField, HField, projection_, p_, singleComponent_);

    oHGrid_.pushDataNext(HField);
    oEGrid_.pushDataNext(EField);
  }
  return 0;
}



CompoundIterator Result::getDominant() {
  CompoundIterator p, q;

  q = 0;

  std::complex<double> TEMax = scatter_coef[0];
  std::complex<double> TMMax = scatter_coef[p.max(nMax)];

  for(p = 0; p < p.max(nMax); p++) {
    if((abs(scatter_coef[p]) > abs(TEMax)) ||
       (abs(scatter_coef[p.max(nMax) + p.compound]) > abs(TMMax))) {
      q = p;

      TEMax = scatter_coef[p];
      TMMax = scatter_coef[p.max(nMax) + p.compound];
    }
  }

  return q;
}

void Result::getEHFieldsContCheck(Spherical<double> R_, SphericalP<std::complex<double>> &EField_,
                                  SphericalP<std::complex<double>> &HField_, bool projection_,
                                  int inside_) {
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

  std::complex<double> iZ = (consCmi / sqrt(geometry->bground.mu / geometry->bground.epsilon));

  int pMax = Tools::iteratorMax(nMax);

  CompoundIterator p;

  // Check for inner point and set to 0
  int intInd = geometry->checkInner(R_);
  intInd = inside_;

  if(intInd < 0) // Outside a sphere
  {
    if(!flagSH) // this a fundamental frequency result - calculate the incoming
                // field
    {
      // Incoming field
      optimet::AuxCoefficients aCoefInc(R_, waveK, 1, nMax);
      for(p = 0; p < p.max(nMax); p++) {
        Einc = Einc + (aCoefInc.M(static_cast<long>(p)) * excitation->dataIncAp[p] +
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
      for(size_t j = 0; j < geometry->objects.size(); j++) {
        Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
        optimet::AuxCoefficients aCoef(Rrel, waveK, 0, nMax);

        for(p = 0; p < pMax; p++) {
          Einc =
              Einc +
              aCoef.M(static_cast<long>(p)) * geometry->objects[j].sourceCoef[static_cast<int>(p)] +
              aCoef.N(static_cast<long>(p)) *
                  geometry->objects[j].sourceCoef[static_cast<int>(p) + pMax];
        }
      }
    }

    // Scattered field
    for(size_t j = 0; j < geometry->objects.size(); j++) {
      SphericalP<std::complex<double>> Efield_local = SphericalP<std::complex<double>>(
          std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
          std::complex<double>(0.0, 0.0));

      Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
      optimet::AuxCoefficients aCoef(Rrel, waveK, 0, nMax);

      for(p = 0; p < p.max(nMax); p++) {
        Efield = Efield + aCoef.M(static_cast<long>(p)) * scatter_coef[j * 2 * pMax + p.compound] +
                 aCoef.N(static_cast<long>(p)) * scatter_coef[pMax + j * 2 * pMax + p.compound];
        Hfield = Hfield +
                 (aCoef.N(static_cast<long>(p)) * scatter_coef[j * 2 * pMax + p.compound] +
                  aCoef.M(static_cast<long>(p)) * scatter_coef[pMax + j * 2 * pMax + p.compound]) *
                     iZ;
      }
    }
  } else // Inside a sphere
  {
    Rrel = Tools::toPoint(R_, geometry->objects[intInd].vR);
    optimet::AuxCoefficients aCoef(Rrel, waveK * sqrt(geometry->objects[intInd].elmag.epsilon_r *
                                                      geometry->objects[intInd].elmag.mu_r),
                                   1, nMax);

    std::complex<double> iZ_object = (consCmi / sqrt(geometry->objects[intInd].elmag.mu /
                                                     geometry->objects[intInd].elmag.epsilon));

    for(p = 0; p < p.max(nMax); p++) {
      Efield =
          Efield +
          (aCoef.M(static_cast<long>(p)) * internal_coef[intInd * 2 * pMax + p.compound] +
           aCoef.N(static_cast<long>(p)) * internal_coef[pMax + intInd * 2 * pMax + p.compound]);
      Hfield =
          Hfield +
          (aCoef.N(static_cast<long>(p)) * internal_coef[intInd * 2 * pMax + p.compound] +
           aCoef.M(static_cast<long>(p)) * internal_coef[pMax + intInd * 2 * pMax + p.compound]) *
              iZ_object;
    }
  }

  if(projection_) {
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
  auto const projection = true; // Spherical projection - True - projection is internally
                                // set to be evaluated w.r.t. object[0]
  int outside = -1;             // Forces result to be outside an object
  int inside = 0;               // Forces result to be inside an object
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
  int max_ii = 180; // theta observation range (increment by 1 degree) - [1, max_ii-1]
  int max_jj = 180; // phi   observation range (increment by 1 degree) - [1, max_jj-1]
  for(int ii = 1; ii <= max_ii - 1; ii++) {
    for(int jj = 1; jj <= max_jj - 1; jj++) {

      APoint =
          Spherical<double>(radius, consPi * (double(ii) / 180.), consPi * (double(jj) / 180.));
      getEHFieldsContCheck(APoint, AnEField_out, AnHField_out, projection, outside);
      getEHFieldsContCheck(APoint, AnEField_in, AnHField_in, projection, inside);
      // op -------------------------------------------------------------------
      std::cout << "Continuity check : computed " << jj + ((ii - 1) * (max_ii - 1))
                << " out of a total of " << (max_ii - 1) * (max_jj - 1) << std::endl;
      // EF
      E1_err_mag << (abs(AnEField_out.rrr) - abs(AnEField_in.rrr * eps_r)) /
                        abs(AnEField_in.rrr * eps_r)
                 << " ";
      E2_err_mag << (abs(AnEField_out.the) - abs(AnEField_in.the)) / abs(AnEField_in.the) << " ";
      E3_err_mag << (abs(AnEField_out.phi) - abs(AnEField_in.phi)) / abs(AnEField_in.phi) << " ";
      // EF
      H1_err_mag << (abs(AnHField_out.rrr) - abs(AnHField_in.rrr * mu_r)) /
                        abs(AnHField_in.rrr * mu_r)
                 << " ";
      H2_err_mag << (abs(AnHField_out.the) - abs(AnHField_in.the)) / abs(AnHField_in.the) << " ";
      H3_err_mag << (abs(AnHField_out.phi) - abs(AnHField_in.phi)) / abs(AnHField_in.phi) << " ";
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
