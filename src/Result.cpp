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
    :result_FF(nullptr) {
     
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
  result_FF = NULL;
  result_SH = NULL;
  nMax = geometry_->nMax();
  nMaxS = geometry_->nMaxS();

  CLGcoeff.resize(9); 
  scatter_coef.resize(2 * Tools::iteratorMax(nMax) * geometry->objects.size());
  internal_coef.resize(2 * Tools::iteratorMax(nMax) * geometry->objects.size());
  
  if (geometry->objects[0].scatterer_type == "sphere"){
  scatter_coef_SH.resize(4 * Tools::iteratorMax(nMaxS) * geometry->objects.size());
  internal_coef_SH.resize(4 * Tools::iteratorMax(nMaxS) * geometry->objects.size());
 }

else if (geometry->objects[0].scatterer_type == "arbitrary.shape"){
  scatter_coef_SH.resize(2 * Tools::iteratorMax(nMaxS) * geometry->objects.size());
  internal_coef_SH.resize(2 * Tools::iteratorMax(nMaxS) * geometry->objects.size());
  }

}

void Result::update(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> excitation_) {
  geometry = geometry_;
  excitation = excitation_;
  waveK = excitation->waveK;
}

void Result::init(std::shared_ptr<Geometry> geometry_, std::shared_ptr<Excitation> excitation_,
                  Result *result_FF_) {

  init(geometry_, excitation_);
  result_FF = result_FF_;
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
      optimet::AuxCoefficients aCoefInc(R_, waveK, 1, nMax); //regular VSWFs
      
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

 
      optimet::AuxCoefficients aCoefFF(Rrel, waveK, 0, nMax); //radiative VSWFs

       
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
  if(excitation->SH_cond){  
  for(size_t j = 0; j < geometry->objects.size(); j++) {


      Rrel = Tools::toPoint(R_, geometry->objects[j].vR);
      
      optimet::AuxCoefficients aCoefSH(Rrel, std::complex<double>(2.0, 0.0) * waveK, 0, nMaxS); //radiative spherical

      for(p = 0; p < p.max(nMaxS); p++) {
      
      if (geometry->objects[j].scatterer_type == "sphere"){   
      bmnSH = scatter_coef_SH[j * 2 * pMaxS + p.compound];
      
      amnSH = scatter_coef_SH[j * 2 * pMaxS + pMaxS + p.compound];
     }

      else if (geometry->objects[j].scatterer_type == "arbitrary.shape"){ 
      bmnSH =  (1.0/waveK_0) * scatter_coef_SH[j * 2 * pMaxS + p.compound];
      
      amnSH =  (1.0/waveK_0) * scatter_coef_SH[j * 2 * pMaxS + pMaxS + p.compound];
      }

        Efield_SH = Efield_SH + (aCoefSH.M(static_cast<long>(p)) * (bmnSH) +
                                 aCoefSH.N(static_cast<long>(p)) * (amnSH)) * waveK_0; 
        
        Hfield_SH = Hfield_SH +
                 (aCoefSH.N(static_cast<long>(p)) * (bmnSH) +
                  aCoefSH.M(static_cast<long>(p)) * (amnSH)) * (iZ * waveK_0);
                      
      }
      
    }
    
   }

  } // if
 
 
  else // Inside a sphere
 
    {
     // this is FF result  
    Rrel = Tools::toPoint(R_, geometry->objects[intInd].vR);
    
    optimet::AuxCoefficients aCoefFF(Rrel, waveK_0 * sqrt(geometry->objects[intInd].elmag.epsilon_r *
                                                      geometry->objects[intInd].elmag.mu_r), 1, nMax); //regular VSWFs
                                                      
                                   
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
  if(excitation->SH_cond){  
    optimet::AuxCoefficients aCoefSH(Rrel, std::complex<double> (2.0 , 0.0) * waveK_0 * sqrt(geometry->objects[intInd].elmag.epsilon_r_SH * geometry->objects                                      [intInd].elmag.mu_r_SH), 1, nMaxS); // regular VSWFs
                                   
                                   
    std::complex<double> iZ_object_SH = (consCmi / sqrt(geometry->objects[intInd].elmag.mu_SH /
                                                     geometry->objects[intInd].elmag.epsilon_SH));
                                                                                                                                                       

    for(p = 0; p < p.max(nMaxS); p++) {

    if (geometry->objects[intInd].scatterer_type == "sphere"){
    cmnSH = internal_coef_SH[intInd * 2 * pMaxS + p.compound];

    dmnSH = internal_coef_SH[pMaxS + intInd * 2 * pMaxS + p.compound];
    }

   if (geometry->objects[intInd].scatterer_type == "arbitrary.shape"){
    cmnSH =  (1.0/waveK_0) * internal_coef_SH[intInd * 2 * pMaxS + p.compound];
    
    dmnSH =  (1.0/waveK_0) * internal_coef_SH[pMaxS + intInd * 2 * pMaxS + p.compound];
    }
    
     Efield_SH =
          Efield_SH +
          (aCoefSH.M(static_cast<long>(p)) * (cmnSH) +
           aCoefSH.N(static_cast<long>(p)) * (dmnSH)) * (waveK_0);
           
      Egamma_SH = Egamma_SH +
             (aCoefSH.Xm(static_cast<long>(p)) * coeffXmn[p] +
            aCoefSH.Xp(static_cast<long>(p)) * coeffXpl[p]);
        
      Hfield_SH =
          Hfield_SH +
          (aCoefSH.N(static_cast<long>(p)) * (cmnSH) +
           aCoefSH.M(static_cast<long>(p)) * (dmnSH)) * iZ_object_SH * waveK_0;
              

   } 
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


double Result::getScatteringCrossSection(int gran1, int gran2) {
 
  CompoundIterator p, q;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);

  double temp1(0.0);
  
  std::complex<double> **T_AB = new std::complex<double> *[2 * (p.max(nMax))];
  
  std::complex<double> *SCcoefpom = new std::complex<double>[2 * pMax];
  
  for(p = 0; p < (int)(2 * p.max(nMax)); p++) {
    T_AB[p] = new std::complex<double>[2 * p.max(nMax)];
  }

  auto const omega = excitation->omega();
  
  int NO = geometry->objects.size();
  
  for(int j = gran1; j < gran2; j++) {

   Spherical<double> point_ = geometry->objects[j].vR;
   
   Spherical<double> Rrel = point_ - Spherical<double>(0.0, 0.0, 0.0);
  
   optimet::Coupling const coupling(Rrel,  waveK, nMax, false); 
   
    for(p = 0; p < pMax; p++)   {
    for(q = 0; q < qMax; q++) {
    
      T_AB[p][q] = coupling.diagonal(q, p);
      T_AB[p + pMax][q + qMax] = coupling.diagonal(q, p);
      T_AB[p + pMax][q] = coupling.offdiagonal(q, p);
      T_AB[p][q + qMax] = coupling.offdiagonal(q, p);
      
    }
  }
  
  for(p = 0; p < 2 * pMax; p++) {
  
  SCcoefpom[ p ] = scatter_coef[j * 2 * pMax + p.compound] ;
  
  }
  
   for(p = 0; p < 2 * pMax; p++) {
   
    for(q = 0; q < 2 * qMax; q++) {


      temp1 += std::real( (T_AB[p][q]) * std::conj(SCcoefpom[q.compound]) * 
      
        std::conj(T_AB[p][q]) * SCcoefpom[q.compound] );
      

    }
  }  
}

   delete[] SCcoefpom;

  for(p = 0; p < 2 * pMax; p++) {
    delete[] T_AB[p];
  }

  delete[] T_AB;
  
  return  (1. / (std::real(waveK) * std::real(waveK)))* temp1 ;
}



double Result::getAbsorptionCrossSection(int gran1, int gran2) {
 
  CompoundIterator p;
  int pMax = p.max(nMax);

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
  int pMax = p.max(nMaxS);
  int qMax = q.max(nMaxS);

  double temp1(0.0), ArbCf;
  
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
  
   optimet::Coupling const coupling(Rrel,  2.0 * waveK, nMaxS, false);  // waveK is doubled
   
    for(p = 0; p < pMax; p++)   {
    for(q = 0; q < qMax; q++) {
    
      T_AB[p][q] = coupling.diagonal(q, p);
      T_AB[p + pMax][q + qMax] = coupling.diagonal(q, p);
      T_AB[p + pMax][q] = coupling.offdiagonal(q, p);
      T_AB[p][q + qMax] = coupling.offdiagonal(q, p);
      
    }
  }
  
  for(p = 0; p < 2*pMax; p++) {
  
  if (geometry->objects[0].scatterer_type == "sphere"){
  SCcoefpom_SH[ p ] = scatter_coef_SH[j * 2 * pMax + p.compound];
  ArbCf = eps_b_r * mu_b_r;
  }

  else if (geometry->objects[0].scatterer_type == "arbitrary.shape"){
  SCcoefpom_SH[ p ] = scatter_coef_SH[j * 2 * pMax + p.compound];
  
  ArbCf = std::real(waveK) * std::real(waveK);
  }
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
  
  return  (1.0 / (4.0 * ArbCf))* temp1 ; // frequency doubled because of SH
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

  int objIndex, brojac(0);
 
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
   
}
