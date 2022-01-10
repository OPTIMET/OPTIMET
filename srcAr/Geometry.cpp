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
#include "AuxCoefficients.h"
#include "Coupling.h"
#include "Geometry.h"
#include "Trian.h"
#include "Algebra.h"
#include "Bessel.h"
#include "CompoundIterator.h"
#include "HarmonicsIterator.h"
#include "Symbol.h"
#include "Tools.h"
#include "Types.h"
#include "constants.h"

#include <cmath>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>
#include <chrono>
using namespace std::chrono;

Geometry::~Geometry() {}
Geometry::Geometry() {}

void Geometry::pushObject(Scatterer const &object_) {
  for(auto const &obj : objects)
    if(Tools::findDistance(obj.vR, object_.vR) <= (object_.radius + obj.radius)) {
      std::ostringstream sstr;
      sstr << "The particle at (" << Tools::toCartesian(object_.vR).x << ", "
           << Tools::toCartesian(object_.vR).y << ", " << Tools::toCartesian(object_.vR).z << ") "
           << "overlaps with the one at (" << Tools::toCartesian(obj.vR).x << ", "
           << Tools::toCartesian(obj.vR).y << ", " << Tools::toCartesian(obj.vR).z << "), "
           << "with circumscribed radii " << object_.radius << " and " << obj.radius;
      throw std::runtime_error(sstr.str());
    }
  objects.emplace_back(object_);
}

bool Geometry::is_valid() const {
  using namespace optimet;
  if(objects.size() == 0)
    return false;
  for(t_uint i(0); i < objects.size(); ++i)
    for(t_uint j(0); j < objects.size(); ++j)
      if(Tools::findDistance(objects[i].vR, objects[j].vR) <= objects[i].radius + objects[j].radius)
        return false;
  return true;
}

void Geometry::initBground(ElectroMagnetic bground_) { bground = bground_; }

optimet::t_uint Geometry::scatterer_size() const {
  auto const object_size = [](optimet::t_uint current, Scatterer const &scatterer) {
    return current + 2 * scatterer.nMax * (scatterer.nMax + 2);
  };
  return std::accumulate(objects.cbegin(), objects.cend(), 0, object_size);
}

bool Geometry::no_overlap(Scatterer const &object_) {
  if(objects.size() < 1)
    return true;

  for(auto const &obj : objects)
    if(Tools::findDistance(obj.vR, object_.vR) <= (object_.radius + obj.radius))
      return false;

  return true;
}


int Geometry::checkInner(Spherical<double> R_) {
  Spherical<double> Rrel;
  Cartesian <double> Rcar, centCar;
  std::vector<double> RtoCp(3);
  double DOT;

  if (objects[0].scatterer_type == "sphere"){
  for(size_t j = 0; j < objects.size(); j++) {
    // Translate to object j
    Rrel = Tools::toPoint(R_, objects[j].vR);

    // Check if R_.rrr is inside radius of object j
    if(Rrel.rrr <= objects[j].radius)
      return j;
  }

  return -1;
}

else if (objects[0].scatterer_type == "arbitrary.shape"){
  Rcar = Tools::toCartesian(R_);
  
  for(size_t j = 0; j < objects.size(); j++) {
  
  int clearPass(1);
  centCar = Tools::toCartesian(objects[j].vR);
  unsigned int Nt = (objects[j].getTopolsize()) / 3;
  
  for (int ele1 = 0; ele1 < Nt; ++ele1){
	        
	        const int* n1 = objects[j].getNOvertex(ele1);
	
		const double* p1 = objects[j].getCoord(n1[0]);

		const double* p2 = objects[j].getCoord(n1[1]);
	
		const double* p3 = objects[j].getCoord(n1[2]);
                
                Trian trian(p1, p2, p3); 
               
		std::vector<double> nvec = trian.getnorm();
		std::vector<double> cp = trian.getcp();

		
		cp[0] = cp[0] + centCar.x;
		cp[1] = cp[1] + centCar.y;
		cp[2] = cp[2] + centCar.z;
		
		RtoCp[0] = Rcar.x - cp[0];
		RtoCp[1] = Rcar.y - cp[1];
		RtoCp[2] = Rcar.z - cp[2];
		
       DOT = Tools::dot(&nvec[0], &RtoCp[0]);
         
         if ( DOT >= 0 ){  // the point is outside the particle
         clearPass = 0;
	 break;	
	 }
	 
	 
    }
      	
      if (clearPass)
      return j;		
  }
  
  return -1;
  
  }
}

#ifdef OPTIMET_MPI
void Geometry::Coefficients(int nMax, int nMaxS, std::vector<double *> CLGcoeff, int gran1, int gran2){ 

   
   optimet::symbol::C_10m1coeff(CLGcoeff[0], nMax, nMaxS, gran1, gran2);
   optimet::symbol::C_11m1coeff(CLGcoeff[1], nMax, nMaxS, gran1, gran2);
   optimet::symbol::C_00m1coeff(CLGcoeff[2], nMax, nMaxS, gran1, gran2);
   optimet::symbol::C_01m1coeff(CLGcoeff[3], nMax, nMaxS, gran1, gran2);
   
   optimet::symbol::W_m1m1coeff(CLGcoeff[4], nMax, nMaxS, gran1, gran2);
   optimet::symbol::W_11coeff(CLGcoeff[5], nMax, nMaxS, gran1, gran2);
   optimet::symbol::W_00coeff(CLGcoeff[6], nMax, nMaxS, gran1, gran2);
   optimet::symbol::W_10coeff(CLGcoeff[7], nMax, nMaxS, gran1, gran2);
   optimet::symbol::W_01coeff(CLGcoeff[8], nMax, nMaxS, gran1, gran2);
}
#endif

int Geometry::COEFFpartSH(int objectIndex_, std::shared_ptr<optimet::Excitation const> incWave_, 
optimet::Vector<optimet::t_complex> &internalCoef_FF_, double r,
        int nMaxS_, std::complex<double> *coefXmn, std::complex<double> *coefXpl, std::vector<double *> CLGcoeff) {
                 
   double *W_m1m1 = CLGcoeff[4];
   double *W_11 = CLGcoeff[5];
   double *W_00 = CLGcoeff[6];                                
                 
   CompoundIterator p;
  
   int pMax = p.max(nMaxS_);
   int nMax_ = this->nMax();   

  auto const omega = incWave_->omega();

 
    
    for(p = 0; p < pMax; p++) {  
         
             
    coefXmn[p] = optimet::symbol::CXm1(p, W_m1m1, W_00, W_11, nMax_,
                  internalCoef_FF_, r,
                   objectIndex_,
                  omega, objects[objectIndex_]);
                  
                  
                 
    coefXpl[p] = optimet::symbol::CXp1(p, W_m1m1, W_00, W_11, nMax_,
                  internalCoef_FF_, r,
                   objectIndex_,
                  omega, objects[objectIndex_]);              
                  
                   
     
    }

  return 0;
}

#ifdef OPTIMET_MPI
void Geometry::getEXCvecSH_ARB3_parall(optimet::Vector<optimet::t_complex>& EXvec, std::shared_ptr<optimet::Excitation const> excitation, optimet::Vector<optimet::t_complex> &externalCoef_FF_, optimet::Vector<optimet::t_complex> &internalCoef_FF_, int gran1, int gran2, int objIndex){
  using namespace optimet;
  
  int nMax = this->nMax();
  int nMaxS = this->nMaxS();
  auto const N = nMaxS * (nMaxS + 2);
  
          CompoundIterator mu1, mup, pp;
         
          optimet::t_real omega_ = excitation->omega();
          Cartesian<double> intpoinCar; // integration point in Cartesian system
          Spherical<double> intpoinSph; // integration point in spherical system
          
          SphericalP<std::complex<double>> Eint_FF;
          
          int muMax =mu1.max(nMaxS);
          int ppMax = pp.max(nMax);
          int brojac(0);
          int size = gran2 - gran1;
          
          std::complex<double>* Etan = new std::complex<double>[3];
          std::complex<double>* Tan = new std::complex<double>[3];
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> Enorm, resDOT1, resDOT2, EE;
          
         const std::complex<double> k_0_SH = 2 * (omega_) * std::sqrt(consEpsilon0 * consMu0); 
         auto const k_s_SH = 2 * omega_ * std::sqrt(objects[objIndex].elmag.epsilon_SH * objects[objIndex].elmag.mu_SH); 
         auto const k_b_SH = 2 * omega_ * std::sqrt(bground.epsilon * bground.mu);
         auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu);
         auto const k_s = omega_ * std::sqrt(objects[objIndex].elmag.epsilon * objects[objIndex].elmag.mu);
        
         const std::complex<double> ksiparppar = objects[objIndex].elmag.ksiparppar; // SH surf tensor coefficient
         const std::complex<double> ksippp = objects[objIndex].elmag.ksippp; // SH surf tensor coefficient
         const std::complex<double> gamma = objects[objIndex].elmag.gamma; // SH bulk tensor coefficient

         auto const Cf = (k_b_SH * consEpsilon0 * gamma)/objects[objIndex].elmag.epsilon_SH;
  
        unsigned int Nt = (objects[objIndex].getTopolsize()) / 3;  //number of triangles
         
        std::vector<std::vector<double>> Points = Tools::getPoints4();

        std::vector<double> Weights = Tools::getWeights4();
        
	for (mu1 = gran1; mu1 < gran2 ; mu1++){
	
	std::complex<double> I11s(0.0,0.0), I21s(0.0,0.0), I31s(0.0,0.0), I12s(0.0,0.0), I22s(0.0,0.0), I32s(0.0,0.0);
        	
	int nn = mu1.first;
	int mm =  mu1.second;
	mup = nn * (nn + 1) + mm - 1;
	
       	
	for (int ele1 = 0; ele1 < Nt; ++ele1){
	        
	        const int* n1 = objects[objIndex].getNOvertex(ele1);
	
		const double* p1 = objects[objIndex].getCoord(n1[0]);

		const double* p2 = objects[objIndex].getCoord(n1[1]);
	
		const double* p3 = objects[objIndex].getCoord(n1[2]);
                
                Trian trian(p1, p2, p3); 
                
		double det = trian.getDeter();
		std::vector<double> nvec = trian.getnorm();
		
	// surface integration	
	for (int ni = 0; ni != Weights.size(); ++ni) {
      
      Eint_FF = SphericalP<std::complex<double>>(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
	
	double N1 = Points[ni][0];
	double N2 = Points[ni][1];
        double N0 = 1.0 - N1 - N2;

	double wi = Weights[ni];
         
	intpoinCar.x = (p1[0])*N1 + (p2[0])*N2 + (p3[0])*N0;
	intpoinCar.y = (p1[1])*N1 + (p2[1])*N2 + (p3[1])*N0;
	intpoinCar.z = (p1[2])*N1 + (p2[2])*N2 + (p3[2])*N0;

	intpoinSph = Tools::toSpherical(intpoinCar); 
	
	// Internal field

      optimet::AuxCoefficients aCoefInt(intpoinSph, k_s, 1, nMax); //internal FF field
      
      
      for(pp = 0; pp < ppMax; pp++) {

        Eint_FF = Eint_FF + (aCoefInt.M(static_cast<long>(pp)) * internalCoef_FF_[objIndex*2*ppMax + pp] +
                       aCoefInt.N(static_cast<long>(pp)) * internalCoef_FF_[objIndex*2*ppMax + ppMax + pp]);
                                                                  
      }
      //tangential component -nxnxEext
      Tools::crossTanTr(Etan, resCR, &nvec[0], Eint_FF);
      
      // inner component of normal field at FF
      Enorm = Tools::dot (&nvec[0], Eint_FF);
                        
       Tools::crossTan(Tan, resCR, &nvec[0], &Etan[0]);
                                 
       optimet::AuxCoefficients aCoefext3(intpoinSph, k_b_SH, 0, nMaxS); //radiative VSWFs (3)
                                                    
       resDOT1 = Tools::dot (&Tan[0], aCoefext3.M(static_cast<long>(mup)));
       resDOT2 = Tools::dot (&Tan[0], aCoefext3.N(static_cast<long>(mup))); 

       I11s = I11s + wi * det * 2.0 * Enorm * ksiparppar * k_0_SH * k_0_SH * resDOT1;
      I12s = I12s + wi * det * 2.0 * Enorm * ksiparppar * k_0_SH * k_0_SH * resDOT2;
      
      resDOT1 = Tools::dot (&nvec[0], aCoefext3.M(static_cast<long>(mup)));
      resDOT2 = Tools::dot (&nvec[0], aCoefext3.N(static_cast<long>(mup)));
      
      I21s = I21s - wi * det * k_b_SH * k_b_SH *  ksippp * Enorm*Enorm * resDOT1;
      I22s = I22s - wi * det * k_b_SH * k_b_SH *  ksippp * Enorm*Enorm * resDOT2 ;
                                                                             
      
      // calculation of the particular solution on the surface 
     
     EE = Etan[0] * Etan[0] + Etan[1] * Etan[1] + Etan[2] * Etan[2] + Enorm*Enorm;
     
     resDOT1 = Tools::dot (&nvec[0], aCoefext3.M(static_cast<long>(mup)));
     resDOT2 = Tools::dot (&nvec[0], aCoefext3.N(static_cast<long>(mup)));

     I31s = I31s - ( wi * det * k_b_SH * Cf) * EE * resDOT1;
     I32s = I32s - ( wi * det * k_b_SH * Cf) * EE * resDOT2; 
  
  }// surface integration over the triangle is OVER
  
  } // all the triangles surface integration
 

  EXvec (brojac) = I11s + I21s + I31s;
  EXvec (brojac + size) = I12s + I22s + I32s;
  brojac++;
  } // outer for loop
 
  }

  void Geometry::getEXCvecSH_ARB1_parall(optimet::Vector<optimet::t_complex>& EXvec, std::shared_ptr<optimet::Excitation const> excitation, optimet::Vector<optimet::t_complex> &externalCoef_FF_, optimet::Vector<optimet::t_complex> &internalCoef_FF_, int gran1, int gran2, int objIndex){
  using namespace optimet;
  
  int nMax = this->nMax();
  int nMaxS = this->nMaxS();
  auto const N = nMaxS * (nMaxS + 2);
  
          CompoundIterator mu1, mup, pp;
         
          optimet::t_real omega_ = excitation->omega();
          Cartesian<double> intpoinCar; // integration point in Cartesian system
          Spherical<double> intpoinSph; // integration point in spherical system
          
          SphericalP<std::complex<double>> Eint_FF;
          
          int muMax =mu1.max(nMaxS);
          int ppMax = pp.max(nMax);
          int brojac(0);
          int size = gran2 - gran1;
          
          std::complex<double>* Etan = new std::complex<double>[3];
          std::complex<double>* Tan = new std::complex<double>[3];
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> Enorm, resDOT1, resDOT2, EE;
          
         const std::complex<double> k_0_SH = 2 * (omega_) * std::sqrt(consEpsilon0 * consMu0); 
         auto const k_s_SH = 2 * omega_ * std::sqrt(objects[objIndex].elmag.epsilon_SH * objects[objIndex].elmag.mu_SH); 
         auto const k_b_SH = 2 * omega_ * std::sqrt(bground.epsilon * bground.mu);
         auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu);
         auto const k_s = omega_ * std::sqrt(objects[objIndex].elmag.epsilon * objects[objIndex].elmag.mu);
 
         const std::complex<double> ksiparppar = objects[objIndex].elmag.ksiparppar; // SH surface tensor coefficient
         const std::complex<double> ksippp = objects[objIndex].elmag.ksippp; // SH surface tensor coefficient
         const std::complex<double> gamma = objects[objIndex].elmag.gamma; // SH surface tensor coefficient

         auto const Cf = (k_b_SH * consEpsilon0 * gamma)/objects[objIndex].elmag.epsilon_SH;

  
        unsigned int Nt = (objects[objIndex].getTopolsize()) / 3;  //number of triangles
        
        std::vector<std::vector<double>> Points = Tools::getPoints4();

	std::vector<double> Weights = Tools::getWeights4();

	for (mu1 = gran1; mu1 < gran2 ; mu1++){
      
	std::complex<double> I11s(0.0,0.0), I21s(0.0,0.0), I31s(0.0,0.0), I12s(0.0,0.0), I22s(0.0,0.0), I32s(0.0,0.0);
        	
	int nn = mu1.first;
	int mm =  mu1.second;
	mup = nn * (nn + 1) + mm - 1;
	
	
	for (int ele1 = 0; ele1 < Nt; ++ele1){
	        
	        const int* n1 = objects[objIndex].getNOvertex(ele1);
	
		const double* p1 = objects[objIndex].getCoord(n1[0]);

		const double* p2 = objects[objIndex].getCoord(n1[1]);
	
		const double* p3 = objects[objIndex].getCoord(n1[2]);
                
                Trian trian(p1, p2, p3); 
                
		double det = trian.getDeter();
		std::vector<double> nvec = trian.getnorm();
		
	// surface integration	
	for (int ni = 0; ni != Weights.size(); ++ni) {
      
      Eint_FF = SphericalP<std::complex<double>>(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
	
	double N1 = Points[ni][0];
	double N2 = Points[ni][1];
        double N0 = 1.0 - N1 - N2;

	double wi = Weights[ni];
         
	intpoinCar.x = (p1[0])*N1 + (p2[0])*N2 + (p3[0])*N0;
	intpoinCar.y = (p1[1])*N1 + (p2[1])*N2 + (p3[1])*N0;
	intpoinCar.z = (p1[2])*N1 + (p2[2])*N2 + (p3[2])*N0;
        
      	
	intpoinSph = Tools::toSpherical(intpoinCar); 
	
	// Internal fields

      optimet::AuxCoefficients aCoefInt(intpoinSph, k_s, 1, nMax); //internal FF field


     
      for(pp = 0; pp < ppMax; pp++) {

        Eint_FF = Eint_FF + (aCoefInt.M(static_cast<long>(pp)) * internalCoef_FF_[objIndex*2*ppMax + pp] +
                       aCoefInt.N(static_cast<long>(pp)) * internalCoef_FF_[objIndex*2*ppMax + ppMax + pp]);
                          
                                          
      }
            
      //tangential component -nxnxEext
      Tools::crossTanTr(Etan, resCR, &nvec[0], Eint_FF);
 
      // inner component of normal field at FF
      Enorm = Tools::dot (&nvec[0], Eint_FF);
      
      Tools::crossTan(Tan, resCR, &nvec[0], &Etan[0]);
      
      optimet::AuxCoefficients aCoefext1(intpoinSph, k_b_SH, 1, nMaxS); //incoming VSWFs (1)
     
      resDOT1 = Tools::dot (&Tan[0], aCoefext1.M(static_cast<long>(mup)));
      resDOT2 = Tools::dot (&Tan[0], aCoefext1.N(static_cast<long>(mup)));
      
      I11s = I11s + wi * det * 2.0 * Enorm * ksiparppar * k_0_SH * k_0_SH * resDOT1;
      I12s = I12s + wi * det * 2.0 * Enorm * ksiparppar * k_0_SH * k_0_SH * resDOT2;
      
      resDOT1 = Tools::dot (&nvec[0], aCoefext1.M(static_cast<long>(mup)));
      resDOT2 = Tools::dot (&nvec[0], aCoefext1.N(static_cast<long>(mup)));
      
       I21s = I21s - wi * det * k_b_SH * k_b_SH * ksippp * Enorm*Enorm * resDOT1;
      I22s = I22s - wi * det * k_b_SH * k_b_SH * ksippp * Enorm*Enorm * resDOT2 ;
                                                                             
      // calculation of the particular solution on the surface 
      EE =  Etan[0] * Etan[0] + Etan[1] * Etan[1] + Etan[2] * Etan[2] + Enorm*Enorm;
      
     resDOT1 = Tools::dot (&nvec[0], aCoefext1.M(static_cast<long>(mup)));
     resDOT2 = Tools::dot (&nvec[0], aCoefext1.N(static_cast<long>(mup)));

     I31s = I31s - ( wi * det * k_b_SH * Cf) * EE * resDOT1;
     I32s = I32s - ( wi * det * k_b_SH * Cf) * EE * resDOT2; 
  
  }// surface integration over the triangle is OVER
  
  } // all the triangles surface integration
  

  EXvec (brojac) = I11s + I21s + I31s;
  EXvec (brojac + size) = I12s + I22s + I32s;
  brojac++;
  }
  }
#endif

void Geometry::update(std::shared_ptr<optimet::Excitation const> incWave_) {
  // Update the ElectroMagnetic properties of each object
  for(auto &object : objects)
    object.elmag.update(incWave_->lambda());
}

