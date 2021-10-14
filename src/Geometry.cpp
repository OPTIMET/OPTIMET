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

Geometry::~Geometry() {}
Geometry::Geometry() {}

void Geometry::pushObject(Scatterer const &object_) {
  for(auto const &obj : objects)
    if(Tools::findDistance(obj.vR, object_.vR) <= (object_.radius + obj.radius)) {
      std::ostringstream sstr;
      sstr << "The sphere at (" << Tools::toCartesian(object_.vR).x << ", "
           << Tools::toCartesian(object_.vR).y << ", " << Tools::toCartesian(object_.vR).z << ") "
           << "overlaps with the one at (" << Tools::toCartesian(obj.vR).x << ", "
           << Tools::toCartesian(obj.vR).y << ", " << Tools::toCartesian(obj.vR).z << "), "
           << "with radii " << object_.radius << " and " << obj.radius;
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



int Geometry::getCabsAux(double omega_, int objectIndex_, int nMax_, double *Cabs_aux_) {

  std::complex<double> temp1(0., 0.);
  double temp2(0.);
  std::complex<double> k_s =
      omega_ * std::sqrt(objects[objectIndex_].elmag.epsilon * objects[objectIndex_].elmag.mu);
  std::complex<double> k_b = omega_ * std::sqrt(bground.epsilon * bground.mu);

  std::complex<double> rho = k_s / k_b;

  std::complex<double> r_0 = k_b * objects[objectIndex_].radius;

  std::complex<double> mu_j = objects[objectIndex_].elmag.mu;
  std::complex<double> mu_0 = bground.mu;

  std::complex<double> psi(0., 0.), ksi(0., 0.);
  std::complex<double> dpsi(0., 0.), dksi(0., 0.);
  std::complex<double> psirho(0., 0.);
  std::complex<double> dpsirho(0., 0.);

  std::vector<std::complex<double>> J_n_data, J_n_ddata;
  std::vector<std::complex<double>> Jrho_n_data, Jrho_n_ddata;
  
 
  
  try {
    std::tie(J_n_data, J_n_ddata) = optimet::bessel<optimet::Bessel>(r_0, nMax_);
    std::tie(Jrho_n_data, Jrho_n_ddata) = optimet::bessel<optimet::Bessel>(rho * r_0, nMax_);
  } catch(std::range_error &e) { 
    std::cerr << e.what() << std::endl;
    return 1;
  }

  CompoundIterator p;

  int pMax = p.max(nMax_);

  for(p = 0; (int)p < pMax; p++) {
    // obtain Riccati-Bessel functions
    psi = r_0 * J_n_data[p.first];
    dpsi = r_0 * J_n_ddata[p.first] + J_n_data[p.first];

    psirho = r_0 * rho * Jrho_n_data[p.first];
    dpsirho = r_0 * rho * Jrho_n_ddata[p.first] + Jrho_n_data[p.first];

    // Stout 2002 - from scattered
    // TE Part
    temp1 = std::complex<double>(0., 1.) * rho * mu_0 * conj(mu_j) * conj(psirho) * dpsirho;
    temp2 = abs((mu_j * psirho * dpsi - mu_0 * rho * dpsirho * psi));
    temp2 *= temp2;
    Cabs_aux_[p] = real(temp1) / temp2;

    // TM part
    temp1 = std::complex<double>(0., 1.) * conj(rho) * mu_0 * mu_j * conj(psirho) * dpsirho;
    temp2 = abs((mu_0 * rho * psirho * dpsi - mu_j * dpsirho * psi));
    temp2 *= temp2;
    Cabs_aux_[(int)p + pMax] = real(temp1) / temp2;
  }
  return 0;
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



// Second harmonic sources

int Geometry::getIncLocalSH(std::vector<double *> CLGcoeff, int objectIndex_, 
                         std::shared_ptr<optimet::Excitation const> incWave_,
                       optimet::Vector<optimet::t_complex> &internalCoef_FF_, int nMaxS_, 
                       std::complex<double> *Inc_local) {
  
  
   double *C_10m1 = CLGcoeff[0];
   double *C_11m1 = CLGcoeff[1];
   double *C_00m1 = CLGcoeff[2];
   double *C_01m1 = CLGcoeff[3];
   
   double *W_m1m1 = CLGcoeff[4];
   double *W_11 = CLGcoeff[5];
   double *W_00 = CLGcoeff[6];
   double *W_10 = CLGcoeff[7];
   double *W_01 = CLGcoeff[8];
  
                               
  CompoundIterator p;
  
  int pMax = p.max(nMaxS_);
  int nMax_ = this->nMax();
  
  auto const omega = incWave_->omega();
  
    
    for(p = 0; p < pMax; p++) {  
     
    
      Inc_local[p] =
                optimet::symbol::vp_mn(p, C_00m1, C_01m1, nMax_,
                                              internalCoef_FF_,
                                              objectIndex_,
                                              omega, objects[objectIndex_], bground);
                                              
                                              
                                              
                                            
      Inc_local[p+pMax] = optimet::symbol::up_mn(p, C_10m1, C_11m1, nMax_,
                                              internalCoef_FF_,
                                              objectIndex_,
                                              omega, objects[objectIndex_], bground);            
       
                                              
                                              

      Inc_local[p + 2*pMax]=  std::complex<double> (0);//<- this last bit is zero for the moment, vpp
              
                                  
                  
      Inc_local[p + 3*pMax] = optimet::symbol::upp_mn(p, W_m1m1, W_00, W_11, W_10, W_01, nMax_, 
                  internalCoef_FF_,
                   objectIndex_,
                  omega, objects[objectIndex_]); 
                                                 
     
    }

  
  return 0;
}


// second harmonic sources adapted for paralellization

int Geometry::getIncLocalSH_parallel(std::vector<double *> CLGcoeff, int gran1, int gran2, std::shared_ptr<optimet::Excitation const> incWave_,
                       optimet::Vector<optimet::t_complex> &internalCoef_FF_, int nMaxS_, std::complex<double> *Inc_local) {
 
   double *C_10m1 = CLGcoeff[0];
   double *C_11m1 = CLGcoeff[1];
   double *C_00m1 = CLGcoeff[2];
   double *C_01m1 = CLGcoeff[3];
   
   double *W_m1m1 = CLGcoeff[4];
   double *W_11 = CLGcoeff[5];
   double *W_00 = CLGcoeff[6];
   double *W_10 = CLGcoeff[7];
   double *W_01 = CLGcoeff[8]; 
                             
  CompoundIterator p, q;
  
  int pMax = p.max(nMaxS_);
  int qMax = q.max(nMaxS_);
  int nMax_ = this->nMax();
  int objectIndex_;
  int brojac (0);
  int size = gran2 - gran1;
  

  auto const omega = incWave_->omega();

    
    for(q = gran1; q < gran2; q++) {  
     
      objectIndex_ = q / qMax;

      p = q - objectIndex_ * qMax;
    
       Inc_local[brojac] = optimet::symbol::vp_mn(p, C_00m1, C_01m1, nMax_,
                                              internalCoef_FF_,
                                              objectIndex_,
                                              omega, objects[objectIndex_], bground);                                        
                                              
                                           
      Inc_local[brojac+size] = optimet::symbol::up_mn(p, C_10m1, C_11m1, nMax_,
                                              internalCoef_FF_,
                                              objectIndex_,
                                              omega, objects[objectIndex_], bground);


       
                                              
      Inc_local[brojac + 2*size]=  std::complex<double> (0);//<- this last bit is zero for the moment, vpp
              
                                  
                  
      Inc_local[brojac + 3*size] = optimet::symbol::upp_mn(p, W_m1m1, W_00, W_11, W_10, W_01, nMax_, 
                  internalCoef_FF_,
                   objectIndex_,
                  omega, objects[objectIndex_]); 
                             
      brojac++;
     
    }

  
  return 0;
}


 int Geometry::AbsCSSHcoeff(std::vector<double *> CLGcoeff, int gran1, int gran2, std::shared_ptr<optimet::Excitation const> incWave_, 
optimet::Vector<optimet::t_complex> &internalCoef_FF_, optimet::Vector<optimet::t_complex> &internalCoef_SH_,
        int nMaxS_, std::complex<double> *coefABS) {
                 
                                

   double *W_m1m1 = CLGcoeff[4];
   double *W_11 = CLGcoeff[5];
   double *W_00 = CLGcoeff[6];
              
   CompoundIterator p, q;
     
   int pMax = p.max(nMaxS_);
   int qMax = q.max(nMaxS_);
   int nMax_ = this->nMax();

   int objectIndex_; 
   
  int NO = this->objects.size();  

  auto const omega = incWave_->omega();

  std::complex<double> cmnSH, dmnSH;
 
    int brojac(0);
   
    for(q = gran1; q < gran2; q++) {

      objectIndex_ = q / qMax;

      p = q - objectIndex_ * qMax;
  
  if (objects[0].scatterer_type == "sphere"){   
  
  cmnSH = internalCoef_SH_[objectIndex_ * 2 * pMax + p.compound];
   
  dmnSH = internalCoef_SH_[pMax + objectIndex_ * 2 * pMax + p.compound];
  }

  else if (objects[0].scatterer_type == "arbitrary.shape"){
 
 cmnSH = (1.0 / incWave_->waveK) * internalCoef_SH_[objectIndex_ * 2 * pMax + p.compound];

 dmnSH = (1.0 / incWave_->waveK) * internalCoef_SH_[objectIndex_ * 2 * pMax + pMax + p.compound];
 
 }
     
    
    coefABS[brojac] = optimet::symbol::ACSshcoeff(p, W_m1m1, W_00, W_11, nMax_, nMaxS_,
                  internalCoef_FF_, cmnSH, dmnSH,
                   objectIndex_,
                  omega, objects[objectIndex_]);

    brojac++;        
    }
   
  return 0;
}


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

void Geometry::getEXCvecSH_ARB3(optimet::Vector<optimet::t_complex>& EXvec, std::shared_ptr<optimet::Excitation const> excitation, optimet::Vector<optimet::t_complex> &externalCoef_FF_, int objIndex){
  using namespace optimet;
  
  int nMax = this->nMax();
  int nMaxS = this->nMaxS();
  auto const N = nMaxS * (nMaxS + 2);
  
          CompoundIterator mu1, mup, pp;
         
          optimet::t_real omega_ = excitation->omega();
          Cartesian<double> intpoinCar; // integration point in Cartesian system
          Spherical<double> intpoinSph; // integration point in spherical system

          SphericalP<std::complex<double>> Eext_FF;
          
          int muMax =mu1.max(nMaxS);
          int ppMax = pp.max(nMax);
          
          std::complex<double>* Etan = new std::complex<double>[3];
          std::complex<double>* Tan = new std::complex<double>[3];
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> Enorm, resDOT1, resDOT2, EE;
          
         const std::complex<double> k_0_SH = 2 * (omega_) * std::sqrt(consEpsilon0 * consMu0); 
         auto const k_s_SH = 2 * omega_ * std::sqrt(objects[objIndex].elmag.epsilon_SH * objects[objIndex].elmag.mu_SH); 
         auto const k_b_SH = 2 * omega_ * std::sqrt(bground.epsilon * bground.mu);
         auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu);
         auto const k_s = omega_ * std::sqrt(objects[objIndex].elmag.epsilon * objects[objIndex].elmag.mu); 
         
         const std::complex<double> ksiparppar = objects[objIndex].elmag.ksiparppar; // SH coefficient
         const std::complex<double> ksippp = objects[objIndex].elmag.ksippp; // SH coefficient
         const std::complex<double> gamma = objects[objIndex].elmag.gamma; // SH coefficient
         
         auto const Cf = (k_b_SH * consEpsilon0 * gamma)/objects[objIndex].elmag.epsilon_SH;
         unsigned int Nt = (objects[objIndex].getTopolsize()) / 3;  //number of triangles
        
        std::vector<std::vector<double>> Points = Tools::getPoints4();

	std::vector<double> Weights = Tools::getWeights4();

	for (mu1 = 0; mu1 < muMax ; mu1++){
	
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

          Eext_FF = SphericalP<std::complex<double>>(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
	
	double N1 = Points[ni][0];
	double N2 = Points[ni][1];
        double N0 = 1.0 - N1 - N2;

	double wi = Weights[ni];
         
	intpoinCar.x = (p1[0])*N1 + (p2[0])*N2 + (p3[0])*N0;
	intpoinCar.y = (p1[1])*N1 + (p2[1])*N2 + (p3[1])*N0;
	intpoinCar.z = (p1[2])*N1 + (p2[2])*N2 + (p3[2])*N0;
	
	intpoinSph = Tools::toSpherical(intpoinCar); 
	
	// Incoming field
	optimet::AuxCoefficients aCoefInc(intpoinSph, k_b, 1, nMax); //incident FF field
      
      optimet::AuxCoefficients aCoefSca(intpoinSph, k_b, 0, nMax); // scattered FF field
      
      // this works just for object centered at (0,0,0), so the global coordinates are equal to local!!!!!/////
      //  calculation of the total exterior field on the boundary of the particle (surface)
      for(pp = 0; pp < ppMax; pp++) {
     
        Eext_FF = Eext_FF + (aCoefInc.M(static_cast<long>(pp)) * excitation->dataIncAp[pp] +
                       aCoefInc.N(static_cast<long>(pp)) * excitation->dataIncBp[pp]); 
                       
        Eext_FF = Eext_FF + (aCoefSca.M(static_cast<long>(pp)) * externalCoef_FF_[pp] +
                       aCoefSca.N(static_cast<long>(pp)) * externalCoef_FF_[ppMax + pp]);    
                                                                  
      }
      
      // tan component -nxnxE
      Tools::crossTanTr(Etan, resCR, &nvec[0], Eext_FF);
      
      // inner component of normal field at FF
      Enorm = (bground.epsilon / objects[objIndex].elmag.epsilon) * Tools::dot (&nvec[0], Eext_FF);
      
     
      Tools::crossTan(Tan, resCR, &nvec[0], &Etan[0]);
      
      optimet::AuxCoefficients aCoefext3(intpoinSph, k_b_SH, 0, nMaxS); //radiative VSWF (3)

      resDOT1 = Tools::dot (&Tan[0], aCoefext3.M(static_cast<long>(mup)));
      resDOT2 = Tools::dot (&Tan[0], aCoefext3.N(static_cast<long>(mup)));
      
      I11s = I11s + wi * det * 2.0 * Enorm * ksiparppar * k_0_SH * k_0_SH * resDOT1;
      I12s = I12s + wi * det * 2.0 * Enorm * ksiparppar * k_0_SH * k_0_SH * resDOT2;
      
      resDOT1 = Tools::dot (&nvec[0], aCoefext3.M(static_cast<long>(mup)));
      resDOT2 = Tools::dot (&nvec[0], aCoefext3.N(static_cast<long>(mup)));
      
      I21s = I21s - wi * det * k_b_SH * k_b_SH *  ksippp * Enorm*Enorm * resDOT1;
      I22s = I22s - wi * det * k_b_SH * k_b_SH *  ksippp * Enorm*Enorm * resDOT2 ;
                                                                             
      
      // calculation of the particular solution on the surface // valid just for one object
      
      EE = Etan[0] * Etan[0] + Etan[1] * Etan[1] + Etan[2] * Etan[2] + Enorm*Enorm;
     resDOT1 = Tools::dot (&nvec[0], aCoefext3.M(static_cast<long>(mup)));
     resDOT2 = Tools::dot (&nvec[0], aCoefext3.N(static_cast<long>(mup)));

     I31s = I31s - ( wi * det * k_b_SH * Cf) * EE * resDOT1;
     I32s = I32s - ( wi * det * k_b_SH * Cf) * EE * resDOT2; 
  
  }// surface integration over the triangle is OVER
  } // all the triangles surface integration

  EXvec (mu1) = I11s + I21s + I31s;
  EXvec (mu1 + muMax) = I12s + I22s + I32s;
  
  }
  }

void Geometry::getEXCvecSH_ARB1(optimet::Vector<optimet::t_complex>& EXvec, std::shared_ptr<optimet::Excitation const> excitation, optimet::Vector<optimet::t_complex> &externalCoef_FF_, int objIndex){
  using namespace optimet;
  
  int nMax = this->nMax();
  int nMaxS = this->nMaxS();
  auto const N = nMaxS * (nMaxS + 2);
  
          CompoundIterator mu1, mup, pp;
         
          optimet::t_real omega_ = excitation->omega();
          Cartesian<double> intpoinCar; // integration point in Cartesian system
          Spherical<double> intpoinSph; // integration point in spherical system

          SphericalP<std::complex<double>> Eext_FF;
          
          int muMax =mu1.max(nMaxS);
          int ppMax = pp.max(nMax);
          
          std::complex<double>* Etan = new std::complex<double>[3];
          std::complex<double>* Tan = new std::complex<double>[3];
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> Enorm, resDOT1, resDOT2, EE;
          
         const std::complex<double> k_0_SH = 2 * (omega_) * std::sqrt(consEpsilon0 * consMu0); 
         auto const k_s_SH = 2 * omega_ * std::sqrt(objects[objIndex].elmag.epsilon_SH * objects[objIndex].elmag.mu_SH); 
         auto const k_b_SH = 2 * omega_ * std::sqrt(bground.epsilon * bground.mu);
         auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu);
         auto const k_s = omega_ * std::sqrt(objects[objIndex].elmag.epsilon * objects[objIndex].elmag.mu); 
         const std::complex<double> ksiparppar = objects[objIndex].elmag.ksiparppar; // SH coefficient
         const std::complex<double> ksippp = objects[objIndex].elmag.ksippp; // SH coefficient
         const std::complex<double> gamma = objects[objIndex].elmag.gamma; // SH coefficient
         
         auto const Cf = (k_b_SH * consEpsilon0 * gamma)/objects[objIndex].elmag.epsilon_SH;
          unsigned int Nt = (objects[objIndex].getTopolsize()) / 3;  //number of triangles
        
        std::vector<std::vector<double>> Points = Tools::getPoints4();

	std::vector<double> Weights = Tools::getWeights4();
	
	for (mu1 = 0; mu1 < muMax ; mu1++){
	
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
	
      Eext_FF = SphericalP<std::complex<double>>(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
	
	double N1 = Points[ni][0];
	double N2 = Points[ni][1];
        double N0 = 1.0 - N1 - N2;

	double wi = Weights[ni];
         
	intpoinCar.x = (p1[0])*N1 + (p2[0])*N2 + (p3[0])*N0;
	intpoinCar.y = (p1[1])*N1 + (p2[1])*N2 + (p3[1])*N0;
	intpoinCar.z = (p1[2])*N1 + (p2[2])*N2 + (p3[2])*N0;
	
	intpoinSph = Tools::toSpherical(intpoinCar); 
	
	// Incoming field
	optimet::AuxCoefficients aCoefInc(intpoinSph, k_b, 1, nMax); //incident FF field
      
      optimet::AuxCoefficients aCoefSca(intpoinSph, k_b, 0, nMax); // scattered FF field
      
      // this works just for object centered at (0,0,0), so the global coordinates are equal to local!!!!!/////
      // calculation of the total exterior field on the boundary of the particle (surface) 
      
      for(pp = 0; pp < ppMax; pp++) {
     
        Eext_FF = Eext_FF + (aCoefInc.M(static_cast<long>(pp)) * excitation->dataIncAp[pp] +
                       aCoefInc.N(static_cast<long>(pp)) * excitation->dataIncBp[pp]); 
                       
        Eext_FF = Eext_FF + (aCoefSca.M(static_cast<long>(pp)) * externalCoef_FF_[pp] +
                       aCoefSca.N(static_cast<long>(pp)) * externalCoef_FF_[ppMax + pp]);   
                                 
                                          
      }
      
      // -nxnxE
      Tools::crossTanTr(Etan, resCR, &nvec[0], Eext_FF);
 
      // inner component of normal field at FF
      Enorm = (bground.epsilon / objects[objIndex].elmag.epsilon) * Tools::dot (&nvec[0], Eext_FF);
      
      Tools::crossTan(Tan, resCR, &nvec[0], &Etan[0]);
      
      optimet::AuxCoefficients aCoefext1(intpoinSph, k_b_SH, 1, nMaxS); //incoming VSWF (1)
      
      resDOT1 = Tools::dot (&Tan[0], aCoefext1.M(static_cast<long>(mup)));
      resDOT2 = Tools::dot (&Tan[0], aCoefext1.N(static_cast<long>(mup)));
      
      I11s = I11s + wi * det * 2.0 * Enorm * ksiparppar * k_0_SH * k_0_SH * resDOT1;
      I12s = I12s + wi * det * 2.0 * Enorm * ksiparppar * k_0_SH * k_0_SH * resDOT2;
      
      resDOT1 = Tools::dot (&nvec[0], aCoefext1.M(static_cast<long>(mup)));
      resDOT2 = Tools::dot (&nvec[0], aCoefext1.N(static_cast<long>(mup)));
      
      I21s = I21s - wi * det * k_b_SH * k_b_SH * ksippp * Enorm*Enorm * resDOT1;
      I22s = I22s - wi * det * k_b_SH * k_b_SH * ksippp * Enorm*Enorm * resDOT2 ;
                                                                             
      
      // calculation of the particular solution on the surface // valid just for one object

      EE =  Etan[0] * Etan[0] + Etan[1] * Etan[1] + Etan[2] * Etan[2] + Enorm*Enorm;
     resDOT1 = Tools::dot (&nvec[0], aCoefext1.M(static_cast<long>(mup)));
     resDOT2 = Tools::dot (&nvec[0], aCoefext1.N(static_cast<long>(mup)));

     I31s = I31s - ( wi * det * k_b_SH * Cf) * EE * resDOT1;
     I32s = I32s - ( wi * det * k_b_SH * Cf) * EE * resDOT2;
  
  }// surface integration over the triangle is OVER
  
  } // all the triangles surface integration

  EXvec (mu1) = I11s + I21s + I31s;
  EXvec (mu1 + muMax) = I12s + I22s + I32s;
  
  }
  }

void Geometry::getEXCvecSH_ARB3_parall(optimet::Vector<optimet::t_complex>& EXvec, std::shared_ptr<optimet::Excitation const> excitation, optimet::Vector<optimet::t_complex> &externalCoef_FF_, int gran1, int gran2){
  using namespace optimet;
  
  int nMax = this->nMax();
  int nMaxS = this->nMaxS();
  auto const N = nMaxS * (nMaxS + 2);
  
          CompoundIterator mu1, mup, pp;
         
          optimet::t_real omega_ = excitation->omega();
          Cartesian<double> intpoinCar; // integration point in Cartesian system
          Spherical<double> intpoinSph; // integration point in spherical system
          
          SphericalP<std::complex<double>> Eext_FF;
          
          int muMax =mu1.max(nMaxS);
          int ppMax = pp.max(nMax);
          int objIndex(0), brojac(0);
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

         const std::complex<double> ksiparppar = objects[objIndex].elmag.ksiparppar; // SH coefficient
         const std::complex<double> ksippp = objects[objIndex].elmag.ksippp; // SH coefficient
         const std::complex<double> gamma = objects[objIndex].elmag.gamma; // SH coefficient   
   
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
	
      Eext_FF = SphericalP<std::complex<double>>(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
	
	double N1 = Points[ni][0];
	double N2 = Points[ni][1];
        double N0 = 1.0 - N1 - N2;

	double wi = Weights[ni];
         
	intpoinCar.x = (p1[0])*N1 + (p2[0])*N2 + (p3[0])*N0;
	intpoinCar.y = (p1[1])*N1 + (p2[1])*N2 + (p3[1])*N0;
	intpoinCar.z = (p1[2])*N1 + (p2[2])*N2 + (p3[2])*N0;
	
	intpoinSph = Tools::toSpherical(intpoinCar); 
	
	// Incoming field
	 optimet::AuxCoefficients aCoefInc(intpoinSph, k_b, 1, nMax); //incident FF field
      
      optimet::AuxCoefficients aCoefSca(intpoinSph, k_b, 0, nMax); // scattered FF field
      
      // this works just for object centered at (0,0,0), so the global coordinates are equal to local!!!!!/////
      // calculation of the total exterior field on the boundary of the particle (surface)
      for(pp = 0; pp < ppMax; pp++) {
     
        Eext_FF = Eext_FF + (aCoefInc.M(static_cast<long>(pp)) * excitation->dataIncAp[pp] +
                       aCoefInc.N(static_cast<long>(pp)) * excitation->dataIncBp[pp]); 
                       
        Eext_FF = Eext_FF + (aCoefSca.M(static_cast<long>(pp)) * externalCoef_FF_[pp] +
                       aCoefSca.N(static_cast<long>(pp)) * externalCoef_FF_[ppMax + pp]);    
                                                                  
      }
      // -nxnxE
      Tools::crossTanTr(Etan, resCR, &nvec[0], Eext_FF);
      
      // inner component of normal field at FF
      Enorm = (bground.epsilon / objects[objIndex].elmag.epsilon) * Tools::dot (&nvec[0], Eext_FF);
                        
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
                                                                             
      
      // calculation of the particular solution on the surface // valid just for one object
             
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
  }
  }

  void Geometry::getEXCvecSH_ARB1_parall(optimet::Vector<optimet::t_complex>& EXvec, std::shared_ptr<optimet::Excitation const> excitation, optimet::Vector<optimet::t_complex> &externalCoef_FF_, int gran1, int gran2){
  using namespace optimet;
  
  int nMax = this->nMax();
  int nMaxS = this->nMaxS();
  auto const N = nMaxS * (nMaxS + 2);
  
          CompoundIterator mu1, mup, pp;
         
          optimet::t_real omega_ = excitation->omega();
          Cartesian<double> intpoinCar; // integration point in Cartesian system
          Spherical<double> intpoinSph; // integration point in spherical system
          
          SphericalP<std::complex<double>> Eext_FF;
          
          int muMax =mu1.max(nMaxS);
          int ppMax = pp.max(nMax);
          int objIndex(0), brojac(0);
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
         const std::complex<double> ksiparppar = objects[objIndex].elmag.ksiparppar; // SH coefficient
         const std::complex<double> ksippp = objects[objIndex].elmag.ksippp; // SH coefficient
         const std::complex<double> gamma = objects[objIndex].elmag.gamma; // SH coefficient
        
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
	
      Eext_FF = SphericalP<std::complex<double>>(std::complex<double>(0.0, 0.0), std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0));
	
	double N1 = Points[ni][0];
	double N2 = Points[ni][1];
        double N0 = 1.0 - N1 - N2;

	double wi = Weights[ni];
         
	intpoinCar.x = (p1[0])*N1 + (p2[0])*N2 + (p3[0])*N0;
	intpoinCar.y = (p1[1])*N1 + (p2[1])*N2 + (p3[1])*N0;
	intpoinCar.z = (p1[2])*N1 + (p2[2])*N2 + (p3[2])*N0;
	
	intpoinSph = Tools::toSpherical(intpoinCar); 
	
	// Incoming field
	optimet::AuxCoefficients aCoefInc(intpoinSph, k_b, 1, nMax); //incident FF field
      
      optimet::AuxCoefficients aCoefSca(intpoinSph, k_b, 0, nMax); // scattered FF field

      // this works just for object centered at (0,0,0), so the global coordinates are equal to local!!!!!/////
      // calculation of the total exterior field on the boundary of the particle (surface)
      for(pp = 0; pp < ppMax; pp++) {
     
        Eext_FF = Eext_FF + (aCoefInc.M(static_cast<long>(pp)) * excitation->dataIncAp[pp] +
                       aCoefInc.N(static_cast<long>(pp)) * excitation->dataIncBp[pp]); 
                       
        Eext_FF = Eext_FF + (aCoefSca.M(static_cast<long>(pp)) * externalCoef_FF_[pp] +
                       aCoefSca.N(static_cast<long>(pp)) * externalCoef_FF_[ppMax + pp]);   
                                 
                                          
      }
      
      //-nxnxE
      Tools::crossTanTr(Etan, resCR, &nvec[0], Eext_FF);
 
      // inner component of normal field at FF
      Enorm = (bground.epsilon / objects[objIndex].elmag.epsilon) * Tools::dot (&nvec[0], Eext_FF);
      
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
                                                                             
      
      // calculation of the particular solution on the surface // valid just for one object

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

void Geometry::update(std::shared_ptr<optimet::Excitation const> incWave_) {
  // Update the ElectroMagnetic properties of each object
  for(auto &object : objects)
    object.elmag.update(incWave_->lambda());
}

