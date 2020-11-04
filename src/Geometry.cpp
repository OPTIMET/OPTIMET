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

#include "Coupling.h"
#include "Geometry.h"

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



Spherical<double> Geometry::translateCoordinates(int firstObject_, int secondObject_) {
  return objects[firstObject_].vR - objects[secondObject_].vR;
}

int Geometry::checkInner(Spherical<double> R_) {
  Spherical<double> Rrel;

  for(size_t j = 0; j < objects.size(); j++) {
    // Translate to object j
    Rrel = Tools::toPoint(R_, objects[j].vR);

    // Check if R_.rrr is inside radius of object j
    if(Rrel.rrr <= objects[j].radius)
      return j;
  }

  return -1;
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
    
  cmnSH = internalCoef_SH_[objectIndex_ * 4 * pMax + p.compound] + internalCoef_SH_[2 * pMax + objectIndex_ * 4 * pMax + p.compound];
   
  dmnSH = internalCoef_SH_[pMax + objectIndex_ * 4 * pMax + p.compound] + internalCoef_SH_[3 * pMax + 
                           objectIndex_ * 4 * pMax + p.compound];
     
    
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


optimet::Vector<optimet::t_complex> Geometry::getTLocal(optimet::t_real omega_,
                                                        optimet::t_int objectIndex_,
                                                        optimet::t_uint nMax_) const {
  if(objectIndex_ >= static_cast<int>(objects.size()))
    std::out_of_range("Object index out of range");
  if(nMax_ != static_cast<optimet::t_uint>(objects[objectIndex_].nMax))
    std::runtime_error("Internal inconsistency");
  return objects[objectIndex_].getTLocal(omega_, bground);
}



void Geometry::update(std::shared_ptr<optimet::Excitation const> incWave_) {
  // Update the ElectroMagnetic properties of each object
  for(auto &object : objects)
    object.elmag.update(incWave_->lambda());
}

void Geometry::updateRadius(double radius_, int object_) { objects[object_].radius = radius_; }

void Geometry::rebuildStructure() {
  if(structureType == 1) {
    // Spiral structure needs to be rebuilt
    double R, d;
    int Np, No;
    double Theta;

    Np = ((objects.size() - 1) / 4) + 1;

    Theta = consPi / (Np - 1); // Calculate the separation angle

    d = 2 * (spiralSeparation + 2 * objects[0].radius);
    R = d / (4 * sin(Theta / 2.0));

    No = objects.size();

    // Create vectors for r, theta, x and y, X and Y
    std::vector<double> X(No - 1);
    std::vector<double> Y(No - 1);

    int j = 0;

    //"Vertical" arms
    for(int i = 0; i < Np; i++) {
      double x_;
      double y_;

      Tools::Pol2Cart(R, i * Theta, x_, y_);

      if(i == 0) {
        X[j] = x_ + R;
        Y[j] = y_;
        j++;
        continue;
      }

      if(i == Np - 1) {
        X[j] = x_ - R;
        Y[j] = -1.0 * y_;
        j++;
        continue;
      }

      X[j] = x_ + R;
      Y[j] = y_;
      j++;

      X[j] = x_ - R;
      Y[j] = -1.0 * y_;
      j++;
    }

    //"Horizontal" arms
    for(int i = 0; i < Np; i++) {
      double x_;
      double y_;
      Tools::Pol2Cart(R, i * Theta - consPi / 2, x_, y_);

      if(i == 0) {
        X[j] = x_;
        Y[j] = y_ - R;
        j++;
        continue;
      }

      if(i == Np - 1) {
        X[j] = -1.0 * x_;
        Y[j] = y_ + R;
        j++;
        continue;
      }

      X[j] = -1.0 * x_;
      Y[j] = y_ + R;
      j++;

      X[j] = x_;
      Y[j] = y_ - R;
      j++;
    }

    // Determine normal, convert to a spherical object and push

    Cartesian<double> auxCar(0.0, 0.0, 0.0);
    Spherical<double> auxSph(0.0, 0.0, 0.0);

    for(int i = 0; i < No - 1; i++) {
      if(normalToSpiral == 0) {
        // x is normal (conversion is x(pol) -> y; y(pol) -> z
        auxCar.x = 0.0;
        auxCar.y = X[i];
        auxCar.z = Y[i];
      }

      if(normalToSpiral == 1) {
        // y is normal (conversion is x(pol) -> z; y(pol) -> x
        auxCar.x = Y[i];
        auxCar.y = 0.0;
        auxCar.z = X[i];
      }

      if(normalToSpiral == 2) {
        // z is normal (conversion is x(pol) -> x; y(pol) -> x
        auxCar.x = X[i];
        auxCar.y = Y[i];
        auxCar.z = 0.0;
      }

      auxSph = Tools::toSpherical(auxCar);
      objects[i + 1].vR = auxSph;
    }
  }
}
