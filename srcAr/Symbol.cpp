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

#include "Symbol.h"

#include <cmath>
#include "constants.h"
#include "Bessel.h"
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


// A numbers
std::complex<double> A_0(int n, double R, const std::complex<double> &waveK_i,
                         const std::complex<double> &cmn, int nMax) {
                                                 
  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = bessel<Bessel, false>(R * waveK_i, nMax);  

  return data[n] * cmn;
}

std::complex<double> A_1(int n, double R, const std::complex<double> &waveK_i,
                         const std::complex<double> &dmn, int nMax) {
  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = bessel<Bessel, false>(R * waveK_i, nMax); 

  return std::complex<double>(0.0, 1.0) * (1.0 / waveK_i) *
         (waveK_i * ddata[n] + data[n] / R) * dmn;   
}

std::complex<double> A_m1(int n, double R, const std::complex<double> &waveK_i,
                          const std::complex<double> &dmn, int nMax) {
  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = bessel<Bessel, false>(R * waveK_i, nMax);


  return std::complex<double>(0.0, 1.0) * std::sqrt(n * (n + 1.0)) *
         (1.0 / waveK_i / R) * data[n] * dmn;  
}

std::complex<double> F_00(int n1, int n2, double r, const std::complex<double> &waveK_i,
                          int nMax) {
 
  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = bessel<Bessel, false>(r * waveK_i, nMax);


  return  data[n1] * ddata[n2];
}

std::complex<double> F_11(int n1, int n2, double r, const std::complex<double> &waveK_i,
                          int nMax) {
 
  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = bessel<Bessel, false>(r * waveK_i, nMax);

  return (waveK_i * ddata[n2] + data[n2]/r) * (waveK_i * ddata[n1] + data[n1]/r);
}

std::complex<double> F_m1m1(int n1, int n2, double r, const std::complex<double> &waveK_i,
                          int nMax) {
 
  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = bessel<Bessel, false>(r * waveK_i, nMax);


  return  (1.0 / (std::pow(r , 2.0))) * (data[n1] * data[n2]);
}

std::complex<double> F_d00(int n1, int n2, double r, const std::complex<double> &waveK_i,
                          int nMax) {
 
  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = bessel<Bessel, false>(r * waveK_i, nMax);


  return  data[n1] * ddata[n2] * waveK_i + waveK_i * ddata[n1] * data[n2];
}

std::complex<double> F_d11(int n1, int n2, double r, const std::complex<double> &waveK_i,
                          int nMax) {
 
  std::vector<std::complex<double>> data, ddata, dddata;
  std::tie(data, ddata) = bessel<Bessel, false>(r * waveK_i, nMax);
  dddata = bessel3der<Bessel, false>(r * waveK_i, nMax);

  return  (std::pow(waveK_i , 2.0) * dddata[n1] - (1.0 / std::pow(r , 2.0)) * data[n1] + (1.0 / r) * waveK_i * ddata[n1]) *
          (waveK_i * ddata[n2] + data[n2]/r) + (std::pow(waveK_i , 2.0) * dddata[n2] - (1.0 / std::pow(r , 2.0)) * data[n2] 
+ (1.0/r) * waveK_i * ddata[n2]) *
          (waveK_i * ddata[n1] + data[n1]/r);
}

std::complex<double> F_dm1m1(int n1, int n2, double r, const std::complex<double> &waveK_i,
                          int nMax) {
 
  std::vector<std::complex<double>> data, ddata;
  std::tie(data, ddata) = bessel<Bessel, false>(r * waveK_i, nMax);


  return  (1.0 / (std::pow(r , 2.0))) * (waveK_i * ddata[n1] * data[n2] + data[n1] * waveK_i * ddata[n2]) - (2.0 / 
(std::pow(r , 3.0))) * (data[n1] * data[n2]);
}

double W(int L1, int J1, int M1, int L2, int J2, int M2, int L, int M) {

  return std::pow(-1.0, J2 + L1 + L) *
         std::sqrt((2.0 * J1 + 1.0) * (2.0 * J2 + 1.0) * (2.0 * L1 + 1.0) *
                   (2.0 * L2 + 1.0) / (4.0 * consPi * (2.0 * L + 1.0))) *
         Wigner6j(L1, L2, L, J2, J1, 1.0) * CleGor(L, 0, L1, 0, L2, 0) *
         CleGor(L, M, J1, M1, J2, M2);  
}

} // namespace

                          
// function for coefficients with Xm1 spherical function, particular solution of diff equations, SH
std::complex<double> CXm1(CompoundIterator &kk, double *W_m1m1, double *W_00, double *W_11, int nMax,
                            optimet::Vector<optimet::t_complex> &internalCoef_FF_,
                            double r,
                            int objectIndex_, double omega,
                            const Scatterer &object) {
                            
   int n = kk.first;
   int m = kk.second;                         

  // Basic relations
  const double R = object.radius;
  const std::complex<double> eps_0 = consEpsilon0;
  const std::complex<double> mu_0 = consMu0;
  const std::complex<double> mu_j = object.elmag.mu;
  const std::complex<double> eps_j = object.elmag.epsilon;
  
  
  // SH Basic relations

  const std::complex<double> gamma = object.elmag.gamma;
  const std::complex<double> eps_j2 = object.elmag.epsilon_SH;
  const std::complex<double> mu_j2 = object.elmag.mu_SH;
  
  // Auxiliary variables
  const std::complex<double> waveK_01 = (omega) * std::sqrt(eps_0 * mu_0); 
  const std::complex<double> waveK_j1 = (omega) * std::sqrt(eps_j * mu_j);


  double W00, W11, Wm1m1;
  
  CompoundIterator p, q;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);
  
  int size1 = pMax * qMax;
 
 std::complex<double> cmn_1, dmn_1, cmn_2, dmn_2; 
 
  
 std::complex<double> COEFFXm1(0.0, 0.0);
 
 int brojac (0);
 
 
 for(p = 0; p <  pMax; p++) {  
 
 cmn_1 = internalCoef_FF_(objectIndex_ * 2 * pMax + p.compound);
 
  dmn_1 = internalCoef_FF_(pMax + objectIndex_ * 2 * pMax + p.compound);
  
  for(q = 0; q <  qMax; q++) {    
  
  cmn_2 = internalCoef_FF_(objectIndex_ * 2 * qMax + q.compound);
 
   dmn_2 = internalCoef_FF_(qMax + objectIndex_ * 2 * qMax + q.compound);
   
   Wm1m1 = W_m1m1[kk*size1 +brojac];
                    
   W00 = W_00[kk*size1 +brojac];                 
   
   W11 =  W_11[kk*size1 +brojac];
                 
   
   COEFFXm1 += (-eps_0 / eps_j2) * gamma * ( cmn_1 * cmn_2 * W00 * F_d00(p.first, q.first, r, waveK_j1, nMax)

      + dmn_1 * dmn_2 * (1.0 / (std::pow(waveK_j1, 2.0))) * (W11 * F_d11(p.first, q.first, r, waveK_j1, nMax)
       
       + Wm1m1 * std::sqrt(p.first * q.first * (p.first + 1) * (q.first + 1)) *  F_dm1m1(p.first, q.first, r, waveK_j1, nMax)));
                          
                        
    brojac ++;          
   }
 
  }
        
  return COEFFXm1 ;    
}


// function for coefficients with Xp1 spherical function, particular solution of diff equations, SH
std::complex<double> CXp1(CompoundIterator &kk, double *W_m1m1, double *W_00, double *W_11, int nMax,
                            optimet::Vector<optimet::t_complex> &internalCoef_FF_,
                            double r,
                            int objectIndex_, double omega,
                            const Scatterer &object) {
                            
  int n = kk.first;
  int m = kk.second;                          

  // Basic relations
  const double R = object.radius;
  const std::complex<double> eps_0 = consEpsilon0;
  const std::complex<double> mu_0 = consMu0;
  const std::complex<double> mu_j = object.elmag.mu;
  const std::complex<double> eps_j = object.elmag.epsilon;
  
  
  // SH Basic relations

  const std::complex<double> gamma = object.elmag.gamma;
  const std::complex<double> eps_j2 = object.elmag.epsilon_SH;
  const std::complex<double> mu_j2 = object.elmag.mu_SH;
  
  // Auxiliary variables
  const std::complex<double> waveK_01 = (omega) * std::sqrt(eps_0 * mu_0); 
  const std::complex<double> waveK_j1 = (omega) * std::sqrt(eps_j * mu_j);

  double  W00, W11, Wm1m1;
  
  CompoundIterator p, q;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);
  
  int size1 = pMax * qMax;
 
 std::complex<double> cmn_1, dmn_1, cmn_2, dmn_2; 
 

 std::complex<double> COEFFXp1(0.0, 0.0);
 
 int brojac (0);
 
 for(p = 0; p <  pMax; p++) {  
 
 cmn_1 = internalCoef_FF_(objectIndex_ * 2 * pMax + p.compound);
 
  dmn_1 = internalCoef_FF_(pMax + objectIndex_ * 2 * pMax + p.compound);
  
  for(q = 0; q <  qMax; q++) {    
  
  cmn_2 = internalCoef_FF_(objectIndex_ * 2 * qMax + q.compound);
 
   dmn_2 = internalCoef_FF_(qMax + objectIndex_ * 2 * qMax + q.compound);
   
   Wm1m1 = W_m1m1[kk*size1 +brojac];
                    
   W00 = W_00[kk*size1 +brojac];                 
   
   W11 = W_11[kk*size1 +brojac];
                 
   
                          
   
   COEFFXp1 += (-eps_0 / eps_j2) * gamma * ( cmn_1 * cmn_2 * W00 * std::sqrt(n * (n + 1)) * (1.0 / r) * F_00(p.first, q.first, r, waveK_j1, nMax) + 
              dmn_1  * dmn_2 * (1.0 / (std::pow(waveK_j1, 2.0))) * (W11 * std::sqrt(n * (n + 1)) * (1.0 / r) * 

               F_11(p.first, q.first, r, waveK_j1, nMax) + Wm1m1 * std::sqrt(p.first * q.first * (p.first + 1) * (q.first + 1)) * 
              
               std::sqrt(n * (n + 1)) * (1.0 / r) * F_m1m1(p.first, q.first, r, waveK_j1, nMax)));      
               
               brojac++;                  
              
   }
 
  }
  
   
  return COEFFXp1 ;    

}

 // C numbers
#ifdef OPTIMET_MPI
void C_10m1coeff (double *C_10m1, int nMax, int nMaxS, int gran1, int gran2){

CompoundIterator p, q, k;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);
  int kMax = k.max(nMaxS);
  
 int J, M, J1, M1, J2, M2, brojac(0);
  
for (int tt = gran1; tt < gran2; tt++) {
 
 q = tt % qMax;
 p = (tt / qMax) % pMax;
 k = tt / (pMax * qMax);


J =k.first;
M = k.second;

J1 =p.first;
M1 = p.second;

J2 =q.first;
M2 = q.second; 



C_10m1[brojac] =std::sqrt(3.0 / 2.0 / consPi) * (2.0 * J1 + 1.0) *
         CleGor(J, M, J1, M1, J2, M2) *
         (std::sqrt(J2 * (2.0 * J2 - 1.0)) *  
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
          std::sqrt((J2 + 1) * (2.0 * J2 + 3.0)) *      
              Wigner9j(J1, J1, 1, J2, J2 + 1, 1, J, J - 1, 1) *
              CleGor(J - 1, 0, J1, 0, J2 + 1, 0) *
              std::sqrt((J + 1) / (2.0 * J + 1.0)));
              
     brojac++;           
    }

}


void C_11m1coeff (double *C_11m1, int nMax, int nMaxS, int gran1, int gran2){

CompoundIterator p, q, k;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);
  int kMax = k.max(nMaxS);
  
  int J, M, J1, M1, J2, M2, brojac(0);
  
for (int tt = gran1; tt < gran2; tt++) {

 q = tt % qMax;
 p = (tt / qMax) % pMax;
 k = tt / (pMax * qMax);

J =k.first;
M = k.second;

J1 =p.first;
M1 = p.second;

J2 =q.first;
M2 = q.second;

  C_11m1[brojac] = std::sqrt(3.0 / 2.0 / consPi) * CleGor(J, M, J1, M1, J2, M2) *
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
              Wigner9j(J1, J1 + 1, 1, J2, J2 + 1, 1, J, J + 1, 1) *         // seems alright
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
              
 brojac++;
}
}

void C_00m1coeff(double *C_00m1, int nMax, int nMaxS, int gran1, int gran2){

CompoundIterator p, q, k;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);
  int kMax = k.max(nMaxS);
  
  int J, M, J1, M1, J2, M2, brojac(0);
  
for (int tt = gran1; tt < gran2; tt++) {

 q = tt % qMax;
 p = (tt / qMax) % pMax;
 k = tt / (pMax * qMax);

J =k.first;
M = k.second;

J1 =p.first;
M1 = p.second;

J2 =q.first;
M2 = q.second;

  C_00m1[brojac] = std::sqrt(3.0 / 2.0 / consPi) * (2.0 * J1 + 1.0) *
         CleGor(J, M, J1, M1, J2, M2) *
         (std::sqrt(J2 * (2.0 * J2 - 1.0)) * 
              Wigner9j(J1, J1, 1, J2, J2 - 1, 1, J, J, 1) *
              CleGor(J, 0, J1, 0, J2 - 1, 0) -
          std::sqrt((J2 + 1.0) * (2.0 * J2 + 3.0)) *
              Wigner9j(J1, J1, 1, J2, J2 + 1, 1, J, J, 1) *
              CleGor(J, 0, J1, 0, J2 + 1, 0));
              
              brojac++;
}
}

void C_01m1coeff(double *C_01m1, int nMax, int nMaxS, int gran1, int gran2){

CompoundIterator p, q, k;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);
  int kMax = k.max(nMaxS);
  
  int J, M, J1, M1, J2, M2, brojac(0);
  
for (int tt = gran1; tt < gran2; tt++) {

 q = tt % qMax;
 p = (tt / qMax) % pMax;
 k = tt / (pMax * qMax);

J =k.first;
M = k.second;

J1 =p.first;
M1 = p.second;

J2 =q.first;
M2 = q.second;

  C_01m1[brojac] = std::sqrt(3.0 / 2.0 / consPi) * CleGor(J, M, J1, M1, J2, M2) *
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
          std::sqrt(J1 * (J2 + 1) * (2.0 * J1 + 3.0) * (2.0 * J2 + 3.0)) * 
              Wigner9j(J1, J1 + 1.0, 1.0, J2, J2 + 1.0, 1.0, J, J, 1.0) *
              CleGor(J, 0, J1 + 1.0, 0, J2 + 1.0, 0));
              
             brojac++;
}
}

// W numbers

void W_m1m1coeff (double *W_m1m1, int nMax, int nMaxS, int gran1, int gran2){

CompoundIterator p, q, k;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);
  int kMax = k.max(nMaxS);
  
  int J, M, J1, M1, J2, M2, brojac(0);
  
for (int tt = gran1; tt < gran2; tt++) {

 q = tt % qMax;
 p = (tt / qMax) % pMax;
 k = tt / (pMax * qMax);

J =k.first;
M = k.second;

J1 =p.first;
M1 = p.second;

J2 =q.first;
M2 = q.second;

W_m1m1[brojac] = std::sqrt(J1 / (2.0 * J1 + 1.0)) *
                  std::sqrt(J2 / (2.0 * J2 + 1.0)) * 
                  W(J1 - 1, J1, M1, J2 - 1, J2,
                    M2, J, M) +
              std::sqrt((J1 + 1) / (2.0 * J1 + 1.0)) *
                  std::sqrt((J2 + 1) / (2.0 * J2 + 1.0)) *
                  W(J1 + 1, J1, M1, J2 + 1, J2,  
                    M2, J, M) -
              std::sqrt(J1 / (2.0 * J1 + 1.0)) *
                  std::sqrt((J2 + 1) / (2.0 * J2 + 1.0)) *
                  W(J1 - 1, J1, M1, J2 + 1, J2,
                    M2, J, M) -
              std::sqrt((J1 + 1) / (2.0 * J1 + 1.0)) *
                  std::sqrt(J2 / (2.0 * J2 + 1.0)) *
                  W(J1 + 1, J1, M1, J2 - 1, J2,
                    M2, J, M);
 
           brojac++;
              
       }
              
  }
     
     
void W_11coeff (double *W_11, int nMax, int nMaxS, int gran1, int gran2){

CompoundIterator p, q, k;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);
  int kMax = k.max(nMaxS);
  
  int J, M, J1, M1, J2, M2, brojac(0);
  
for (int tt = gran1; tt < gran2; tt++) {

 q = tt % qMax;
 p = (tt / qMax) % pMax;
 k = tt / (pMax * qMax);

J =k.first;
M = k.second;

J1 =p.first;
M1 = p.second;

J2 =q.first;
M2 = q.second;

  W_11[brojac] =  std::sqrt((J1 + 1.0) / (2.0 * J1 + 1.0)) *
               std::sqrt((J2 + 1.0) / (2.0 * J2 + 1.0)) *
               W(J1 - 1, J1, M1, J2 - 1, J2, M2,
                 J, M) +
           std::sqrt(J1 / (2.0 * J1 + 1.0)) *
               std::sqrt(J2 / (2.0 * J2 + 1.0)) *
               W(J1 + 1.0, J1, M1, J2 + 1.0, J2, M2,
                 J, M) +                                                             
           std::sqrt((J1 + 1.0) / (2.0 * J1 + 1.0)) *
               std::sqrt(J2 / (2.0 * J2 + 1.0)) *
               W(J1 - 1.0, J1, M1, J2 + 1.0, J2, M2,
                 J, M) +
           std::sqrt(J1 / (2.0 * J1 + 1.0)) *
               std::sqrt((J2 + 1.0) / (2.0 * J2 + 1.0)) *
               W(J1 + 1.0, J1, M1, J2 - 1.0, J2, M2,
                 J, M);

         brojac++;  
                 
        }
  }
  
void W_00coeff (double *W_00, int nMax, int nMaxS, int gran1, int gran2){

CompoundIterator p, q, k;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);
  int kMax = k.max(nMaxS);
  
  int J, M, J1, M1, J2, M2, brojac(0);

for (int tt = gran1; tt < gran2; tt++) {

 q = tt % qMax;
 p = (tt / qMax) % pMax;
 k = tt / (pMax * qMax);

J =k.first;
M = k.second;

J1 =p.first;
M1 = p.second;

J2 =q.first;
M2 = q.second;

W_00[brojac] = W(J1, J1, M1, J2, J2, M2, J, M);

brojac++;

}
}

void W_10coeff (double *W_10, int nMax, int nMaxS, int gran1, int gran2){

CompoundIterator p, q, k;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);
  int kMax = k.max(nMaxS);
  
  int J, M, J1, M1, J2, M2, brojac(0);

for (int tt = gran1; tt < gran2; tt++) {

 q = tt % qMax;
 p = (tt / qMax) % pMax;
 k = tt / (pMax * qMax);

J =k.first;
M = k.second;

J1 =p.first;
M1 = p.second;

J2 =q.first;
M2 = q.second;

 W_10[brojac] = std::sqrt((J1 + 1.0) / (2.0 * J1 + 1.0)) *
                    W(J1 - 1.0, J1, M1, J2, J2,
                      M2, J, M) +
                std::sqrt(J1 / (2.0 * J1 + 1.0)) *
                    W(J1 + 1.0, J1, M1, J2, J2,
                      M2, J, M);

    brojac++;

   }
  }          


void W_01coeff (double *W_01, int nMax, int nMaxS, int gran1, int gran2){

CompoundIterator p, q, k;
  int pMax = p.max(nMax);
  int qMax = q.max(nMax);
  int kMax = k.max(nMaxS);
  
  int J, M, J1, M1, J2, M2, brojac(0);
  
for (int tt = gran1; tt < gran2; tt++) {

 q = tt % qMax;
 p = (tt / qMax) % pMax;
 k = tt / (pMax * qMax);

J =k.first;
M = k.second;

J1 =p.first;
M1 = p.second;

J2 =q.first;
M2 = q.second;

   W_01[brojac] = std::sqrt((J2 + 1.0) / (2.0 * J2 + 1.0)) *
                    W(J1, J1, M1, J2 - 1.0, J2,
                      M2, J, M) +
                  std::sqrt(J2 / (2.0 * J2 + 1.0)) *
                    W(J1, J1, M1, J2 + 1.0, J2,
                      M2, J, M);
  
         brojac++;
                      
     }
   }                   
#endif
         
} // namespace symbol
} // namespace optimet
