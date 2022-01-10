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

#include "ElectroMagnetic.h"
#include "constants.h"

ElectroMagnetic::ElectroMagnetic() {
  init_r(std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0), std::complex<double>(1.0, 0.0));
}

ElectroMagnetic::ElectroMagnetic(std::complex<double> epsilon_r_,
                                 std::complex<double> mu_r_, std::complex<double> epsilon_r_SH, std::complex<double> ksippp, std::complex<double> ksiparppar, std::complex<double> gamma) {
  init_r(epsilon_r_, mu_r_, epsilon_r_SH, ksippp, ksiparppar, gamma);
}

ElectroMagnetic::~ElectroMagnetic() {
  //
}


void ElectroMagnetic::init_r(std::complex<double> epsilon_r_,
                             std::complex<double> mu_r_, std::complex<double> epsilon_r_SH_, std::complex<double> ksippp_, std::complex<double> ksiparppar_, std::complex<double> gamma_) {
  epsilon_r = epsilon_r_;
  mu_r = mu_r_;

  epsilon = epsilon_r * consEpsilon0;
  mu = mu_r * consMu0;
  
  epsilon_r_SH = epsilon_r_SH_;
  mu_r_SH = mu_r_;

  epsilon_SH = epsilon_r_SH * consEpsilon0;
  mu_SH = mu_r_SH * consMu0;
  
  ksippp = ksippp_;
  ksiparppar = ksiparppar_;
  gamma = gamma_;


  modelType = 0;
}


void ElectroMagnetic::initHydrodynamicModel_r(std::complex<double> a_,
                                       std::complex<double> b_,
                                       std::complex<double> d_,
                                       std::complex<double> mu_r_) {
  a_SH = a_;
  b_SH = -b_;
  d_SH = d_;

  mu_r = mu_r_;

  modelType = 3;
}



void ElectroMagnetic::populateHydrodynamicModel() {

  // interpolation of Johnson_Christy
  double input_freq = consC / lambda;
  
  double input_omega = 2*consPi*input_freq;
  
  double a0[5], a1[5], b0[5], b1[5], b2[5];
  
  a0[0]=2.000003399882560;
  a0[1]=1.782388034422510e+32;
  a0[2]=9.571140818411450e+26;
  a0[3]=3.141025290600320e+24;
  a0[4]=5.056282927859510e+31;
  
  a1[0]=0.0;
  a1[1]=0.0;
  a1[2]=8.034165109695690e+15;
  a1[3]=1.060027902520820e+14;
  a1[4]=2.176317566053640e+16;
  
  b0[0]=1.0;
  b0[1]=0.0;
  b0[2]=1.398566311205070e+26;
  b0[3]=5.984581206741880e+23;
  b0[4]=1.707510287416960e+31;
  
  b1[0]=1.326291192399820e-15;
  b1[1]=1.122727361975370e+14;
  b1[2]=7.280057361739550e+15;
  b1[3]=4.393809682455200e+15;
  b1[4]=3.256258123271410e+15;
  
  b2[0]=0.0;
  b2[1]=1.0;
  b2[2]=1.0;
  b2[3]=1.0;
  b2[4]=1.0;
  
  std::complex<double> sumFF, sumSH;
  
  for(int i = 0; i <  5; i++) {
  
  sumFF += (a0[i] + input_omega * std::complex<double>(0.0, -1.0)*a1[i]) / (b0[i] + input_omega * std::complex<double>(0.0, -1.0)*b1[i] + std::pow(input_omega * std::complex<double>(0.0, -1.0), 2) *b2[i]); 
  
  
  sumSH += (a0[i] + 2.0*input_omega * std::complex<double>(0.0, -1.0)*a1[i]) / (b0[i] + 2.0*input_omega * std::complex<double>(0.0, -1.0)*b1[i] + std::pow(2.0*input_omega * std::complex<double>(0.0, -1.0), 2)     *b2[i]);
  
  }
  
  
  epsilon_r = 1. + sumFF;
                    
  epsilon_r_SH = 1. + sumSH;

  epsilon_SH = epsilon_r_SH * consEpsilon0;
  
  epsilon = epsilon_r * consEpsilon0;
  
 
 
  double mele = 9.10938356e-31 ;
  
  double charge = 1.602176e-19 ;
  
  
  ksippp = - (a_SH / 4.0) * (epsilon_r - 1.0) * (charge) / (mele * pow(2.0 * consPi * input_freq, 2.0));
  ksiparppar = - (b_SH / 2.0) * (epsilon_r - 1.0) * (charge) / (mele * pow(2.0 * consPi * input_freq, 2.0));
  gamma = - (d_SH / 8.0) * (epsilon_r - 1.0) * (charge) / (mele * pow(2.0 * consPi * input_freq, 2.0));
  
}


void ElectroMagnetic::initSiliconModel_r(std::complex<double> mu_r_) {
  
  refInd = {1.6370, 1.7370, 2.0300, 2.8400,4.1850,5.0490,5.0910,5.0850,5.1350,5.2450,5.4230,5.9140,6.8200,6.5870,6.0250,
5.6230,5.3410,5.1100,4.9320,4.7900,4.6730,4.5720,4.4850,4.4120,4.3490,4.2890,4.2350,
4.1870,4.1450,4.1030,4.0730,4.0380,4.0060,3.9770,3.9540,3.9310,3.9080,3.8880,3.8690,
3.8510,3.8350,3.8170,3.8050,3.7910,3.7760,3.7650,3.7530,3.7410,3.7300,3.7190,3.7120,
3.7010,3.6930,3.6840,3.6770,3.6690,3.6620,3.6550,3.6460,3.6410,3.6360,3.6280,3.6220,
3.6170,3.6130,3.6100,3.6040,3.5980,3.5970,3.5900,3.5840,3.5840,3.5780,3.5820,3.5790,
3.5750,3.5720,3.5680,3.5650,3.5620,3.5590,3.5560,3.5530,3.5490,3.5470,3.5450,3.5420,
3.5400,3.5370,3.5340,3.5330,3.5300,3.5270,3.5260,3.5240,3.5220,3.5200,3.5180,3.5170,
3.5150,3.5130,3.5120,3.5090,3.5090,3.5060,3.5050,3.5030,3.5020,3.5010,3.5000,3.4990,
3.4970,3.4960,3.4960,3.4960,3.4930,3.4920,3.4920,3.4900,3.4880,3.4870
};
  
  Extcoeff = {3.5889,3.9932,4.5958,5.1961,5.3124,4.2900,3.6239,3.2824,3.0935,2.9573,
  2.9078,2.9135,2.1403,0.9840,0.5031,0.3263,0.2413,0.1769,0.1377,0.1120,0.0954,0.0791,
  0.0702,0.0598,0.0538,0.0485,0.0438,0.0395,0.0348,0.0299,0.0280,0.0266,0.0237,0.0219,
  0.0201,0.0185,0.0173,0.0168,0.0163,0.0147,0.0144,0.0136,0.0128,0.0120,0.0113,0.0106,
  0.0100,0.0093,0.0087,0.0082,0.0076,0.0071,0.0066,0.0061,0.0057,0.0053,0.0049,0.0045,
  0.0041,0.0038,0.0035,0.0032,0.0029,0.0026,0.0023,0.0021,0.0019,0.0017,0.0015,0.0013,
  0.0011,9.8243e-04,8.4060e-04,7.1334e-04,5.9638e-04,4.9020e-04,3.9616e-04,3.1437e-04,
  2.4048e-04,1.7959e-04,1.3043e-04,9.2450e-05,6.7820e-05,5.2168e-05,3.9770e-05,3.0217e-05,
  2.2913e-05,1.7068e-05,1.2382e-05,8.6210e-06,5.6876e-06,3.4275e-06,1.7653e-06,5.5561e-07,
  2.3153e-07,1.3904e-07,8.0863e-08,4.7940e-08,2.7132e-08,1.4318e-08,5.8798e-09,2.3352e-09,
  1.2714e-09,7.5284e-10,4.4799e-10,2.7228e-10,1.5856e-10,8.7196e-11,4.2039e-11,1.8128e-11,
  1.0428e-11,6.2911e-12,3.9030e-12,2.6367e-12,1.7377e-12,1.0428e-12,6.0422e-13,4.2895e-13,
  2.0381e-13,1.3785e-13,1.0901e-13
  };

  mu_r = mu_r_;

  modelType = 4;
}

void ElectroMagnetic::populateSiliconModel() {

double nFF, nSH, kFF, kSH, lambdaumFF, lambdaumSH, lambdastep1, lambdastep2;

double epsrFF_real, epsrFF_imag, epsrSH_real, epsrSH_imag;

int brojac;

lambdaumFF = lambda * 1e6;
lambdaumSH = lambdaumFF / 2.0 ;

 for(brojac = 0; brojac < Extcoeff.size(); brojac++) {
     
     lambdastep1 = 0.25 + brojac * 0.01;   
     lambdastep2 = 0.25 + (brojac+1) * 0.01; 
          if ((lambdaumFF>=lambdastep1)&&(lambdaumFF<=lambdastep2))
                break;                 
     }
     
       nFF = refInd[brojac] + 
                        ((refInd[brojac+1] - refInd[brojac])/(lambdastep2 - lambdastep1)) * 
                         (lambdaumFF-lambdastep1);
                         
       kFF = Extcoeff[brojac] + 
                        ((Extcoeff[brojac+1] - Extcoeff[brojac])/(lambdastep2 - lambdastep1)) * 
                         (lambdaumFF-lambdastep1);       
                         
                         
  for(brojac = 0; brojac < Extcoeff.size(); brojac++) {
     
     lambdastep1 = 0.25 + brojac * 0.01;   
     lambdastep2 = 0.25 + (brojac+1) * 0.01; 
          if ((lambdaumSH>=lambdastep1)&&(lambdaumSH<=lambdastep2))
                break;                 
     }
     
       nSH = refInd[brojac] + 
                        ((refInd[brojac+1] - refInd[brojac])/(lambdastep2 - lambdastep1)) * 
                         (lambdaumSH-lambdastep1);
                         
       kSH = Extcoeff[brojac] + 
                        ((Extcoeff[brojac+1] - Extcoeff[brojac])/(lambdastep2 - lambdastep1)) * 
                         (lambdaumSH-lambdastep1);                               
                                                      

  epsrFF_real = std::pow (nFF, 2) - std::pow (kFF, 2);
  epsrFF_imag = 2.0 * nFF * kFF;
  
  epsrSH_real = std::pow (nSH, 2) - std::pow (kSH, 2);
  epsrSH_imag = 2.0 * nSH * kSH; 
  
  epsilon_r =  epsrFF_real + std::complex<double>(0.0, 1.0)*epsrFF_imag;
                    
  epsilon_r_SH =  epsrSH_real + std::complex<double>(0.0, 1.0)*epsrSH_imag;

  epsilon_SH = epsilon_r_SH * consEpsilon0;
  
  epsilon = epsilon_r * consEpsilon0;

  ksippp = 65e-19;
  ksiparppar = 3.5e-19;
  gamma = 1.3e-19; 
}


void ElectroMagnetic::update(double lambda_) {
  lambda = lambda_;

  if (modelType == 0) // Static model
  {
    // Do nothing as we are in the static case.
  }

  if (modelType == 3) // Hydrodynamic model
  {
    populateHydrodynamicModel();
  }
  if (modelType == 4) // Silicon model
  {
    populateSiliconModel();
  }  
}
