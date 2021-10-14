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
#include "Bessel.h"
#include "HarmonicsIterator.h"
#include "Scatterer.h"
#include "Tools.h"
#include "Trian.h"
#include <Eigen/LU> 

Scatterer::Scatterer(Spherical<double> vR_, ElectroMagnetic elmag_, double radius_, int nMax_, int nMaxS_)
    : vR(vR_), elmag(elmag_), radius(radius_), nMax(nMax_), nMaxS(nMaxS_) {}

Scatterer::Scatterer(int nMax, int nMaxS) : Scatterer({0, 0, 0}, {0, 0, 0, 0, 0, 0}, 0e0, nMax, nMaxS) {}

Scatterer::~Scatterer() {}

void Scatterer::Mesh(std::vector<double> co, std::vector<int> top)	 
{
coord = co;
topol = top;

}

void Scatterer::getTLocal(optimet::Matrix<optimet::t_complex>& Tmatrix,
 optimet::t_real omega_, ElectroMagnetic const &bground) const {
  using namespace optimet;
 
  auto const N = HarmonicsIterator::max_flat(nMax) - 1;
  
  if (this->scatterer_type == "sphere"){
 
  // Fundamental frequency coefficients 
  auto const k_s = omega_ * std::sqrt(elmag.epsilon * elmag.mu);
  auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu);
  

  auto const rho = k_s / k_b;
  auto const r_0 = k_b * radius;
  
  
  auto const mu_sob = elmag.mu / bground.mu;

  auto const Jn = bessel<Bessel>(r_0, nMax);
  auto const Jrho = bessel<Bessel>(rho * r_0, nMax);
  auto const Hn = bessel<Hankel1>(r_0, nMax);
  auto const Hnrho = bessel<Hankel1>(rho * r_0, nMax);
   

  std::complex<double> zeta_b2 = std::sqrt(bground.mu / bground.epsilon);
  std::complex<double> zeta_j2 = std::sqrt(elmag.mu / elmag.epsilon); 
  std::complex<double> zeta_boj2 = zeta_b2 / zeta_j2;  

  auto const N = HarmonicsIterator::max_flat(nMax) - 1; 
  
  Vector<t_complex> result = Vector<t_complex>::Zero(2 * N);
  std::vector<std::complex<double>> data, ddata;
  
  for(t_uint n(1), current(0); n <= nMax; current += 2 * n + 1, ++n) {
    
    // fundamental frequency functions
    auto const psi = r_0 * std::get<0>(Jn)[n];
    
    auto const dpsi = r_0 * std::get<1>(Jn)[n] + std::get<0>(Jn)[n];
    
    auto const ksi = r_0 * std::get<0>(Hn)[n];
    auto const dksi = r_0 * std::get<1>(Hn)[n] + std::get<0>(Hn)[n];

    auto const psirho = r_0 * rho * std::get<0>(Jrho)[n];
    auto const dpsirho = r_0 * rho * std::get<1>(Jrho)[n] + std::get<0>(Jrho)[n];

    // TE Part  b_n coefficients fundamental frequency
    auto const TE = (psi / ksi) * (mu_sob * dpsi / psi - rho * dpsirho / psirho) /
                    (rho * dpsirho / psirho - mu_sob * dksi / ksi);
                  
                                        
     
    result.segment(current, 2 * n + 1).fill(TE);
    

    // TM part  a_n coefficients fundamental frequency
    auto const TM = (psi / ksi) * (mu_sob * dpsirho / psirho - rho * dpsi / psi) /
                    (rho * dksi / ksi - mu_sob * dpsirho / psirho);  
                    
                                                
                    
    result.segment(current + N, 2 * n + 1).fill(TM);
    
}
  	
  Tmatrix = result.asDiagonal();
  }
  
  // calculate the Tmatrix for arbitrary shaped scatterer (evaluation of surface integrals)
  
  else if (this->scatterer_type == "arbitrary.shape"){
  
          CompoundIterator mu1;
          CompoundIterator mup;
          CompoundIterator nu1;
          
          Cartesian<double> intpoinCar; // integration points in Cartesian system
          Spherical<double> intpoinSph; // integration points in spherical system
          
          const std::complex<double> k_0 = (omega_) * std::sqrt(consEpsilon0 * consMu0);
          
          int muMax =mu1.max(nMax);
          int nuMax =nu1.max(nMax);
          
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> resDOT;
          
          Matrix<t_complex> Q = Matrix<t_complex>::Zero(2 * N, 2 * N);
          Matrix<t_complex> RgQ = Matrix<t_complex>::Zero(2 * N, 2 * N);
          
          auto const k_s = omega_ * std::sqrt(elmag.epsilon * elmag.mu);
          auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu);
  
        unsigned int Nt = (topol.size()) / 3;  //number of triangles

        std::vector<std::vector<double>> Points = Tools::getPoints4();

	std::vector<double> Weights = Tools::getWeights4();
       
	for (mu1 = 0; mu1 < muMax ; mu1++){
	
	int nn = mu1.first;
	int mm =  mu1.second;
	mup = nn * (nn + 1) + mm - 1;
		
	for (nu1 = 0; nu1 < nuMax ; nu1++){
	
	std::complex<double> I11(0.0,0.0), I12(0.0,0.0), RgI11(0.0,0.0), RgI12(0.0,0.0), I21(0.0,0.0), I22(0.0,0.0), RgI21(0.0,0.0), RgI22(0.0,0.0);
        std::complex<double> I31(0.0,0.0), I32(0.0,0.0), RgI31(0.0,0.0), RgI32(0.0,0.0), I41(0.0,0.0), I42(0.0,0.0), RgI41(0.0,0.0), RgI42(0.0,0.0);
	
	// surface integration
	 for (int ele1 = 0; ele1 < Nt; ++ele1){

        const int* n1 = this->getNOvertex(ele1);
	
		const double* p1 = this->getCoord(n1[0]);
                
		const double* p2 = this->getCoord(n1[1]);
	
		const double* p3 = this->getCoord(n1[2]);
                
                Trian trian(p1, p2, p3); 
                
		double det = trian.getDeter();
		std::vector<double> nvec = trian.getnorm();
	
	for (int ni = 0; ni != Weights.size(); ++ni) {
	
	double N1 = Points[ni][0];
	double N2 = Points[ni][1];
        double N0 = 1.0 - N1 - N2;

	double wi = Weights[ni];
         
	intpoinCar.x = (p1[0])*N1 + (p2[0])*N2 + (p3[0])*N0;
	intpoinCar.y = (p1[1])*N1 + (p2[1])*N2 + (p3[1])*N0;
	intpoinCar.z = (p1[2])*N1 + (p2[2])*N2 + (p3[2])*N0;
	
	intpoinSph = Tools::toSpherical(intpoinCar); 
	

        optimet::AuxCoefficients aCoefext3(intpoinSph, k_b, 0, nMax); //radiative VSWF (3)
         optimet::AuxCoefficients aCoefext1(intpoinSph, k_b, 1, nMax); //regular VSWF (1)
        optimet::AuxCoefficients aCoefint1(intpoinSph, k_0 * sqrt(elmag.epsilon_r * elmag.mu_r), 1, nMax); //regular VSWF (1)

        //*****************************************
        Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
		
        resDOT = Tools::dot (&nvec[0], resCR);

        I11 = I11 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I12 = I12 + wi * det *(consCi * k_b*k_s)*resDOT;
	//*****************************************
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI11 = RgI11 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI12 = RgI12 + wi * det *(consCi * k_b*k_s)*resDOT;
	
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I21 = I21 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I22 = I22 + wi * det *(consCi * k_b*k_s)*resDOT;
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI21 = RgI21 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI22 = RgI22 + wi * det *(consCi * k_b*k_s)*resDOT;
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I31 = I31 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I32 = I32 + wi * det *(consCi * k_b*k_s)*resDOT;
	//*****************************************
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI31 = RgI31 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI32 = RgI32 + wi * det *(consCi * k_b*k_s)*resDOT;
	
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I41 = I41 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I42 = I42 + wi * det *(consCi * k_b*k_s)*resDOT;
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI41 = RgI41 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI42 = RgI42 + wi * det *(consCi * k_b*k_s)*resDOT;
	
	
	}// Gauss Legendre integration of spherical harmonics over a triangle 
        }// iteration over all triangles, complete surface integration	

       Q (mu1 , nu1) = (I11 + I12);
	Q (mu1 , nu1 + nuMax) = (I21 + I22);
	Q (mu1 + muMax , nu1) = (I31 + I32);
	Q (mu1 + muMax , nu1 + nuMax) = (I41 + I42);
	
	RgQ (mu1 , nu1) = (RgI11 + RgI12);
	RgQ (mu1 , nu1 + nuMax) = (RgI21 + RgI22);
	RgQ (mu1 + muMax , nu1) = (RgI31 + RgI32);
	RgQ (mu1 + muMax , nu1 + nuMax) = (RgI41 + RgI42);
         
	}
	
	}// iteration over all spherical harmonics
           
  
  Tmatrix = - RgQ * Q.inverse();
  
  }
   
}


void Scatterer::getQLocal(optimet::Matrix<optimet::t_complex>& Intrmatrix, optimet::t_real omega_, ElectroMagnetic const &bground) const {
  using namespace optimet;
  
  auto const N = HarmonicsIterator::max_flat(nMax) - 1;

// calculate the Rg(Q) matrix for arbitrary shaped scatterer, needed for conversion of scattered to internal coefficients
          CompoundIterator mu1;
          CompoundIterator mup;
          CompoundIterator nu1;
          
          Cartesian<double> intpoinCar; // integration points in Cartesian system
          Spherical<double> intpoinSph; // integration points in spherical system
         const std::complex<double> k_0 = (omega_) * std::sqrt(consEpsilon0 * consMu0);
          
          int muMax =mu1.max(nMax);
          int nuMax =nu1.max(nMax);
          
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> resDOT;
          
          Matrix<t_complex> RgQ = Matrix<t_complex>::Zero(2 * N, 2 * N);
          
          auto const k_s = omega_ * std::sqrt(elmag.epsilon * elmag.mu);  
          auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu); 
  
        unsigned int Nt = (topol.size()) / 3;  //number of triangles

       std::vector<std::vector<double>> Points = Tools::getPoints4();

	std::vector<double> Weights = Tools::getWeights4();

	for (mu1 = 0; mu1 < muMax ; mu1++){
	
	int nn = mu1.first;
	int mm =  mu1.second;
	mup = nn * (nn + 1) + mm - 1;
		
	for (nu1 = 0; nu1 < nuMax ; nu1++){
	
	std::complex<double> RgI11(0.0,0.0), RgI12(0.0,0.0), RgI21(0.0,0.0), RgI22(0.0,0.0);
        std::complex<double> RgI31(0.0,0.0), RgI32(0.0,0.0), RgI41(0.0,0.0), RgI42(0.0,0.0);
	
	// surface integration
	 for (int ele1 = 0; ele1 < Nt; ++ele1){

         const int* n1 = this->getNOvertex(ele1);
	
		const double* p1 = this->getCoord(n1[0]);

		const double* p2 = this->getCoord(n1[1]);
	
		const double* p3 = this->getCoord(n1[2]);
                
                Trian trian(p1, p2, p3); 
                
		double det = trian.getDeter();
		std::vector<double> nvec = trian.getnorm();
		
	for (int ni = 0; ni != Weights.size(); ++ni) {
	
	double N1 = Points[ni][0];
	double N2 = Points[ni][1];
        double N0 = 1.0 - N1 - N2;

	double wi = Weights[ni];
         
	intpoinCar.x = (p1[0])*N1 + (p2[0])*N2 + (p3[0])*N0;
	intpoinCar.y = (p1[1])*N1 + (p2[1])*N2 + (p3[1])*N0;
	intpoinCar.z = (p1[2])*N1 + (p2[2])*N2 + (p3[2])*N0;
	
	intpoinSph = Tools::toSpherical(intpoinCar); 
	

        optimet::AuxCoefficients aCoefext1(intpoinSph, k_b, 1, nMax); //regular VSWF (1)
        optimet::AuxCoefficients aCoefint1(intpoinSph, k_0 * sqrt(elmag.epsilon_r * elmag.mu_r), 1, nMax); //regular VSWF (1)
	

	//*****************************************
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI11 = RgI11 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI12 = RgI12 + wi * det *(consCi * k_b*k_s)*resDOT;
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI21 = RgI21 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI22 = RgI22 + wi * det *(consCi * k_b*k_s)*resDOT;
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI31 = RgI31 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI32 = RgI32 + wi * det *(consCi * k_b*k_s)*resDOT;
	
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI41 = RgI41 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI42 = RgI42 + wi * det *(consCi * k_b*k_s)*resDOT;
	
	
	}// Gauss Legendre integration of spherical harmonics over a triangle 
        }// iteration over all triangles, complete surface integration

       RgQ (mu1 , nu1) = (RgI11 + RgI12);
	RgQ (mu1 , nu1 + nuMax) = (RgI21 + RgI22);
	RgQ (mu1 + muMax , nu1) = (RgI31 + RgI32);
	RgQ (mu1 + muMax , nu1 + nuMax) = (RgI41 + RgI42);
         
	}
	
	}// iteration over all spherical harmonics

    Intrmatrix = RgQ.inverse();
  
  }
                  
void Scatterer::getTLocalSH(optimet::Matrix<optimet::t_complex>& Tmatrix, optimet::t_real omega_, ElectroMagnetic const &bground) const {
  using namespace optimet;
  
  auto const N = HarmonicsIterator::max_flat(nMaxS) - 1;
 
  if (this->scatterer_type == "sphere"){
  
  auto const k_s = 2.0 * omega_ * std::sqrt(elmag.epsilon_SH * elmag.mu_SH);
  auto const k_b = 2.0 * omega_ * std::sqrt(bground.epsilon * bground.mu);
  auto const rho = k_s / k_b;
  auto const r_0 = k_b * radius;
  
  auto const mu_sob = elmag.mu_SH / bground.mu;

  auto const Jn = bessel<Bessel>(r_0, nMaxS);
  auto const Jrho = bessel<Bessel>(rho * r_0, nMaxS);
  auto const Hn = bessel<Hankel1>(r_0, nMaxS);
  auto const Hnrho = bessel<Hankel1>(rho * r_0, nMaxS);
   

  std::complex<double> zeta_b2 = std::sqrt(bground.mu / bground.epsilon);
  std::complex<double> zeta_j2 = std::sqrt(elmag.mu_SH / elmag.epsilon_SH); 
  std::complex<double> zeta_boj2 = zeta_b2 / zeta_j2;   
  
  Vector<t_complex> result = Vector<t_complex>::Zero(2 * N);
 
  std::vector<std::complex<double>> data, ddata;
  
  for(t_uint n(1), current(0); n <= nMaxS; current += 2 * n + 1, ++n) {
    
    
    auto const psi = r_0 * std::get<0>(Jn)[n];
    auto const dpsi = r_0 * std::get<1>(Jn)[n] + std::get<0>(Jn)[n];
    
    auto const ksi = r_0 * std::get<0>(Hn)[n];
    auto const dksi = r_0 * std::get<1>(Hn)[n] + std::get<0>(Hn)[n];
    
    auto const psirho = r_0 * rho * std::get<0>(Jrho)[n];
    auto const dpsirho = r_0 * rho * std::get<1>(Jrho)[n] + std::get<0>(Jrho)[n];
    

    // TE Part  b_n coefficients fundamental frequency
    auto const TE = ((psi / ksi) * (mu_sob * dpsi / psi - rho * dpsirho / psirho) /
                    (rho * dpsirho / psirho - mu_sob * dksi / ksi));
                                  
    result.segment(current, 2 * n + 1).fill(TE);

    // TM part  a_n coefficients fundamental frequency
    auto const TM = ((psi / ksi) * (mu_sob * dpsirho / psirho - rho * dpsi / psi) /
                    (rho * dksi / ksi - mu_sob * dpsirho / psirho));  
                                                             
    result.segment(current + N, 2 * n + 1).fill(TM);

  }
  	
  Tmatrix = result.asDiagonal();
  
  	
  }
  
  }

optimet::Vector<optimet::t_complex>
Scatterer::getTLocalSH1_outer(optimet::t_real omega_, ElectroMagnetic const &bground) const {
  using namespace optimet;
  
  
    // SH frequency coefficients 
  auto const k_s_SH = 2 * omega_ * std::sqrt(elmag.epsilon_SH * elmag.mu_SH);
  auto const k_b_SH = 2 * omega_ * std::sqrt(bground.epsilon * bground.mu);
  

   
  auto const rho_SH = k_s_SH / k_b_SH;
  auto const r_0_SH = k_b_SH * radius;
  
 
  auto const Jn_SH = bessel<Bessel>(r_0_SH, nMaxS);
  auto const Jrho_SH = bessel<Bessel>(rho_SH * r_0_SH, nMaxS);
  auto const Hn_SH = bessel<Hankel1>(r_0_SH, nMaxS);
  auto const Hnrho_SH = bessel<Hankel1>(rho_SH * r_0_SH, nMaxS);
  
  std::complex<double> x_b2 = k_b_SH * radius;
  std::complex<double> zeta_b2 = std::sqrt(bground.mu / bground.epsilon);
  std::complex<double> zeta_j2 = std::sqrt(elmag.mu_SH / elmag.epsilon_SH); 
  std::complex<double> zeta_boj2 = zeta_b2 / zeta_j2;  

  auto const N = HarmonicsIterator::max_flat(nMaxS) - 1;
  
  
  Vector<t_complex> resultSH1 = Vector<t_complex>::Zero(2 * N);
 
  for(t_uint n(1), current(0); n <= nMaxS; current += 2 * n + 1, ++n) {
    
    
      // SH frequency functions
    auto const psi_SH = r_0_SH * std::get<0>(Jn_SH)[n];
    auto const dpsi_SH = r_0_SH * std::get<1>(Jn_SH)[n] + std::get<0>(Jn_SH)[n];

    auto const ksi_SH = r_0_SH * std::get<0>(Hn_SH)[n];
    auto const dksi_SH = r_0_SH * std::get<1>(Hn_SH)[n] + std::get<0>(Hn_SH)[n];

    auto const psirho_SH = r_0_SH * rho_SH * std::get<0>(Jrho_SH)[n];
    auto const dpsirho_SH = r_0_SH * rho_SH * std::get<1>(Jrho_SH)[n] + std::get<0>(Jrho_SH)[n];
    
     
    // TE1 Part  b_n' coefficients SH frequency                
    auto const TE_SH1 = - x_b2 * psirho_SH / (zeta_boj2 * ksi_SH * dpsirho_SH - psirho_SH * dksi_SH); 
    

    resultSH1.segment(current, 2 * n + 1).fill(TE_SH1);
 
                
    // TM1 Part  a_n' coefficients SH frequency  
    
     auto const TM_SH1 = - x_b2 * dpsirho_SH / (zeta_boj2 * psirho_SH * dksi_SH - ksi_SH * dpsirho_SH);
    

    resultSH1.segment(current + N, 2 * n + 1).fill(TM_SH1);
    
    
  }
  
  return resultSH1;
}

 
optimet::Vector<optimet::t_complex>
Scatterer::getTLocalSH2_outer(optimet::t_real omega_, ElectroMagnetic const &bground) const {
  using namespace optimet;
  
  
    // SH frequency coefficients 
  auto const k_s_SH = 2 * omega_ * std::sqrt(elmag.epsilon_SH * elmag.mu_SH); 
  auto const k_b_SH = 2 * omega_ * std::sqrt(bground.epsilon * bground.mu);

  auto const rho_SH = k_s_SH / k_b_SH;
  auto const r_0_SH = k_b_SH * radius;
 
  auto const Jn_SH = bessel<Bessel>(r_0_SH, nMaxS);
  auto const Jrho_SH = bessel<Bessel>(rho_SH * r_0_SH, nMaxS);
  auto const Hn_SH = bessel<Hankel1>(r_0_SH, nMaxS);
  auto const Hnrho_SH = bessel<Hankel1>(rho_SH * r_0_SH, nMaxS);
  
  std::complex<double> x_b2 = k_b_SH * radius;
  std::complex<double> zeta_b2 = std::sqrt(bground.mu / bground.epsilon);
  std::complex<double> zeta_j2 = std::sqrt(elmag.mu_SH / elmag.epsilon_SH);
  std::complex<double> zeta_boj2 = zeta_b2 / zeta_j2;  

  auto const N = HarmonicsIterator::max_flat(nMaxS) - 1;
 
  Vector<t_complex> resultSH2 = Vector<t_complex>::Zero(2 * N);
 
  for(t_uint n(1), current(0); n <= nMaxS; current += 2 * n + 1, ++n) {
    
    
      // SH frequency functions
    auto const psi_SH = r_0_SH * std::get<0>(Jn_SH)[n];
    auto const dpsi_SH = r_0_SH * std::get<1>(Jn_SH)[n] + std::get<0>(Jn_SH)[n];

    auto const ksi_SH = r_0_SH * std::get<0>(Hn_SH)[n];
    auto const dksi_SH = r_0_SH * std::get<1>(Hn_SH)[n] + std::get<0>(Hn_SH)[n];

    auto const psirho_SH = r_0_SH * rho_SH * std::get<0>(Jrho_SH)[n];
    auto const dpsirho_SH = r_0_SH * rho_SH * std::get<1>(Jrho_SH)[n] + std::get<0>(Jrho_SH)[n];
    
     
    // TE 2 Part  b_n'' coefficients SH frequency 
    auto const TE_SH2 = zeta_boj2 * x_b2 * dpsirho_SH / (zeta_boj2 * ksi_SH * dpsirho_SH - psirho_SH * dksi_SH);
              
    resultSH2.segment(current, 2 * n + 1).fill(TE_SH2);
    
  
    // TM 2 Part  a_n'' coefficients SH frequency  
    
     auto const TM_SH2 =   zeta_boj2 * x_b2 * psirho_SH / (zeta_boj2 * psirho_SH * dksi_SH - ksi_SH * dpsirho_SH);
                
    resultSH2.segment(current + N, 2 * n + 1).fill(TM_SH2);
    
  }
  return resultSH2;
}

void Scatterer::getTLocalSH_ARB(optimet::Matrix<optimet::t_complex>& Tmatrix, optimet::t_real omega_, ElectroMagnetic const &bground) const {
  using namespace optimet;
  
  auto const N = HarmonicsIterator::max_flat(nMaxS) - 1;
  
  // SH frequency coefficients 
  auto const k_s_SH = 2 * omega_ * std::sqrt(elmag.epsilon_SH * elmag.mu_SH);
  auto const k_b_SH = 2 * omega_ * std::sqrt(bground.epsilon * bground.mu);
  
          CompoundIterator mu1;
          CompoundIterator mup;
          CompoundIterator nu1;
          
          Cartesian<double> intpoinCar; // integration points in Cartesian system
          Spherical<double> intpoinSph; // integration points in spherical system 
          const std::complex<double> k_0_SH = (2 * omega_) * std::sqrt(consEpsilon0 * consMu0);
           
          int muMax =mu1.max(nMaxS);
          int nuMax =nu1.max(nMaxS);
          
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> resDOT;
          
          Matrix<t_complex> Q = Matrix<t_complex>::Zero(2 * N, 2 * N);
          Matrix<t_complex> RgQ = Matrix<t_complex>::Zero(2 * N, 2 * N);
       
        unsigned int Nt = (topol.size()) / 3;  //number of triangles

        std::vector<std::vector<double>> Points = Tools::getPoints4();

	std::vector<double> Weights = Tools::getWeights4();
       
	for (mu1 = 0; mu1 < muMax ; mu1++){
         
	int nn = mu1.first;
	int mm =  mu1.second;
	mup = nn * (nn + 1) + mm - 1;
        		
 	for (nu1 = 0; nu1 < nuMax ; nu1++){
	
	std::complex<double> I11(0.0,0.0), I12(0.0,0.0), RgI11(0.0,0.0), RgI12(0.0,0.0), I21(0.0,0.0), I22(0.0,0.0), RgI21(0.0,0.0), RgI22(0.0,0.0);
        std::complex<double> I31(0.0,0.0), I32(0.0,0.0), RgI31(0.0,0.0), RgI32(0.0,0.0), I41(0.0,0.0), I42(0.0,0.0), RgI41(0.0,0.0), RgI42(0.0,0.0);
	
	// surface integration
	for (int ele1 = 0; ele1 < Nt; ++ele1){
        
        const int* n1 = this->getNOvertex(ele1);
	
		const double* p1 = this->getCoord(n1[0]);
                
		const double* p2 = this->getCoord(n1[1]);
	
		const double* p3 = this->getCoord(n1[2]);
                
                Trian trian(p1, p2, p3); 
                
		double det = trian.getDeter();
		std::vector<double> nvec = trian.getnorm();
		
	for (int ni = 0; ni != Weights.size(); ++ni) {
	
	double N1 = Points[ni][0];
	double N2 = Points[ni][1];
        double N0 = 1.0 - N1 - N2;

	double wi = Weights[ni];
         
	intpoinCar.x = (p1[0])*N1 + (p2[0])*N2 + (p3[0])*N0;
	intpoinCar.y = (p1[1])*N1 + (p2[1])*N2 + (p3[1])*N0;
	intpoinCar.z = (p1[2])*N1 + (p2[2])*N2 + (p3[2])*N0;
	
	intpoinSph = Tools::toSpherical(intpoinCar); 
	

        optimet::AuxCoefficients aCoefext3(intpoinSph, k_b_SH, 0, nMaxS); //radiative VSWF (3)
        optimet::AuxCoefficients aCoefext1(intpoinSph, k_b_SH, 1, nMaxS); //regular VSWF (1)
	optimet::AuxCoefficients aCoefint1(intpoinSph, k_0_SH * sqrt(elmag.epsilon_r_SH * elmag.mu_r_SH), 1, nMaxS); //regular VSWF (1)
        
        //*****************************************
        Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I11 = I11 + wi * det *(k_b_SH)*resDOT;
	///***********
        Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I12 = I12 + wi * det *(k_s_SH)*resDOT;
	//*****************************************
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI11 = RgI11 + wi * det *(k_b_SH)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI12 = RgI12 + wi * det *(k_s_SH)*resDOT;
	
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I21 = I21 + wi * det *(k_b_SH)*resDOT;
	///***********
        Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I22 = I22 + wi * det *(k_s_SH)*resDOT;
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI21 = RgI21 + wi * det *(k_b_SH)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI22 = RgI22 + wi * det *(k_s_SH)*resDOT;
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I31 = I31 + wi * det *(k_b_SH)*resDOT;
	///***********
        Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I32 = I32 + wi * det *(k_s_SH)*resDOT;
	//*****************************************
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI31 = RgI31 + wi * det *(k_b_SH)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI32 = RgI32 + wi * det *(k_s_SH)*resDOT;
	
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I41 = I41 + wi * det *(k_b_SH)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I42 = I42 + wi * det *(k_s_SH)*resDOT;
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI41 = RgI41 + wi * det *(k_b_SH)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI42 = RgI42 + wi * det *(k_s_SH)*resDOT;
	
	
	}// Gauss Legendre integration of spherical harmonics over a triangle 
        }// iteration over all triangles, complete surface integration	

       Q (mu1 , nu1) = (I11 + I12);
	Q (mu1 , nu1 + nuMax) = (I21 + I22);
	Q (mu1 + muMax , nu1) = (I31 + I32);
	Q (mu1 + muMax , nu1 + nuMax) = (I41 + I42);
	
	RgQ (mu1 , nu1) = (RgI11 + RgI12);
	RgQ (mu1 , nu1 + nuMax) = (RgI21 + RgI22);
	RgQ (mu1 + muMax , nu1) = (RgI31 + RgI32);
	RgQ (mu1 + muMax , nu1 + nuMax) = (RgI41 + RgI42);
         
	}
	
	}// iteration over all spherical harmonics

       Tmatrix = RgQ * Q.inverse();
  
  }


void Scatterer::getQLocalSH(optimet::Matrix<optimet::t_complex>& Intrmatrix, optimet::t_real omega_, ElectroMagnetic const &bground) const {
  using namespace optimet;
  
  auto const N = HarmonicsIterator::max_flat(nMaxS) - 1;

  // calculate the Rg(Q) matrix for arbitrary shaped scatterer, needed for conversion of scattered to internal coefficients
  
  CompoundIterator mu1;
          CompoundIterator mup;
          CompoundIterator nu1;
          
          Cartesian<double> intpoinCar; // integration points in Cartesian system
          Spherical<double> intpoinSph; // integration points in spherical system
          const std::complex<double> k_0_SH = 2 * (omega_) * std::sqrt(consEpsilon0 * consMu0);
          
          int muMax =mu1.max(nMaxS);
          int nuMax =nu1.max(nMaxS);
          
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> resDOT;
          
         auto const k_s_SH = 2 * omega_ * std::sqrt(elmag.epsilon_SH * elmag.mu_SH); 
         auto const k_b_SH = 2 * omega_ * std::sqrt(bground.epsilon * bground.mu);
  
        unsigned int Nt = (topol.size()) / 3;  //number of triangles
        
        std::vector<std::vector<double>> Points = Tools::getPoints4();

	std::vector<double> Weights = Tools::getWeights4();

	for (mu1 = 0; mu1 < muMax ; mu1++){
	
	int nn = mu1.first;
	int mm =  mu1.second;
	mup = nn * (nn + 1) + mm - 1;
		
	for (nu1 = 0; nu1 < nuMax ; nu1++){
	
	std::complex<double> RgI11(0.0,0.0), RgI12(0.0,0.0), RgI21(0.0,0.0), RgI22(0.0,0.0);
        std::complex<double> RgI31(0.0,0.0), RgI32(0.0,0.0), RgI41(0.0,0.0), RgI42(0.0,0.0);
	
	// surface integration
	for (int ele1 = 0; ele1 < Nt; ++ele1){

        const int* n1 = this->getNOvertex(ele1);
	
		const double* p1 = this->getCoord(n1[0]);

		const double* p2 = this->getCoord(n1[1]);
	
		const double* p3 = this->getCoord(n1[2]);
                
                Trian trian(p1, p2, p3); 
                
		double det = trian.getDeter();
		std::vector<double> nvec = trian.getnorm();
		
	for (int ni = 0; ni != Weights.size(); ++ni) {
	
	double N1 = Points[ni][0];
	double N2 = Points[ni][1];
        double N0 = 1.0 - N1 - N2;

	double wi = Weights[ni];
         
	intpoinCar.x = (p1[0])*N1 + (p2[0])*N2 + (p3[0])*N0;
	intpoinCar.y = (p1[1])*N1 + (p2[1])*N2 + (p3[1])*N0;
	intpoinCar.z = (p1[2])*N1 + (p2[2])*N2 + (p3[2])*N0;
	
	intpoinSph = Tools::toSpherical(intpoinCar); 
	

        optimet::AuxCoefficients aCoefext1(intpoinSph, k_b_SH, 1, nMaxS); //regular VSWF (1)
        optimet::AuxCoefficients aCoefint1(intpoinSph, k_0_SH * sqrt(elmag.epsilon_r_SH * elmag.mu_r_SH), 1, nMaxS); //regular VSWF (1)

       //*****************************************
       Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI11 = RgI11 + wi * det *(k_b_SH)*resDOT;
       ///***********
       Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI12 = RgI12 + wi * det *(k_s_SH)*resDOT;
        //*****************************************
        Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI21 = RgI21 + wi * det *(k_b_SH)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI22 = RgI22 + wi * det *(k_s_SH)*resDOT;
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI31 = RgI31 + wi * det *(k_b_SH)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI32 = RgI32 + wi * det *(k_s_SH)*resDOT;
	
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext1.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI41 = RgI41 + wi * det *(k_b_SH)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext1.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	RgI42 = RgI42 + wi * det *(k_s_SH)*resDOT;
	
	
	}// Gauss Legendre integration of spherical harmonics over a triangle 
        }// iteration over all triangles, complete surface integration	

       Intrmatrix (mu1 , nu1) = (RgI11 + RgI12);
	Intrmatrix (mu1 , nu1 + nuMax) = (RgI21 + RgI22);
	Intrmatrix (mu1 + muMax , nu1) = (RgI31 + RgI32);
	Intrmatrix (mu1 + muMax , nu1 + nuMax) = (RgI41 + RgI42);
         
	}
	
	}// iteration over all spherical harmonics

     }
         

optimet::Vector<optimet::t_complex>
Scatterer::getIaux(optimet::t_real omega_, ElectroMagnetic const &bground) const {

  auto const k_s = omega_ * std::sqrt(elmag.epsilon * elmag.mu);
  auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu);
  auto const rho = k_s / k_b;
  auto const r_0 = k_b * radius;
  auto const mu_j = elmag.mu;
  auto const mu_0 = bground.mu;
  
  std::complex<double> zeta_b2 = std::sqrt(bground.mu / bground.epsilon);
  std::complex<double> zeta_j2 = std::sqrt(elmag.mu / elmag.epsilon); 
  std::complex<double> zeta_boj2 = zeta_b2 / zeta_j2;  

  auto Jdata = optimet::bessel<optimet::Bessel>(r_0, nMax);
  auto const Jrho = optimet::bessel<optimet::Bessel>(rho * r_0, nMax);

  optimet::Vector<optimet::t_complex> result(2 * nMax * (nMax + 2));
  auto TE = result.head(nMax * (nMax + 2));
  auto TM = result.tail(nMax * (nMax + 2));
  
  for(auto n = 1, i = 0; n <= nMax; ++n) {
    // obtain Riccati-Bessel functions
    auto const psi = r_0 * std::get<0>(Jdata)[n];
    auto const dpsi = r_0 * std::get<1>(Jdata)[n] + std::get<0>(Jdata)[n];
    auto const psirho = r_0 * rho * std::get<0>(Jrho)[n];
    auto const dpsirho = r_0 * rho * std::get<1>(Jrho)[n] + std::get<0>(Jrho)[n];

    for(auto m = -n; m <= n; ++m, ++i) {
    
    
      TE(i) = (mu_j * rho) / (mu_0 * rho * dpsirho * psi - mu_j * psirho * dpsi) *
              std::complex<double>(0., 1.);
                  
              
      TM(i) = (mu_j * rho) / (mu_j * psi * dpsirho - mu_0 * rho * psirho * dpsi) *
              std::complex<double>(0., 1.);
               
            
              
    }
  }
  return result;
}


optimet::Vector<optimet::t_complex>
Scatterer::getIauxSH1(optimet::t_real omega_, ElectroMagnetic const &bground) const {

  auto const k_s_SH = 2.0 * omega_ * std::sqrt(elmag.epsilon_SH * elmag.mu_SH);
  auto const k_b_SH = 2.0 * omega_ * std::sqrt(bground.epsilon * bground.mu);
  auto const rho_SH = k_s_SH / k_b_SH;
  auto const r_0_SH = k_b_SH * radius;
  auto const mu_0 = bground.mu;
  
  std::complex<double> x_b2 = k_b_SH * radius;
  std::complex<double> x_i2 = k_s_SH * radius;
  std::complex<double> zeta_b2 = std::sqrt(bground.mu / bground.epsilon);
  std::complex<double> zeta_j2 = std::sqrt(elmag.mu_SH / elmag.epsilon_SH);
  std::complex<double> zeta_boj2 = zeta_b2 / zeta_j2;  

  auto Jdata_SH = optimet::bessel<optimet::Bessel>(r_0_SH, nMaxS);
  auto const Jrho_SH = optimet::bessel<optimet::Bessel>(rho_SH * r_0_SH, nMaxS);
  auto const Hn_SH = optimet::bessel<optimet::Hankel1>(r_0_SH, nMaxS);

  optimet::Vector<optimet::t_complex> result(2 * nMaxS * (nMaxS + 2));
  auto TE = result.head( nMaxS * (nMaxS + 2));
  auto TM = result.tail( nMaxS * (nMaxS + 2));
  
  for(auto n = 1, i = 0; n <= nMaxS; ++n) {
    // obtain Riccati-Bessel functions
    auto const psi_SH = r_0_SH * std::get<0>(Jdata_SH)[n];
    auto const dpsi_SH = r_0_SH * std::get<1>(Jdata_SH)[n] + std::get<0>(Jdata_SH)[n];
    auto const psirho_SH = r_0_SH * rho_SH * std::get<0>(Jrho_SH)[n];
    auto const dpsirho_SH = r_0_SH * rho_SH * std::get<1>(Jrho_SH)[n] + std::get<0>(Jrho_SH)[n];
    auto const ksi_SH = r_0_SH * std::get<0>(Hn_SH)[n];
    auto const dksi_SH = r_0_SH * std::get<1>(Hn_SH)[n] + std::get<0>(Hn_SH)[n];

    for(auto m = -n; m <= n; ++m, ++i) {
    
      TE(i) = (- x_i2 * ksi_SH) / (x_b2 * psirho_SH) ; // result related to the bmn'
             
      TM(i) = (- x_i2 * dksi_SH) / (x_b2 * dpsirho_SH); // result related to the amn'       
                     
    }
  }
  
  
  return result;
}


optimet::Vector<optimet::t_complex>
Scatterer::getIauxSH2(optimet::t_real omega_, ElectroMagnetic const &bground) const {

  auto const k_s_SH = 2.0 * omega_ * std::sqrt(elmag.epsilon_SH * elmag.mu_SH);
  auto const k_b_SH = 2.0 * omega_ * std::sqrt(bground.epsilon * bground.mu);
  auto const rho_SH = k_s_SH / k_b_SH;
  auto const r_0_SH = k_b_SH * radius;
  auto const mu_0 = bground.mu;
  
  std::complex<double> bnpp, anpp;
  std::complex<double> x_b2 = k_b_SH * radius;
  std::complex<double> x_i2 = k_s_SH * radius;
  std::complex<double> zeta_b2 = std::sqrt(bground.mu / bground.epsilon);
  std::complex<double> zeta_j2 = std::sqrt(elmag.mu_SH / elmag.epsilon_SH);
  std::complex<double> zeta_boj2 = zeta_b2 / zeta_j2;  

  auto Jdata_SH = optimet::bessel<optimet::Bessel>(r_0_SH, nMaxS);
  auto const Jrho_SH = optimet::bessel<optimet::Bessel>(rho_SH * r_0_SH, nMaxS);
  auto const Hn_SH = optimet::bessel<optimet::Hankel1>(r_0_SH, nMaxS);

  optimet::Vector<optimet::t_complex> result(2 * nMaxS * (nMaxS + 2));
  auto TE = result.head( nMaxS * (nMaxS + 2));
  auto TM = result.tail( nMaxS * (nMaxS + 2));
  
  for(auto n = 1, i = 0; n <= nMaxS; ++n) {
    // obtain Riccati-Bessel functions
    auto const psi_SH = r_0_SH * std::get<0>(Jdata_SH)[n];
    auto const dpsi_SH = r_0_SH * std::get<1>(Jdata_SH)[n] + std::get<0>(Jdata_SH)[n];
    auto const psirho_SH = r_0_SH * rho_SH * std::get<0>(Jrho_SH)[n];
    auto const dpsirho_SH = r_0_SH * rho_SH * std::get<1>(Jrho_SH)[n] + std::get<0>(Jrho_SH)[n];
    auto const ksi_SH = r_0_SH * std::get<0>(Hn_SH)[n];
    auto const dksi_SH = r_0_SH * std::get<1>(Hn_SH)[n] + std::get<0>(Hn_SH)[n];

    for(auto m = -n; m <= n; ++m, ++i) {

      bnpp =   (zeta_boj2 * x_b2 * dpsirho_SH) / (zeta_boj2 * ksi_SH * dpsirho_SH - psirho_SH * dksi_SH);
      
      anpp =    (zeta_boj2 * x_b2 * psirho_SH) / (zeta_boj2 * psirho_SH * dksi_SH - ksi_SH * dpsirho_SH);       
     
      TE(i) =  bnpp * ((- x_i2 * ksi_SH) / (x_b2 * psirho_SH) + (x_i2 * dksi_SH) / (zeta_boj2 * x_b2 * dpsirho_SH)) ;  //related to bmn and cmn 
      
      TM(i) = anpp * ((- x_i2 * dksi_SH) / (x_b2 * dpsirho_SH) + (x_i2 * ksi_SH) / (zeta_boj2 * x_b2 * psirho_SH));  //related to amn and dmn
                     
    }
  }
  return result;
}





