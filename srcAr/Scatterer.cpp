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

#ifdef OPTIMET_MPI
void Scatterer::getQLocal(optimet::Vector<optimet::t_complex>& Qmatrix,
 optimet::t_real omega_, ElectroMagnetic const &bground, int gran1, int gran2) const {
  using namespace optimet;

  auto const N = HarmonicsIterator::max_flat(nMax) - 1;
 
  // calculate the Qmatrix for arbitrary shaped scatterer (evaluation of surface integrals)
  
          CompoundIterator mu1;
          CompoundIterator mup;
          CompoundIterator nu1;
          
          Cartesian<double> intpoinCar; // integration points in Cartesian system
          Spherical<double> intpoinSph; // integration points in spherical system
          
          const std::complex<double> k_0 = (omega_) * std::sqrt(consEpsilon0 * consMu0);
          
          int muMax =mu1.max(nMax);
          int nuMax =nu1.max(nMax);
          int brojac(0);
          int size = nuMax*(gran2 - gran1);
          
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> resDOT;

          Vector<t_complex> Qmatrix_a(size), Qmatrix_b(size), Qmatrix_c(size), Qmatrix_d(size);
         Vector<t_complex> Qmatrix_1(2*size), Qmatrix_2(2*size);
          
          auto const k_s = omega_ * std::sqrt(elmag.epsilon * elmag.mu);
          auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu);
  
        unsigned int Nt = (topol.size()) / 3;  //number of triangles

        std::vector<std::vector<double>> Points = Tools::getPoints4();

	std::vector<double> Weights = Tools::getWeights4();
       
	for (mu1 = gran1; mu1 < gran2 ; mu1++){
	
	int nn = mu1.first;
	int mm =  mu1.second;
	mup = nn * (nn + 1) + mm - 1;
		
	for (nu1 = 0; nu1 < nuMax ; nu1++){
	
	std::complex<double> I11(0.0,0.0), I12(0.0,0.0), I21(0.0,0.0), I22(0.0,0.0);
        std::complex<double> I31(0.0,0.0), I32(0.0,0.0), I41(0.0,0.0), I42(0.0,0.0);
	
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
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I21 = I21 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I22 = I22 + wi * det *(consCi * k_b*k_s)*resDOT;
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I31 = I31 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I32 = I32 + wi * det *(consCi * k_b*k_s)*resDOT;
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I41 = I41 + wi * det *(consCi * k_b*k_b)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I42 = I42 + wi * det *(consCi * k_b*k_s)*resDOT;
	//*****************************************
	
	
	}// Gauss Legendre integration of spherical harmonics over a triangle 
        }// iteration over all triangles, complete surface integration	

        Qmatrix_a (brojac) = (I11 + I12);
        Qmatrix_b (brojac) = (I21 + I22);
        Qmatrix_c (brojac) = (I31 + I32);
        Qmatrix_d (brojac) = (I41 + I42);

        brojac++;         
	}

	}// iteration over all spherical harmonics
  
      Qmatrix_1.head(size) = Qmatrix_a;
      Qmatrix_1.tail(size) = Qmatrix_b;

      Qmatrix_2.head(size) = Qmatrix_c;
      Qmatrix_2.tail(size) = Qmatrix_d;

      Qmatrix.head(2*size) = Qmatrix_1;
      Qmatrix.tail(2*size) = Qmatrix_2;
   
}


void Scatterer::getRgQLocal(optimet::Vector<optimet::t_complex>& RgQmatrix, optimet::t_real omega_, ElectroMagnetic const &bground, int gran1, int gran2) const {
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
          int brojac(0);
          int size = nuMax *(gran2 - gran1);
          
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> resDOT;

          Vector<t_complex> RgQmatrix_a(size), RgQmatrix_b(size), RgQmatrix_c(size), RgQmatrix_d(size);
         Vector<t_complex> RgQmatrix_1(2*size), RgQmatrix_2(2*size);
          
          auto const k_s = omega_ * std::sqrt(elmag.epsilon * elmag.mu);  
          auto const k_b = omega_ * std::sqrt(bground.epsilon * bground.mu); 
  
        unsigned int Nt = (topol.size()) / 3;  //number of triangles

       std::vector<std::vector<double>> Points = Tools::getPoints4();

	std::vector<double> Weights = Tools::getWeights4();

	for (mu1 = gran1; mu1 < gran2 ; mu1++){
	
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

        RgQmatrix_a (brojac) = (RgI11 + RgI12);
        RgQmatrix_b (brojac) = (RgI21 + RgI22);
        RgQmatrix_c (brojac) = (RgI31 + RgI32);
        RgQmatrix_d (brojac) = (RgI41 + RgI42);

        brojac++;        
	}	
	}// iteration over all spherical harmonics

      RgQmatrix_1.head(size) = RgQmatrix_a;
      RgQmatrix_1.tail(size) = RgQmatrix_b;

      RgQmatrix_2.head(size) = RgQmatrix_c;
      RgQmatrix_2.tail(size) = RgQmatrix_d;

      RgQmatrix.head(2*size) = RgQmatrix_1;
      RgQmatrix.tail(2*size) = RgQmatrix_2;
 
  }
                   

void Scatterer::getQLocalSH(optimet::Vector<optimet::t_complex>& QmatrixSH, optimet::t_real omega_, ElectroMagnetic const &bground, int gran1, int gran2) const {
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
          int brojac(0);
          int size = nuMax *(gran2 - gran1);
          
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> resDOT;

         Vector<t_complex> QmatrixSH_a(size), QmatrixSH_b(size), QmatrixSH_c(size), QmatrixSH_d(size);
         Vector<t_complex> QmatrixSH_1(2*size), QmatrixSH_2(2*size);

        unsigned int Nt = (topol.size()) / 3;  //number of triangles

        std::vector<std::vector<double>> Points = Tools::getPoints4();

	std::vector<double> Weights = Tools::getWeights4();
       
	for (mu1 = gran1; mu1 < gran2 ; mu1++){
         
	int nn = mu1.first;
	int mm =  mu1.second;
	mup = nn * (nn + 1) + mm - 1;
        		
 	for (nu1 = 0; nu1 < nuMax ; nu1++){
	
	std::complex<double> I11(0.0,0.0), I12(0.0,0.0), I21(0.0,0.0), I22(0.0,0.0);
        std::complex<double> I31(0.0,0.0), I32(0.0,0.0), I41(0.0,0.0), I42(0.0,0.0);
	
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
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I21 = I21 + wi * det *(k_b_SH)*resDOT;
	///***********
        Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I22 = I22 + wi * det *(k_s_SH)*resDOT;
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I31 = I31 + wi * det *(k_b_SH)*resDOT;
	///***********
        Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I32 = I32 + wi * det *(k_s_SH)*resDOT;
	
	//*****************************************
	Tools::cross(resCR, aCoefint1.N(static_cast<long>(nu1)), aCoefext3.M(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I41 = I41 + wi * det *(k_b_SH)*resDOT;
	///***********
	Tools::cross(resCR, aCoefint1.M(static_cast<long>(nu1)), aCoefext3.N(static_cast<long>(mup)));
	
	resDOT = Tools::dot (&nvec[0], resCR);
	
	I42 = I42 + wi * det *(k_s_SH)*resDOT;
	//*****************************************
	
	
	}// Gauss Legendre integration of spherical harmonics over a triangle 
        }// iteration over all triangles, complete surface integration	

        QmatrixSH_a (brojac) = (I11 + I12);
        QmatrixSH_b (brojac) = (I21 + I22);
        QmatrixSH_c (brojac) = (I31 + I32);
        QmatrixSH_d (brojac) = (I41 + I42);
 
        brojac++;

	}
	
	}// iteration over all spherical harmonics

      QmatrixSH_1.head(size) = QmatrixSH_a;
      QmatrixSH_1.tail(size) = QmatrixSH_b;

      QmatrixSH_2.head(size) = QmatrixSH_c;
      QmatrixSH_2.tail(size) = QmatrixSH_d;

      QmatrixSH.head(2*size) = QmatrixSH_1;
      QmatrixSH.tail(2*size) = QmatrixSH_2;
  
  }


void Scatterer::getRgQLocalSH(optimet::Vector<optimet::t_complex>& RgQmatrixSH, optimet::t_real omega_, ElectroMagnetic const &bground, int gran1, int gran2) const {
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
          int brojac(0);
          int size = nuMax *(gran2 - gran1);
   
          std::complex<double>* resCR = new std::complex<double>[3];
          std::complex<double> resDOT;

         Vector<t_complex> RgQmatrixSH_a(size), RgQmatrixSH_b(size), RgQmatrixSH_c(size), RgQmatrixSH_d(size);
         Vector<t_complex> RgQmatrixSH_1(2*size), RgQmatrixSH_2(2*size);
          
         auto const k_s_SH = 2 * omega_ * std::sqrt(elmag.epsilon_SH * elmag.mu_SH); 
         auto const k_b_SH = 2 * omega_ * std::sqrt(bground.epsilon * bground.mu);
  
        unsigned int Nt = (topol.size()) / 3;  //number of triangles
        
        std::vector<std::vector<double>> Points = Tools::getPoints4();

	std::vector<double> Weights = Tools::getWeights4();

	for (mu1 = gran1; mu1 < gran2 ; mu1++){
	
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

        RgQmatrixSH_a (brojac) = (RgI11 + RgI12);
	RgQmatrixSH_b (brojac) = (RgI21 + RgI22);
	RgQmatrixSH_c (brojac) = (RgI31 + RgI32);
	RgQmatrixSH_d (brojac) = (RgI41 + RgI42);
        
        brojac++; 
	}

	}// iteration over all spherical harmonics

      RgQmatrixSH_1.head(size) = RgQmatrixSH_a;
      RgQmatrixSH_1.tail(size) = RgQmatrixSH_b;

      RgQmatrixSH_2.head(size) = RgQmatrixSH_c;
      RgQmatrixSH_2.tail(size) = RgQmatrixSH_d;

      RgQmatrixSH.head(2*size) = RgQmatrixSH_1;
      RgQmatrixSH.tail(2*size) = RgQmatrixSH_2;

     }
 #endif                         
    
