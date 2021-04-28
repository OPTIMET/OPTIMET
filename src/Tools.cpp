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

#include "Tools.h"
#include <cmath>
#include <assert.h>
#include <iostream>

Tools::Tools() {
  //
}

Tools::~Tools() {
  //
}

double Tools::findDistance(Spherical<double> point1, Spherical<double> point2) {
  return Tools::findDistance(toCartesian(point1), toCartesian(point2));
}

double Tools::findDistance(Cartesian<double> point1, Cartesian<double> point2) {
  return std::sqrt(std::pow((point2.x - point1.x), 2.0) +
                   std::pow((point2.y - point1.y), 2.0) +
                   std::pow((point2.z - point1.z), 2.0));
}

Cartesian<double> Tools::toCartesian(Spherical<double> point) {
  return Cartesian<double>(point.rrr * sin(point.the) * std::cos(point.phi),
                           point.rrr * sin(point.the) * sin(point.phi),
                           point.rrr * std::cos(point.the));
}

// dot product (double)
double Tools::dot(double *uvec,  double *vvec)
   {
	   return uvec[0] * vvec[0] + uvec[1] * vvec[1] + uvec[2] * vvec[2];
   }
   
 // dot product (complex)
 std::complex<double> Tools::dot(double* uvec,  std::complex<double>* vvec)
   {
	   return uvec[0] * vvec[0] + uvec[1] * vvec[1] + uvec[2] * vvec[2];
   }

  // dot product (double, SphericalP)
  std::complex<double> Tools::dot(double* uvec,  SphericalP<std::complex<double>> vvec)
   {
	   return uvec[0] * vvec.rrr + uvec[1] * vvec.the + uvec[2] * vvec.phi;
   }
   
   // dot product (cmplex double, SphericalP)
   std::complex<double> Tools::dot(std::complex<double>* uvec,  SphericalP<std::complex<double>> vvec)
   {
	   return uvec[0] * vvec.rrr + uvec[1] * vvec.the + uvec[2] * vvec.phi;
   }
  
  // cross product of two vectors (double)
  void Tools::cross(double* res, double* u,
	 double* v)
{

	res[0] = u[1] * v[2] - u[2] *  v[1];
	res[1] = u[2] * v[0] - u[0] *  v[2];
	res[2] = u[0] * v[1] - u[1] *  v[0];
	
}

 void Tools::cross(std::complex<double>* res, SphericalP<std::complex<double>> u, SphericalP<std::complex<double>> v){
        
        res[0] = u.the * v.phi - u.phi *  v.the;
	res[1] = u.phi * v.rrr - u.rrr *  v.phi;
	res[2] = u.rrr * v.the - u.the *  v.rrr;

}

void Tools::cross(std::complex<double>* res, double* u, SphericalP<std::complex<double>> v){
        
        res[0] = u[1] * v.phi - u[2] *  v.the;
	res[1] = u[2] * v.rrr - u[0] *  v.phi;
	res[2] = u[0] * v.the - u[1] *  v.rrr;

}

// -nxnxE
void Tools::crossTanTr(std::complex<double>* res, std::complex<double>* resCR, double* u, SphericalP<std::complex<double>> v){
        
        resCR[0] = u[1] * v.phi - u[2] *  v.the;
	resCR[1] = u[2] * v.rrr - u[0] *  v.phi;
	resCR[2] = u[0] * v.the - u[1] *  v.rrr;
	
	res[0] = -(u[1] * resCR[2] - u[2] *  resCR[1]);
	res[1] = -(u[2] * resCR[0] - u[0] *  resCR[2]);
	res[2] = -(u[0] * resCR[1] - u[1] *  resCR[0]);

}

// nxnx
void Tools::crossTan(std::complex<double>* res, std::complex<double>* resCR, double* u, std::complex<double>* v){
        
        resCR[0] = u[1] * v[2] - u[2] *  v[1];
	resCR[1] = u[2] * v[0] - u[0] *  v[2];
	resCR[2] = u[0] * v[1] - u[1] *  v[0];
	
	res[0] = u[1] * resCR[2] - u[2] *  resCR[1];
	res[1] = u[2] * resCR[0] - u[0] *  resCR[2];
	res[2] = u[0] * resCR[1] - u[1] *  resCR[0];

}

// different norms
double Tools::norm2(double* vec)
   {
	   return dot(vec, vec);
   }

   
   double Tools::norm(double* vec)
   {
	   return sqrt(norm2(vec));
   }
   
   // Gauss-Legendre points on a triangle
 std::vector<double> Tools::getWeights4() {

	std::vector<double> w(4);
	w[0] = -2.812500000000000e-01;
	w[1] = w[2] = w[3] =2.604166666666670e-01;
	
	return w;

}

   std::vector<std::vector<double>> Tools::getPoints4() {

	std::vector<std::vector<double>>P;
	double intvalue = 0;
	P.resize(4, std::vector<double>(2, intvalue));

	P[0][0] = P[0][1] = (1.0 / 3.0);
	P[1][0] = P[2][1] = P[3][0] = P[3][1] = 0.2;
	P[1][1] = P[2][0] = 0.6;
	
	return P;

}

      
   std::vector<double> Tools::getWeights7() {

	std::vector<double> w(7);
	
	w[0] = 0.112500000000000;
	
	w[1] = w[2] = w[3] = 0.066197076394255;
	
	w[4] = w[5] = w[6] = 0.062969590272415;
	
	return w;

}

   std::vector<std::vector<double>> Tools::getPoints7() {

	std::vector<std::vector<double>>P;
	double intvalue = 0;
	P.resize(7, std::vector<double>(2, intvalue));

	P[0][0] = P[0][1] = (1.0 / 3.0);

	P[1][0] = P[1][1] = P[2][0] = P[3][1] = 0.47014206410511;
	
	P[2][1] = P[3][0] = 0.05971587178977;
	
	P[4][0]= P[4][1]= P[5][0]= P[6][1]= 0.10128650732346;
	
	P[5][1] = P[6][0] = 0.79742698535309;
	
	return P;

}

// Gauss Legendre points on a line
std::vector<double> Tools::getLineWghts4() {
	std::vector<double> w(4);
	w[0] = 0.6521451548625461;
	w[1] = 0.6521451548625461;
	w[2] = 0.3478548451374538;
	w[3] = 0.3478548451374538;
	
	return w;

}

std::vector<double> Tools::getLineWghts6() {
	std::vector<double> w(6);
	
	w[0] = 0.3607615730481386;
	w[1] = 0.3607615730481386;
	w[2] = 0.4679139345726910;
	w[3] = 0.4679139345726910;
	w[4] = 0.1713244923791704;
	w[5] = 0.1713244923791704;
	
	
	return w;

}

 std::vector<double> Tools::getLinePts4() {
	std::vector<double> P(4);

	P[0] = -0.3399810435848563;
	P[1] = 0.3399810435848563;
	P[2] = -0.8611363115940526;
	P[3] = 0.8611363115940526;
	
	return P;

}

std::vector<double> Tools::getLinePts6() {
	std::vector<double> P(6);

	P[0] = 0.6612093864662645;
	P[1] = -0.6612093864662645;
	P[2] = -0.2386191860831969;
	P[3] = 0.2386191860831969;
	P[4] = -0.9324695142031521;
	P[5] = 0.9324695142031521;
	
	return P;

}



Spherical<double> Tools::toSpherical(Cartesian<double> point) {
  double r_l =
      std::sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
  if (r_l > 0.0) {
    return Spherical<double>(r_l, std::acos(point.z / r_l),
                             std::atan2(point.y, point.x));
  }
  return Spherical<double>(0.0, 0.0, 0.0);
}

long Tools::iteratorMax(long n) { return n * n + 2 * n; }


SphericalP<double> Tools::toSphericalP(Spherical<double> point) {
  return SphericalP<double>(point.rrr * sin(point.the) * std::cos(point.phi),
                            point.rrr * sin(point.the) * sin(point.phi),
                            point.rrr * std::cos(point.the));
}

SphericalP<std::complex<double>>
Tools::toSphericalP(Spherical<std::complex<double>> point) {
  return SphericalP<std::complex<double>>(
      point.rrr * sin(point.the) * std::cos(point.phi),
      point.rrr * sin(point.the) * sin(point.phi),
      point.rrr * std::cos(point.the));
}

SphericalP<std::complex<double>>
Tools::toProjection(Spherical<double> point,
                    SphericalP<std::complex<double>> vector) {
  return SphericalP<std::complex<double>>(
      sin(point.the) * std::cos(point.phi) * vector.rrr +
          std::cos(point.the) * std::cos(point.phi) * vector.the -
          sin(point.phi) * vector.phi,
      sin(point.the) * sin(point.phi) * vector.rrr +
          std::cos(point.the) * sin(point.phi) * vector.the +
          std::cos(point.phi) * vector.phi,
      std::cos(point.the) * vector.rrr - sin(point.the) * vector.the);
}

SphericalP<std::complex<double>>
Tools::fromProjection(Spherical<double> point,
                      SphericalP<std::complex<double>> vector) {
  return SphericalP<std::complex<double>>(
      sin(point.the) * std::cos(point.phi) * vector.rrr +
          sin(point.the) * sin(point.phi) * vector.the +
          std::cos(point.the) * vector.phi,
      std::cos(point.the) * std::cos(point.phi) * vector.rrr +
          std::cos(point.the) * sin(point.phi) * vector.the -
          sin(point.the) * vector.phi,
      std::cos(point.phi) * vector.the - sin(point.phi) * vector.rrr);
}

Spherical<double> Tools::toPoint(Spherical<double> R, Spherical<double> P) {
  Cartesian<double> R_cart = toCartesian(R);
  Cartesian<double> P_cart = toCartesian(P);

  return toSpherical(R_cart - P_cart);
}

