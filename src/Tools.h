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

#ifndef TOOLS_H_
#define TOOLS_H_

#include "Types.h"
#include "Spherical.h"
#include "Cartesian.h"
#include <complex>
#include <vector>


/**
 * The Tools class implements several static tools for common use.
 */
class Tools {
public:
  /**
   * Default constructor for the Tools class.
   */
  Tools();

  /**
   * Default destructor for the Tools class.
   */
  virtual ~Tools();

  /**
   * Calculates the distance between two points in spherical coordinates.
   * @param point1 the coordinates of point1 in Spherical<double>.
   * @param point2 the coordinates of point2 in Spherical<double>.
   * @return the distance between the two points.
   */
  static double findDistance(Spherical<double> point1,
                             Spherical<double> point2);

  /**
   * Calculates the distance between two points in cartesian coordinates.
   * @param point1 the coordinates of point1 in Cartesian<double>.
   * @param point2 the coordinates of point2 in Cartesian<double>.
   * @return the distance between the two points.
   */
  static double findDistance(Cartesian<double> point1,
                             Cartesian<double> point2);

  /**
   * Projects a spherical point onto a SphericalP vector (cartesian projection).
   * @param point the point to be projected.
   * @param vector the vector to be used as basis for projection.
   * @return the projected vector in SphericalP.
   */
  static SphericalP<std::complex<double>>
  toProjection(Spherical<double> point,
               SphericalP<std::complex<double>> vector);

  /**
   * Projects a spherical point onto a SphericalP vector (spherical projection).
   * @param point the point to be projected.
   * @param vector the vector to be used as basis for projection.
   * @return the projected vector in SphericalP.
   */
  static SphericalP<std::complex<double>>
  fromProjection(Spherical<double> point,
                 SphericalP<std::complex<double>> vector);

  /**
   * Converts a point from spherical to SphericalP coordinates.
   * @param point - the coordinates of the point in
   * Spherical<std::complex<double> >.
   * @return the coordinates of the point in SphericalP<std::complex<double> >.
   */
  static SphericalP<std::complex<double>>
  toSphericalP(Spherical<std::complex<double>> point);

  /**
   * Converts a point from spherical to SphericalP coordinates.
   * @param point - the coordinates of the point in Spherical<double>.
   * @return the coordinates of the point in SphericalP<double>.
   */
  static SphericalP<double> toSphericalP(Spherical<double> point);

  /**
   * Converts a point from spherical to cartesian coordinates.
   * @param point - the coordinates of the point in Spherical<double>.
   * @return the coordinates of the point in Cartesian<double>.
   */
  static Cartesian<double> toCartesian(Spherical<double> point);

  /**
   * Converts a point from cartesian to spherical coordinates.
   * @param point - the coordinates of the point in Cartesian<double>.
   * @return the coordinates of the point in Spherical<double>.
   */
  static Spherical<double> toSpherical(Cartesian<double> point);

  /**
   * Calculates the maximum value of a compound iterator with n.
   * @param n the maximum value of n.
   * @return the maximum value of a compound iterator with n.
   */
  static long iteratorMax(long n);

  // scalar product between two vectors (double)
   static double dot(double *uvec,  double *vvec);
  
  // dot product (complex)
   static std::complex<double> dot(double* uvec,  std::complex<double>* vvec);

   // dot product SphericalP
   static std::complex<double> dot(double* uvec,  SphericalP<std::complex<double>> vvec);

    // dot product (cmplex double, SphericalP)
    static std::complex<double> dot(std::complex<double>* uvec,  SphericalP<std::complex<double>> vvec);

   // cross product of two vectors (double) 
    static void cross(double* res, double* u, double* v);

  // cross product of two vectors (SphericalP) 
  static void cross(std::complex<double>* res, SphericalP<std::complex<double>> u, SphericalP<std::complex<double>> v);

  // cross product of two vectors (double and SphericalP) 
   static void cross(std::complex<double>* res, double* u, SphericalP<std::complex<double>> v);

   // tangential component of a vector (-nxnxE)
   static void crossTanTr(std::complex<double>* res, std::complex<double>* resCR, double* u, SphericalP<std::complex<double>> v);

   // minus tangential component of a vector (nxnx)
   static void crossTan(std::complex<double>* res, std::complex<double>* resCR, double* u, std::complex<double>* v);

  // different norms
        static double norm2(double* vec);
  
        static double norm(double* vec);

   // Gauss-Legendre integration points on a triangle 
   static std::vector<std::vector<double>> getPoints4();
   
   static std::vector<double> getWeights4();
   
    static std::vector<std::vector<double>> getPoints7();
   
   static std::vector<double> getWeights7();

   // Gauss-Legendre integration points on a line
       
    static std::vector<double> getLineWghts4();
             
    static std::vector<double> getLinePts4();
                   
    static std::vector<double> getLineWghts6();
                         
    static std::vector<double> getLinePts6();
 
  /**
   * Translate from one coordinate system to another.
   * @param R the point to be translated
   * @param P the center of the new coordinate system.
   * @return the new coordinate of point R in center P.
   */
  static Spherical<double> toPoint(Spherical<double> R, Spherical<double> P);
};

#endif /* TOOLS_H_ */
