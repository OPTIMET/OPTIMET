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

#ifndef SCATTERER_H_
#define SCATTERER_H_

#include "ElectroMagnetic.h"
#include "Spherical.h"
#include "Types.h"
#include "CompoundIterator.h"
#include <vector>
#include <tuple>
#include <cstring>


/**
 * The Scatterer class is the highest level element of a geometry.
 * This is the only class used to create the geometry. The radius property
 * will define a virtual sphere encompassing the entire structure. Different
 * behaviors will be created based on the type of scatterer. 
 */
class Scatterer {
private:

        std::vector<double> coord;
        std::vector<int> topol;
       
public:
        void Mesh(std::vector<double> co, std::vector<int> top);
	
	const double* getCoord(int vertex_number)const{
	return &(coord[vertex_number * 3]);
	}
	
	 const int* getNOvertex(int trian_number) const {
	
	return &(topol[trian_number*3]);
        }

  const int getTopolsize() const {
	
	return topol.size();
        }
  /**
   * Default Scatterer constructor.
   * Does NOT initialize the object.
   */
  Scatterer(int nMax_, int nMaxS_);

  /**
   * Initialization constructor for Scatterer.
   * @param vR_ the coordinates of the center of the scatterer.
   * @param elmag_ the electromagnetic properties of the scatterer.
   * @param radius_ the radius of the virtual sphere.
   * @param nMax_ the maximum value of the n iterator.
   * @see init()
   */
  Scatterer(Spherical<double> vR_, ElectroMagnetic elmag_, double radius_, int nMax_, int nMaxS_);

  /**
   * Initialization constructor for Scatterer.
   * @param pos the position in cartesian coordinates
   * @param elmag_ the electromagnetic properties of the scatterer.
   * @param radius_ the radius of the virtual sphere.
   * @param nMax_ the maximum value of the n iterator.
   */
  template <class T>
  Scatterer(Eigen::MatrixBase<T> const &pos, ElectroMagnetic elmag, double radius, int nMax)
      : Scatterer(Spherical<optimet::t_real>::toSpherical(
                      Cartesian<optimet::t_real>(pos[0], pos[1], pos[2])),
                  elmag, radius, nMax) {}

  /**
   * Default scatterer destructor.
   */
  virtual ~Scatterer();

  Spherical<double> vR;  /**< The coordinates of the center of the scatterer.*/
  ElectroMagnetic elmag; /**< The electromagnetic properties of the scatterer.*/
  double radius;         /**< The radius of a sphere encompassing the scatterer.*/
  double side;           // the length of the side of Platonic solid used in meshes
  int nMax;              /**< Maximum value of the n iterator. */
  int nMaxS;              /**< Maximum value of the n iterator SH */  
  std::string scatterer_type; // type of scatterer, sphere or arbitrary shaped
  #ifdef OPTIMET_MPI
  // FF Q matrix
  void getQLocal(optimet::Vector<optimet::t_complex>& Qmatrix, optimet::t_real omega_, ElectroMagnetic const &bground, int gran1, int gran2) const;

  // FF RgQmatrix
  void getRgQLocal(optimet::Vector<optimet::t_complex>& RgQmatrix, optimet::t_real omega_, ElectroMagnetic const &bground, int gran1, int gran2) const;  
  // SH Qmatrix for arbitrary shaped objects
  void getQLocalSH(optimet::Vector<optimet::t_complex>& QmatrixSH, optimet::t_real omega_, ElectroMagnetic const &bground, int gran1, int gran2) const;
  
  // SH RgQmatrix for arbitrary shaped objects
  void getRgQLocalSH(optimet::Vector<optimet::t_complex>& RgQmatrixSH, optimet::t_real omega_, ElectroMagnetic const &bground, int gran1, int gran2) const;
#endif  
};

#endif /* SCATTERER_H_ */
