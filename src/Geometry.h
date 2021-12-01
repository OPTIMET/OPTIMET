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

#ifndef GEOMETRY_H_
#define GEOMETRY_H_

#include "Excitation.h"
#include "Scatterer.h"
#include "Types.h"
#include <memory>
#include <numeric>
#include <vector>

/**
 * The Geometry class implements a list of objects and properties of the medium.
 * To use the class, first initialize with init or the constructor by
 * declaring an expected capacity. Add objects with pushObject (Scatterer type)
 * until geometry is full. Finally, always validate to make sure the number
 * of objects are properly stored and objects do not intersect.
 * @warning Never use this class without initializing and validating.
 */
class Geometry {
private:
  //
public:
  std::vector<Scatterer> objects; /**< The list of scatterers. */

  ElectroMagnetic bground; /**< The properties of the background. */
  
  bool ACA_cond_;
  /**
   * Default constructor for the Geometry class. Does not initialize.
   */
  Geometry();

  /**
   * Default destructor for the Geometry class.
   */
  virtual ~Geometry();

  /**
   * Initialize the background. Default is vacuum.
   * @param bground_ the ElectroMagnetic properties of the background.
   */
  void initBground(ElectroMagnetic bground_);

  /**
   * Add an object to the Geometry.
   * @param object_ the Scatterer type object to be added.
   * @return 0 if add successful, 1 otherwise
   */
  void pushObject(Scatterer const &object_);
  
  // conditions for ACA compression
  void ACAcompression(bool ACA_cond){ACA_cond_ = ACA_cond;}
  bool get_ACAcond()const{return ACA_cond_;}

  //! \brief Validate geometry
  //! \details Fails if no objects, or if two objects overlap.
  bool is_valid() const;
  
  // vector used in the transformation of outer to inner expansion coefficients Analytical results FF
  int getCabsAux(double omega_, int objectIndex_, int nMax_, double *Cabs_aux_);

  /**
   * Checks if a point R is located inside an object.
   * @param R_ coordinates of the point to check.
   * @return the index of the object in which this point is or -1 if outside.
   */
  int checkInner(Spherical<double> R_);
  // Clebsch Gordan series coeff
  void Coefficients(int nMax, int nMaxS, std::vector<double *> CLGcoeff, int gran1, int gran2);

   // Incident coefficients for the second harmonic case                  
  int getIncLocalSH(std::vector<double *> CLGcoeff, int objectIndex_, std::shared_ptr<optimet::Excitation const> incWave_,
         optimet::Vector<optimet::t_complex> &internalCoef_FF_, int nMaxS_, std::complex<double> *Inc_local); 
  
  // Incident coefficients for the second harmonic case for parallelization                 
  int getIncLocalSH_parallel(std::vector<double *> CLGcoeff, int gran1, int gran2, std::shared_ptr<optimet::Excitation const> incWave_,
         optimet::Vector<optimet::t_complex> &internalCoef_FF_, int nMaxS_, std::complex<double> *Inc_local);     
           
 // Coefficient for absorption cross section second harmonic      
  int AbsCSSHcoeff(std::vector<double *> CLGcoeff, int gran1, int gran2, std::shared_ptr<optimet::Excitation const> incWave_, 
            optimet::Vector<optimet::t_complex> &internalCoef_FF_, optimet::Vector<optimet::t_complex> &internalCoef_SH_, int 
            nMaxS_, std::complex<double> *coefABS); 
            
   // Coefficients for the particular solution of SH differential equations                                     
  int COEFFpartSH(int objectIndex_, std::shared_ptr<optimet::Excitation const> incWave_, optimet::Vector<optimet::t_complex> 
                  &internalCoef_FF_, double r, int nMaxS_, std::complex<double> *coefXmn, std::complex<double> *coefXpl, 
                   std::vector<double *> CLGcoeff);

 // vectors needed for the SH arbitrary shapes
  void getEXCvecSH_ARB3(optimet::Vector<optimet::t_complex>& EXvec, std::shared_ptr<optimet::Excitation const> excitation, optimet::Vector<optimet::t_complex> &externalCoef_FF_, int objIndex);       
  
 void getEXCvecSH_ARB1(optimet::Vector<optimet::t_complex>& EXvec, std::shared_ptr<optimet::Excitation const> excitation, optimet::Vector<optimet::t_complex> &externalCoef_FF_, int objIndex);                                            
// vectors needed for the SH arbitrary shape in parallel
void getEXCvecSH_ARB3_parall(optimet::Vector<optimet::t_complex>& EXvec, std::shared_ptr<optimet::Excitation const> excitation, optimet::Vector<optimet::t_complex> &externalCoef_FF_, int gran1, int gran2);

void getEXCvecSH_ARB1_parall(optimet::Vector<optimet::t_complex>& EXvec, std::shared_ptr<optimet::Excitation const> excitation, optimet::Vector<optimet::t_complex> &externalCoef_FF_, int gran1, int gran2);

  /**
   * Updates the Geometry object to a new Excitation.
   * @param lambda_ the new wavelength.
   */
  void update(std::shared_ptr<optimet::Excitation const> incWave_);

  //! Size of the scattering vector
  optimet::t_uint scatterer_size() const;

  //! Maximum nMax across objects
  optimet::t_uint nMax() const {
    return std::accumulate(objects.begin(), objects.end(), optimet::t_uint(0u),
                           [](optimet::t_uint prior, Scatterer const &current) {
                             return std::max<optimet::t_uint>(prior, current.nMax);
                           });
  }

 //! Maximum nMaxS across objects
  optimet::t_uint nMaxS() const {
    return std::accumulate(objects.begin(), objects.end(), optimet::t_uint(0u),
                           [](optimet::t_uint prior, Scatterer const &current) {
                             return std::max<optimet::t_uint>(prior, current.nMaxS);
                           });
  }

  //! Minimum nMax across objects
  optimet::t_uint nMin() const {
    return std::accumulate(objects.begin(), objects.end(), nMax(),
                           [](optimet::t_uint prior, Scatterer const &current) {
                             return std::min<optimet::t_uint>(prior, current.nMax);
                           });
  }

protected:
  //! Validate last added sphere
  bool no_overlap(Scatterer const &object);
};

#endif /* GEOMETRY_H_ */
