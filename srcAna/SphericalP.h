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

#ifndef SPHERICALP_H_
#define SPHERICALP_H_

#include <complex>

/**
 * The SphericalP class implements spherical projection coordinates.
 * This is a template class so any standard type may be used to create
 * spherical coordinates.
 * @warning Do not use without initialization.
 */
template <class carType> class SphericalP {
public:
  carType rrr; /**< The uR coordinate. */
  carType the; /**< The uTheta coordinate. */
  carType phi; /**< The uPhi coordinate. */

  /**
   * Initialization constructor for SphericalP.
   * @param rrr_ the carType value for urrr
   * @param the_ the carType value for uthe
   * @param phi_ the carType value for uphi
   * @see init()
   */
  SphericalP(carType rrr_, carType the_, carType phi_) {
    init(rrr_, the_, phi_);
  }

  /**
   * Default SphericalP constructor.
   */
  SphericalP(void) {
    //
  }

  /**
   * Default SphericalP destructor.
   */
  ~SphericalP(void) {
    //
  }

  /**
   * Initialize a SphericalP object.
   * @param rrr_ the carType value for x
   * @param the_ the carType value for y
   * @param phi_ the carType value for z
   * @see Cartesian()
   */
  void init(carType rrr_, carType the_, carType phi_) {
    rrr = rrr_;
    the = the_;
    phi = phi_;
  }

  /**
   * Implements the dot-product of two spherical projection vectors.
   * @param argument_ the vector to be multiplied with.
   * @return the dot-product of this and argument_.
   */
  carType operator*(const SphericalP<carType> &argument) const {
    return rrr * argument.rrr + the * argument.the + phi * argument.phi;
  }

  SphericalP<carType> operator+(const SphericalP<carType> &argument) const {
    return SphericalP<carType>(rrr + argument.rrr, the + argument.the,
                               phi + argument.phi);
  }

  SphericalP<carType> operator-(const SphericalP<carType> &argument) const {
    return SphericalP<carType>(rrr - argument.rrr, the - argument.the,
                               phi - argument.phi);
  }

  SphericalP<carType> operator*(const carType &argument) const {
    return SphericalP<carType>(rrr * argument, the * argument, phi * argument);
  }
};

template <class carType>
SphericalP<carType> conj(const SphericalP<carType> &sp) {
  using std::conj;
  return SphericalP<carType>(conj(sp.rrr), conj(sp.the), conj(sp.phi));
}

#endif /* SPHERICALP_H_ */
