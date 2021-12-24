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

#ifndef SPHERICAL_H_
#define SPHERICAL_H_

#include "Cartesian.h"
#include "SphericalP.h"
#include "Types.h"
#include <Eigen/Core>
#include <cmath>

/**
 * The Spherical class implements spherical coordinates.
 * This is a template class so any standard type may be used to create
 * spherical coordinates.
 * @warning Do not use without initialization.
 */
template <class T> class Spherical {
public:
  /**
   * Returns the point in Cartesian format.
   * @return the Cartesian coordinate object.
   */
  Cartesian<T> toCartesian() const {
    return Cartesian<T>(rrr * std::sin(the) * std::cos(phi), rrr * std::sin(the) * std::sin(phi),
                        rrr * std::cos(the));
  }

  //! \brief Transforms to cartesian vector, but not home-brewed version
  Eigen::Matrix<optimet::t_real, 3, 1> toEigenCartesian() const {
    return Eigen::Matrix<optimet::t_real, 3, 1>(rrr * std::sin(the) * std::cos(phi),
                                                rrr * std::sin(the) * std::sin(phi),
                                                rrr * std::cos(the));
  }

  /**
   * Returns the point in spherical projection format.
   * @return the SphericalP coordinate object.
   */
  SphericalP<T> toSphericalP() const {
    return SphericalP<T>(rrr * std::sin(the) * std::cos(phi), rrr * std::sin(the) * std::sin(phi),
                         rrr * std::cos(the));
  }

  /**
   * Converts a Cartesian vector into a Spherical vector.
   * @param point the Cartesian vector to be converted.
   * @return the vector in Spherical format.
   */
  static Spherical<T> toSpherical(Cartesian<T> point) {
    double r_l = std::sqrt(point.x * point.x + point.y * point.y + point.z * point.z);
    if(r_l > 0.0) {
      return Spherical<double>(r_l, std::acos(point.z / r_l), std::atan2(point.y, point.x));
    }

    return Spherical<double>(0.0, 0.0, 0.0);
  }
  template <class DERIVED> static Spherical<T> toSpherical(Eigen::MatrixBase<DERIVED> const &x) {
    return toSpherical(Cartesian<T>(x(0), x(1), x(2)));
  }

  T rrr; /**< The R spherical coordinate. */
  T the; /**< The Theta spherical coordinate. */
  T phi; /**< The Phi spherical coordinates. */

  /**
   * Initialization constructor for Spherical.
   * @param rrr_ the T value for R
   * @param the_ the T value for Theta
   * @param phi_ the T value for Phi
   * @see init()
   */
  Spherical(T rrr, T the, T phi) : rrr(rrr), the(the), phi(phi) {}

  /**
   * Default Spherical constructor.
   */
  Spherical() : Spherical(0, 0, 0) {}

  /**
   * Implements the dot-product of two spherical vectors.
   * @param argument_ the vector to be multiplied with.
   * @return the dot-product of this and argument_.
   */
  T operator*(Spherical<T> argument_) { return toCartesian() * argument_.toCartesian(); }

  /**
   * Implements the dot-product of a spherical and cartesian vector.
   * @param argument_ the cartesian vector to be multiplied with.
   * @return
   */
  T operator*(Cartesian<T> argument_) { return toCartesian() * argument_; }

  /**
   * Implements the dot-product of a spherical and spherical projection vector.
   * @param argument_ the SphericalP vector to be multiplied with.
   * @return
   */
  T operator*(SphericalP<T> argument_) { return toSphericalP() * argument_; }

  /**
   * Implements the vector*scalar product.
   * @param argument_ the scalar to be multiplied by.
   * @return the vector*scalar product.
   */
  Spherical<T> operator*(T argument_) {
    return Spherical<T>(rrr * argument_, the * argument_, phi * argument_);
  }

  /**
   * Implements the difference between this vector and the argument_.
   * @param argument_ the vector to be subtracted.
   * @return a Spherical vector as the difference.
   */
  Spherical<T> operator-(Spherical<T> argument_) const {
    return toSpherical(toCartesian() - argument_.toCartesian());
  }
};

#endif /* SPHERICAL_H_ */
