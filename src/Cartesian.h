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

#ifndef CARTESIAN_H_
#define CARTESIAN_H_

/**
 * The Cartesian class implements Cartesian coordinates.
 * This is a template class so any standard type may be used to create
 * spherical coordinates.
 * @warning Do not use without initialization.
 */
template <class carType> class Cartesian {
public:
  carType x; /**< The x cartesian coordinate. */
  carType y; /**< The y cartesian coordinate. */
  carType z; /**< The z cartesian coordinate. */

  /**
   * Initialization constructor for Cartesian.
   * @param x_ the carType value for x
   * @param y_ the carType value for y
   * @param z_ the carType value for z
   * @see init()
   */
  Cartesian(carType x_, carType y_, carType z_) { init(x_, y_, z_); }

  /**
   * Default Cartesian constructor.
   */
  Cartesian(void) {
    //
  }

  /**
   * Default Cartesian destructor.
   */
  ~Cartesian(void) {
    //
  }

  /**
   * Initialize a Cartesian object.
   * @param x_ the carType value for x
   * @param y_ the carType value for y
   * @param z_ the carType value for z
   * @see Cartesian()
   */
  void init(carType x_, carType y_, carType z_) {
    x = x_;
    y = y_;
    z = z_;
  }

  /**
   * Implements the dot-product of two cartesian vectors.
   * @param argument_ the vector to be multiplied with.
   * @return the dot-product of this and argument_.
   */
  carType operator*(Cartesian<carType> argument_) {
    return x * argument_.x + y * argument_.y + z * argument_.z;
  }

  /**
   * Implements the difference between two cartesian vectors.
   * @param argument_ the vector to be subtracted.
   * @return
   */
  Cartesian<carType> operator-(Cartesian<carType> argument_) {
    return Cartesian<carType>(x - argument_.x, y - argument_.y,
                              z - argument_.z);
  }
};

#endif /* CARTESIAN_H_ */
