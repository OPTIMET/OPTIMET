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

#ifndef OUTPUT_GRID_H_
#define OUTPUT_GRID_H_

#include "Spherical.h"
#include "SphericalP.h"
#include "Types.h"
#include <array>
#include <complex>
#include <hdf5.h>

namespace optimet {
/**
 * The OutputGrid class implements the field output grid and processes the HDF5
 * datasets.
 * The groupID will specify a group in which the sub-groups (r, theta, phi) and
 * (real, imag)
 * will be created.
 * The grid types and their associated parameters are:
 *
 * 1. O3DCartesianRegular (int value 911)
 *    A regular 3D cartesian grid with equal spacing between the points.
 *    Parameters: parameters[0] - the starting value for the X axis.
 *          parameters[1] - the end value for the X axis.
 *          parameters[2] - the number of points on the X axis.
 *          parameters[3] - the starting value for the Y axis.
 *          parameters[4] - the end value for the Y axis.
 *          parameters[5] - the number of points on the Y axis.
 *          parameters[6] - the starting value for the Z axis.
 *          parameters[7] - the end value for the Z axis.
 *          parameters[8] - the number of points on the Z axis.
 */
class OutputGrid {
private:
  hid_t vecGroupId[4];  /**< The Group IDs for each of the four subgroups.
                           Internal. */
  hid_t vecDataId[7];   /**< The Data IDs for each of the seven datasets. Internal. */
  hid_t vecDSpaceId[7]; /**< The Data IDs for each of the seven dataspaces.
                           Internal. */
  //! Auxiliary data used internally
  std::array<t_real, 3> aux;
  //! Cursor in the current grid mapped to the HDF5 dataset
  std::array<t_int, 3> cursor;
  bool initDone;        /**< Specifies the initialization state of the object. */
public:
  bool gridDone;  /**< Specifies if the grid has been fully parsed. */
  int iterator;   /**< The current position of the iterator in the grid. */
  int type;       /**< The grid type. */
  int gridPoints; /**< The number of grid points. */
  hid_t groupID;  /**< ID of the HDF5 group. */
  //! The gridParameters (size and values depend on grid type
  std::array<t_real, 9> gridParameters;

  /**
   * Default constructor for the OutputGrid class.
   * Does NOT initialize the object.
   */
  OutputGrid();

  /**
   * Default destructor for the OutputGrid class.
   */
  virtual ~OutputGrid();

  /**
   * Initialization constructor for the OutputGrid class.
   * @param type_ the grid type as defined in "Aliases.h".
   * @param parameters_ the parameters (depending on the grid type). See class
   * docs.
   * @param groupID_ the HDF5 data space.
   */
  OutputGrid(int type_, std::array<t_real, 9> const & parameters_, hid_t groupID_);

  OutputGrid(int type_, std::array<t_real, 9> const & parameters_);
  /**
   * Initialization method for the OutputGrid class.
   * @param type_ the grid type as defined in "Aliases.h".
   * @param parameters_ the parameters (depending on the grid type). See class
   * docs.
   * @param groupID_ the HDF5 data space.
   */
  void init(int type_, std::array<t_real, 9> const & parameters_, hid_t groupID_);

   void init(int type_, std::array<t_real, 9> const & parameters_);

  /**
   * Sets the iterator to the first point.
   */
  void gotoStart();

  /**
   * Moves the iterator to the next position.
   */
  void gotoNext();

  /**
   * Returns the current point of the grid in Spherical vector format.
   * If overflow, iterator is set back to start position.
   * @return the current grid point in Spherical<double> format.
   */
  Spherical<double> getPoint();

  /**
   * Push data to current iterator position in HDF5 file.
   * @param data_ the data in SphericalP<complex <double> > format.
   */
  void pushData(SphericalP<std::complex<double>> data_);

  /**
   * Push data to current iterator position in HDF5 file and move to next
   * positon.
   * @param data_ the data in SphericalP<complex <double> > format.
   */
  void pushDataNext(SphericalP<std::complex<double>> data_);

  /**
   * Close the grid to close all HDF5 structures.
   */
  void close();
};
}

#endif /* OUTPUT_GRID_H_ */
