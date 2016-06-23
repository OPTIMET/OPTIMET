#ifndef OUTPUT_GRID_H_
#define OUTPUT_GRID_H_

#include "Types.h"
#include "Spherical.h"
#include "SphericalP.h"
#include <complex>
#include <hdf5.h>

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
 *
 * At the moment, output will be E fields for all values.
 */
class OutputGrid {
private:
  hid_t vecGroupId[3]; /**< The Group IDs for each of the three subgroups.
                          Internal. */
  hid_t
      vecDataId[6]; /**< The Data IDs for each of the six datasets. Internal. */
  hid_t vecDSpaceId[6]; /**< The Data IDs for each of the six dataspaces.
                           Internal. */
  double *aux;          /**< Auxiliary data used internally. */
  int *cursor;   /**< Cursor in the current grid mapped to the HDF5 dataset. */
  bool initDone; /**< Specifies the initialization state of the object. */
public:
  bool gridDone;  /**< Specifies if the grid has been fully parsed. */
  int iterator;   /**< The current position of the iterator in the grid. */
  int type;       /**< The grid type. */
  int gridPoints; /**< The number of grid points. */
  hid_t groupID;  /**< ID of the HDF5 group. */
  double *gridParameters; /**< The gridParameters (size and values depend on
                             grid type. */

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
  OutputGrid(int type_, double *parameters_, hid_t groupID_);

  /**
   * Initialization method for the OutputGrid class.
   * @param type_ the grid type as defined in "Aliases.h".
   * @param parameters_ the parameters (depending on the grid type). See class
   * docs.
   * @param groupID_ the HDF5 data space.
   */
  void init(int type_, double *parameters_, hid_t groupID_);

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

#endif /* OUTPUT_GRID_H_ */
