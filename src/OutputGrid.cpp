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

#include "OutputGrid.h"

#include "Aliases.h"
#include "Cartesian.h"
#include "Tools.h"
#include <cmath>

namespace optimet {
OutputGrid::OutputGrid() : initDone(false), gridDone(false), iterator(0), gridPoints(0) {}

OutputGrid::~OutputGrid() {
  // Destructor does nothing. Close
}

OutputGrid::OutputGrid(int type_, std::array<t_real, 9> const &parameters_, hid_t groupID_)
    : OutputGrid() {
  init(type_, parameters_, groupID_);
}

OutputGrid::OutputGrid(int type_, std::array<t_real, 9> const &parameters_)
    : OutputGrid() {
  init(type_, parameters_);
}

void OutputGrid::init(int type_, std::array<t_real, 9> const &parameters_, hid_t groupID_) {
  type = type_;
  gridParameters = parameters_;
  groupID = groupID_;

  // Cartesian Regular
  if(type == O3DCartesianRegular) {

    iterator = 0;
    gridPoints = (int)(gridParameters[2] * gridParameters[5] * gridParameters[8]);

    cursor = {{0, 0, 0}};
    // Step sizes (X, Y, Z)
    aux = {{std::abs(gridParameters[1] - gridParameters[0]) / (gridParameters[2] - 1),
            std::abs(gridParameters[4] - gridParameters[3]) / (gridParameters[5] - 1),
            std::abs(gridParameters[7] - gridParameters[6]) / (gridParameters[8] - 1)}};

    // Create the HDF5 SubGroups
    vecGroupId[0] = H5Gcreate(groupID, "X", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    vecGroupId[1] = H5Gcreate(groupID, "Y", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    vecGroupId[2] = H5Gcreate(groupID, "Z", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    vecGroupId[3] = H5Gcreate(groupID, "ABS", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // Create the associated datasets and dataspaces
    hsize_t dims[3];

    dims[0] = gridParameters[2];
    dims[1] = gridParameters[5];
    dims[2] = gridParameters[8];

    // Create the dataspaces for real and imag
    for(int i = 0; i < 7; i++) {
      vecDSpaceId[i] = H5Screate_simple(3, dims, NULL);
    }

    // Create the associated datasets
    vecDataId[0] = H5Dcreate(vecGroupId[0], "real", H5T_NATIVE_DOUBLE, vecDSpaceId[0], H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    vecDataId[1] = H5Dcreate(vecGroupId[0], "imag", H5T_NATIVE_DOUBLE, vecDSpaceId[1], H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    vecDataId[2] = H5Dcreate(vecGroupId[1], "real", H5T_NATIVE_DOUBLE, vecDSpaceId[2], H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    vecDataId[3] = H5Dcreate(vecGroupId[1], "imag", H5T_NATIVE_DOUBLE, vecDSpaceId[3], H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    vecDataId[4] = H5Dcreate(vecGroupId[2], "real", H5T_NATIVE_DOUBLE, vecDSpaceId[4], H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    vecDataId[5] = H5Dcreate(vecGroupId[2], "imag", H5T_NATIVE_DOUBLE, vecDSpaceId[5], H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
    vecDataId[6] = H5Dcreate(vecGroupId[3], "abs", H5T_NATIVE_DOUBLE, vecDSpaceId[6], H5P_DEFAULT,
                             H5P_DEFAULT, H5P_DEFAULT);
  }

  gridDone = false;
  initDone = true;
}


void OutputGrid::init(int type_, std::array<t_real, 9> const &parameters_) {
  type = type_;
  gridParameters = parameters_;
  

  // Cartesian Regular
  if(type == O3DCartesianRegular) {

    iterator = 0;
    gridPoints = (int)(gridParameters[2] * gridParameters[5] * gridParameters[8]);

    cursor = {{0, 0, 0}};
    // Step sizes (X, Y, Z)
    aux = {{std::abs(gridParameters[1] - gridParameters[0]) / (gridParameters[2] - 1),
            std::abs(gridParameters[4] - gridParameters[3]) / (gridParameters[5] - 1),
            std::abs(gridParameters[7] - gridParameters[6]) / (gridParameters[8] - 1)}};
}
gridDone = false;
  initDone = true;

}


void OutputGrid::gotoStart() {
  iterator = 0;
  gridDone = false;
}

void OutputGrid::gotoNext() {
  if(gridDone)
    iterator = 0;
  else
    iterator++;
}

Spherical<double> OutputGrid::getPoint() {
  if(gridDone)
    iterator = 0;

  if(iterator == gridPoints - 1) {
    gridDone = true;
  }

  if(type == O3DCartesianRegular) // Regular Cartesian grid
  {
    // Get the cursor coordinates from the iterator.
    cursor[0] = iterator % ((int)gridParameters[2]);
    cursor[1] = (iterator / ((int)gridParameters[2])) % ((int)gridParameters[5]);
    cursor[2] = iterator / (((int)gridParameters[2]) * ((int)gridParameters[5]));

    double local_x = gridParameters[0] + cursor[0] * aux[0] + 1e-12; // Correct for 0
    double local_y = gridParameters[3] + cursor[1] * aux[1] + 1e-12; // Correct for 0
    double local_z = gridParameters[6] + cursor[2] * aux[2] + 1e-12; // Correct for 0
    // diagonal plane just for zincblende local_z=local_y
    // Create a spherical vector from Cartesian and return it with the local
    // coordinates
    return Tools::toSpherical(Cartesian<double>(local_x, local_y, local_z));
  }

  // Default return is 0
  return Spherical<double>(0.0, 0.0, 0.0);
}

void OutputGrid::pushData(SphericalP<std::complex<double>> data_) {
  if(type == O3DCartesianRegular) // Cartesian Regular grid
  {
    hsize_t dim2[] = {1}; // Dimension size of the auxiliary dataset (in memory)
    hsize_t coord[1][3];  // Array to store selected points from the file dataspace

    // Create dataspace for the auxiliary dataset.
    hid_t mid2 = H5Screate_simple(1, dim2, NULL);

    // Set cursor position
    coord[0][0] = cursor[0];
    coord[0][1] = cursor[1];
    coord[0][2] = cursor[2];

    // Select sequence of points in the file dataspace(fid).
    // Write new selection of points to the dataset(dataset).
    // Repeat for all six components.
    double real, imag, data_abs;
    real = data_.rrr.real();

     data_abs = std::sqrt( std::pow((std::abs(data_.rrr)), 2) + std::pow((std::abs(data_.the)), 2)                    + std::pow((std::abs(data_.phi)), 2) );

    H5Sselect_elements(vecDSpaceId[0], H5S_SELECT_SET, 1, (const hsize_t *)coord);
    H5Dwrite(vecDataId[0], H5T_NATIVE_DOUBLE, mid2, vecDSpaceId[0], H5P_DEFAULT, &real);

    imag = data_.rrr.imag();
    H5Sselect_elements(vecDSpaceId[1], H5S_SELECT_SET, 1, (const hsize_t *)coord);
    H5Dwrite(vecDataId[1], H5T_NATIVE_DOUBLE, mid2, vecDSpaceId[1], H5P_DEFAULT, &imag);

    real = data_.the.real();
    H5Sselect_elements(vecDSpaceId[2], H5S_SELECT_SET, 1, (const hsize_t *)coord);
    H5Dwrite(vecDataId[2], H5T_NATIVE_DOUBLE, mid2, vecDSpaceId[2], H5P_DEFAULT, &real);

    imag = data_.the.imag();
    H5Sselect_elements(vecDSpaceId[3], H5S_SELECT_SET, 1, (const hsize_t *)coord);
    H5Dwrite(vecDataId[3], H5T_NATIVE_DOUBLE, mid2, vecDSpaceId[3], H5P_DEFAULT, &imag);

    real = data_.phi.real();
    H5Sselect_elements(vecDSpaceId[4], H5S_SELECT_SET, 1, (const hsize_t *)coord);
    H5Dwrite(vecDataId[4], H5T_NATIVE_DOUBLE, mid2, vecDSpaceId[4], H5P_DEFAULT, &real);

    imag = data_.phi.imag();
    H5Sselect_elements(vecDSpaceId[5], H5S_SELECT_SET, 1, (const hsize_t *)coord);
    H5Dwrite(vecDataId[5], H5T_NATIVE_DOUBLE, mid2, vecDSpaceId[5], H5P_DEFAULT, &imag);
    
    H5Sselect_elements(vecDSpaceId[6], H5S_SELECT_SET, 1, (const hsize_t *)coord);
    H5Dwrite(vecDataId[6], H5T_NATIVE_DOUBLE, mid2, vecDSpaceId[6], H5P_DEFAULT, &data_abs);

    // Close the auxiliary dataspace
    H5Sclose(mid2);
  }
}

void OutputGrid::pushDataNext(SphericalP<std::complex<double>> data_) {
  pushData(data_);
  gotoNext();
}

void OutputGrid::close() {
  if(initDone) {
    // Close datasets
    for(int i = 0; i < 7; i++)
      H5Dclose(vecDataId[i]);

    // Close dataspaces
    for(int i = 0; i < 7; i++)
      H5Sclose(vecDSpaceId[i]);

    // Close sub-groups
    for(int i = 0; i < 4; i++)
      H5Gclose(vecGroupId[i]);

    // Close the Group from Output
    H5Gclose(groupID);
  }
}
}
