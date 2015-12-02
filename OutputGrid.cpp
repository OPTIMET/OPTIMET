#include "OutputGrid.h"

#include "aliases.h"
#include "Cartesian.h"
#include "Tools.h"
#include <cmath>

using std::abs;

OutputGrid::OutputGrid()
{
  iterator = 0;
  gridPoints = 0;
  gridDone = false;
  initDone = false;
}

OutputGrid::~OutputGrid()
{
  //Destructor does nothing. Close
}

OutputGrid::OutputGrid(int type_, double *parameters_, hid_t groupID_)
{
  init(type_, parameters_, groupID_);
}

void OutputGrid::init(int type_, double *parameters_, hid_t groupID_)
{
  type = type_;
  gridParameters = parameters_;
  groupID = groupID_;

  //Cartesian Regular
  if(type == O3DCartesianRegular)
  {
    iterator = 0;
    gridPoints = (int)(gridParameters[2] * gridParameters[5] * gridParameters[8]);

    aux = new double[3];
    cursor = new int[3]; //cursor_x, cursor_y, cursor_z
    cursor[0] = 0;
    cursor[1] = 0;
    cursor[2] = 0;

    //Step sizes (X, Y, Z)
    aux[0] = abs(gridParameters[1] - gridParameters[0]) / (gridParameters[2] - 1);
    aux[1] = abs(gridParameters[4] - gridParameters[3]) / (gridParameters[5] - 1);
    aux[2] = abs(gridParameters[7] - gridParameters[6]) / (gridParameters[8] - 1);

    //Create the HDF5 SubGroups
     vecGroupId[0] = H5Gcreate(groupID, "X", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     vecGroupId[1] = H5Gcreate(groupID, "Y", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     vecGroupId[2] = H5Gcreate(groupID, "Z", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

     //Create the associated datasets and dataspaces
     hsize_t dims[3];

     dims[0] = gridParameters[2];
     dims[1] = gridParameters[5];
     dims[2] = gridParameters[8];

     //Create the dataspaces for real and imag
     for(int i=0; i<6; i++)
     {
       vecDSpaceId[i] = H5Screate_simple(3, dims, NULL);
     }

     //Create the associated datasets
     vecDataId[0] = H5Dcreate(vecGroupId[0], "real", H5T_NATIVE_DOUBLE, vecDSpaceId[0],
         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     vecDataId[1] = H5Dcreate(vecGroupId[0], "imag", H5T_NATIVE_DOUBLE, vecDSpaceId[1],
         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     vecDataId[2] = H5Dcreate(vecGroupId[1], "real", H5T_NATIVE_DOUBLE, vecDSpaceId[2],
         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     vecDataId[3] = H5Dcreate(vecGroupId[1], "imag", H5T_NATIVE_DOUBLE, vecDSpaceId[3],
         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     vecDataId[4] = H5Dcreate(vecGroupId[2], "real", H5T_NATIVE_DOUBLE, vecDSpaceId[4],
         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
     vecDataId[5] = H5Dcreate(vecGroupId[2], "imag", H5T_NATIVE_DOUBLE, vecDSpaceId[5],
         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  }

  gridDone = false;
  initDone = true;
}

void OutputGrid::gotoStart()
{
  iterator = 0;
  gridDone = false;
}

void OutputGrid::gotoNext()
{
  if(gridDone)
    iterator = 0;
  else
    iterator++;
}

Spherical<double> OutputGrid::getPoint()
{
  if (gridDone)
    iterator = 0;

  if(iterator == gridPoints-1)
  {
    gridDone = true;
  }

  if(type == O3DCartesianRegular) //Regular Cartesian grid
  {
    //Get the cursor coordinates from the iterator.
    cursor[0] = iterator % ((int)gridParameters[2]);
    cursor[1] = (iterator / ((int)gridParameters[2])) % ((int)gridParameters[5]);
    cursor[2] = iterator / ( ((int)gridParameters[2]) * ((int)gridParameters[5]));

    double local_x = gridParameters[0] + cursor[0] * aux[0] + 1e-12; //Correct for 0
    double local_y = gridParameters[3] + cursor[1] * aux[1] + 1e-12; //Correct for 0
    double local_z = gridParameters[6] + cursor[2] * aux[2] + 1e-12; //Correct for 0

    //Create a spherical vector from Cartesian and return it with the local coordinates
    return Tools::toSpherical(Cartesian<double>(local_x, local_y, local_z));
  }

  //Default return is 0
  return Spherical<double>(0.0, 0.0, 0.0);
}

void OutputGrid::pushData(SphericalP<complex<double> > data_)
{
  if(type == O3DCartesianRegular) //Cartesian Regular grid
  {
    hsize_t dim2[] = {1};   // Dimension size of the auxiliary dataset (in memory)
    hsize_t coord[1][3];      // Array to store selected points from the file dataspace

    // Create dataspace for the auxiliary dataset.
    hid_t mid2 = H5Screate_simple(1, dim2, NULL);

    // Set cursor position
    coord[0][0] = cursor[0];
    coord[0][1] = cursor[1];
    coord[0][2] = cursor[2];

    // Select sequence of points in the file dataspace(fid).
    // Write new selection of points to the dataset(dataset).
    // Repeat for all six components.
    double real, imag;
    real = data_.rrr.real();
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

    //Close the auxiliary dataspace
    H5Sclose(mid2);
  }
}

void OutputGrid::pushDataNext(SphericalP<complex<double> > data_)
{
  pushData(data_);
  gotoNext();
}

void OutputGrid::close()
{
  if(initDone)
  {
    //Close datasets
    for(int i=0; i<6; i++)
      H5Dclose(vecDataId[i]);

    //Close dataspaces
    for(int i=0; i<6; i++)
      H5Sclose(vecDSpaceId[i]);

    //Close sub-groups
    for(int i=0; i<3; i++)
      H5Gclose(vecGroupId[i]);

    //Close the Group from Output
    H5Gclose(groupID);

    delete [] aux;
    delete [] cursor;
  }
}
