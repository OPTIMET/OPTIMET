#include "Output.h"

Output::Output()
{
  initDone = false;
}

Output::Output(std::string const& outputFileName_)
{
  init(outputFileName_);
}

Output::~Output()
{
//
}

hid_t Output::init(std::string const& outputFileName_)
{
  outputFile = H5Fcreate(outputFileName_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
  initDone = true;

  //Create base GroupIds and Description Attribute
  hid_t auxGroupID;//, auxDSpaceID, auxAttrID, auxType;
//  hsize_t auxDims[1] = {1};
//  string auxSData;


//  auxType = H5Tcopy(H5T_C_S1);

  auxGroupID = H5Gcreate(outputFile, "Field_E", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//  auxDSpaceID = H5Screate_simple(1, auxDims, NULL);
//  auxAttrID = H5Acreate(auxGroupID, "Description", auxType, auxDSpaceID, H5P_DEFAULT, H5P_DEFAULT);
//  H5Awrite(auxAttrID, H5T_NATIVE_INT, &auxData);
//  H5Aclose(auxAttrID);
//  H5Sclose(auxDSpaceID);
//  auxDSpaceID = H5Screate_simple(1, auxDims, NULL);
//  auxAttrID = H5Acreate(auxGroupID, "Status", auxType, auxDSpaceID, H5P_DEFAULT, H5P_DEFAULT);
//  auxSData = "Not Requested      ";
//  H5Awrite(auxAttrID, H5T_NATIVE_INT, auxSData.c_str());
//  H5Aclose(auxAttrID);
//  H5Sclose(auxDSpaceID);

  auxGroupID = H5Gcreate(outputFile, "Field_H", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//  auxDSpaceID = H5Screate_simple(1, auxDims, NULL);
//  auxAttrID = H5Acreate(auxGroupID, "Description", auxType, auxDSpaceID, H5P_DEFAULT, H5P_DEFAULT);
//  auxSData = "Magnetic Field     ";
//  H5Awrite(auxAttrID, H5T_NATIVE_INT, auxSData.c_str());
//  H5Aclose(auxAttrID);
//  H5Sclose(auxDSpaceID);
//  auxDSpaceID = H5Screate_simple(1, auxDims, NULL);
//  auxAttrID = H5Acreate(auxGroupID, "Status", auxType, auxDSpaceID, H5P_DEFAULT, H5P_DEFAULT);
//  auxSData = "Not Requested      ";
//  H5Awrite(auxAttrID, H5T_NATIVE_INT, auxSData.c_str());
//  H5Aclose(auxAttrID);
//  H5Sclose(auxDSpaceID);

  auxGroupID = H5Gcreate(outputFile, "CS_Sca", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//  auxDSpaceID = H5Screate_simple(1, auxDims, NULL);
//  auxAttrID = H5Acreate(auxGroupID, "Description", auxType, auxDSpaceID, H5P_DEFAULT, H5P_DEFAULT);
//  auxSData = "Scattering CS      ";
//  H5Awrite(auxAttrID, H5T_NATIVE_INT, auxSData.c_str());
//  H5Aclose(auxAttrID);
//  H5Sclose(auxDSpaceID);
//  auxDSpaceID = H5Screate_simple(1, auxDims, NULL);
//  auxAttrID = H5Acreate(auxGroupID, "Status", auxType, auxDSpaceID, H5P_DEFAULT, H5P_DEFAULT);
//  auxSData = "Not Requested      ";
//  H5Awrite(auxAttrID, H5T_NATIVE_INT, auxSData.c_str());
//  H5Aclose(auxAttrID);
//  H5Sclose(auxDSpaceID);

  auxGroupID = H5Gcreate(outputFile, "CS_Abs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//  auxDSpaceID = H5Screate_simple(1, auxDims, NULL);
//  auxAttrID = H5Acreate(auxGroupID, "Description", auxType, auxDSpaceID, H5P_DEFAULT, H5P_DEFAULT);
//  auxSData = "Absorption CS      ";
//  H5Awrite(auxAttrID, H5T_NATIVE_INT, auxSData.c_str());
//  H5Aclose(auxAttrID);
//  H5Sclose(auxDSpaceID);
//  auxDSpaceID = H5Screate_simple(1, auxDims, NULL);
//  auxAttrID = H5Acreate(auxGroupID, "Status", auxType, auxDSpaceID, H5P_DEFAULT, H5P_DEFAULT);
//  auxSData = "Not Requested      ";
//  H5Awrite(auxAttrID, H5T_NATIVE_INT, auxSData.c_str());
//  H5Aclose(auxAttrID);
//  H5Sclose(auxDSpaceID);

  auxGroupID = H5Gcreate(outputFile, "CS_Ext", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
//  auxDSpaceID = H5Screate_simple(1, auxDims, NULL);
//  auxAttrID = H5Acreate(auxGroupID, "Description", auxType, auxDSpaceID, H5P_DEFAULT, H5P_DEFAULT);
//  auxSData = "Extinction CS      ";
//  H5Awrite(auxAttrID, H5T_NATIVE_INT, auxSData.c_str());
//  H5Aclose(auxAttrID);
//  H5Sclose(auxDSpaceID);
//  auxDSpaceID = H5Screate_simple(1, auxDims, NULL);
//  auxAttrID = H5Acreate(auxGroupID, "Status", auxType, auxDSpaceID, H5P_DEFAULT, H5P_DEFAULT);
//  auxSData = "Not Requested      ";
//  H5Awrite(auxAttrID, H5T_NATIVE_INT, auxSData.c_str());
//  H5Aclose(auxAttrID);
//  H5Sclose(auxDSpaceID);

  return outputFile;
}

hid_t Output::getHandle(std::string code_)
{
  if(initDone)
    return H5Gopen(outputFile, code_.c_str(), H5P_DEFAULT);
  return -1;
}

void Output::close()
{
  if(initDone)
    H5Fclose(outputFile);
}


