#include "Output.h"

Output::Output() { initDone = false; }

Output::Output(std::string const &outputFileName_) { init(outputFileName_); }

Output::~Output() {}

hid_t Output::init(std::string const &outputFileName_) {
  outputFile = H5Fcreate(outputFileName_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                         H5P_DEFAULT);
  initDone = true;

  // Create base GroupIds and Description Attribute
  hid_t auxGroupID; //, auxDSpaceID, auxAttrID, auxType;

  auxGroupID =
      H5Gcreate(outputFile, "Field_E", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(auxGroupID);

  auxGroupID =
      H5Gcreate(outputFile, "Field_H", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(auxGroupID);

  auxGroupID =
      H5Gcreate(outputFile, "CS_Sca", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(auxGroupID);

  auxGroupID =
      H5Gcreate(outputFile, "CS_Abs", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(auxGroupID);

  auxGroupID =
      H5Gcreate(outputFile, "CS_Ext", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Gclose(auxGroupID);

  return outputFile;
}

hid_t Output::getHandle(std::string code_) {
  if (initDone)
    return H5Gopen(outputFile, code_.c_str(), H5P_DEFAULT);
  return -1;
}

void Output::close() {
  if (initDone)
    H5Fclose(outputFile);
}
