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

#ifndef OUTPUT_H_
#define OUTPUT_H_

#include <string>
#include <hdf5.h>

/**
 * The Output class implements HDF5 output.
 * This class does not do any data processing. It simply initializes,
 * prepares and finally closes the HDF5 file.
 *
 * Initialization assumes that new "base" Groups are created for all possible
 * output data and a corresponding description attribute. The description will
 * be originally set to "Data not requested.".
 *
 * The current base Groups are:
 *    Field_E - Total Electric Field
 *    Field_H - Total Magnetic Field
 *    CS_Sca  - Scattering Cross Section
 *    CS_Abs  - Absorption Cross Section
 *    CS_Ext  - Extinction Cross Section
 *
 * The groupIds are then passed to the various OutputX objects
 * (e.g. OutputGrid, OutputScan, etc.) by the Result object.
 */
class Output {
private:
  bool initDone;

public:
  hid_t outputFile; /**< Handle to the HDF5 output file. */

  /**
   * Default constructor for the Output class.
   * Does NOT initialize the Output object.
   */
  Output();

  /**
   * Initialization constructor for the Output class.
   * Sets the iterator to the starting position.
   * @param outputFileName_ the name of the hdf5 output file.
   */
  Output(std::string const &outputFileName_);

  /**
   * Default destructor for the Output class.
   * Closes the HDF5 file if it was opened.
   */
  virtual ~Output();

  /**
   * Initialization method for the Output class.
   * @param outputFileName_ the name of the hdf5 output file.
   * @return the handle to the HDF5 file.
   */
  hid_t init(std::string const &outputFileName_);

  /**
   * Returns the handle to the base GroupID.
   * @param code_ the GroupID code (see class documentation).
   * @return the handle to the the base GroupID.
   */
  hid_t getHandle(std::string code_);

  void close();
};

#endif /* OUTPUT_H_ */
