# (C) University College London 2017
# This file is part of Optimet, licensed under the terms of the GNU Public License
#
# Optimet is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Optimet is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Optimet. If not, see <http:#www.gnu.org/licenses/>.

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)


include(CMakeParseArguments)

if(NOT CMAKE_PREFIX_PATH AND NOT "$ENV{CMAKE_PREFIX_PATH}" STREQUAL "")
  string(REGEX REPLACE ":" ";" CMAKE_PREFIX_PATH $ENV{CMAKE_PREFIX_PATH})
endif()
find_or_add_hunter_package(Boost)
find_or_add_hunter_package(Eigen PACKAGE Eigen3)
find_or_add_hunter_package(hdf5 PACKAGE HDF5 COMPONENTS C)
find_or_add_hunter_package(GSL)
find_or_add_hunter_package(F2C)

if(dobenchmarks)
  find_or_add_hunter_package(GBenchmark)
endif()

set(OPTIMET_SCALAPACK FALSE)
if(dompi AND "$ENV{CRAYOS_VERSION}" STREQUAL "")
  find_package(MPI REQUIRED)

  if(NOT DEFINED BLA_VENDOR OR BLA_VENDOR STREQUAL "All" OR BLA_VENDOR MATCHES "Intel")
    find_package(MKL COMPONENTS ScaLAPACK)
    set(scalapack_FOUND ${MKL_ScaLAPACK_FOUND})
    set(SCALAPACK_LIBRARIES -Wl,--start-group ${MKL_LIBRARIES} -Wl,--end-group ${MPI_CXX_LIBRARIES})
  endif()

  if(NOT MKL_FOUND)
    find_package(scalapack)
    if(scalapack_FOUND)
      set(SCALAPACK_LIBRARIES scalapack)
    endif()
  endif()

  if (scalapack_FOUND)
    message(STATUS "Using ScaLAPACK libraries: ${SCALAPACK_LIBRARIES}")
    set(OPTIMET_SCALAPACK TRUE)
  endif()
elseif(NOT "$ENV{CRAYOS_VERSION}" STREQUAL "")
  unset(SCALAPACK_LIBRARIES)
  set(OPTIMET_SCALAPACK TRUE)
endif()
if(dompi AND NOT MPIEXEX_MAX_NUMPROCS)
  set(MPIEXEC_MAX_NUMPROCS 6)
endif()

# GMRes and other solvers
find_package(Belos)
set(OPTIMET_BELOS ${Belos_FOUND})
if(dompi AND NOT Belos_FOUND)
  message(STATUS "Compiling without Belos solvers")
endif()
