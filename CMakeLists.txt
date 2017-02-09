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

cmake_minimum_required(VERSION 3.0)

# Location of extra cmake includes for the project
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
# Downloads and installs GreatCMakeCookOff
# It contains a number of cmake recipes
include(LookUp-GreatCMakeCookOff)

# Hunter helps us download stuff automatically
include(AddHunter)
set(HUNTER_ROOT ${CMAKE_BINARY_DIR}/hunter_root)
initialize_hunter(HUNTER_ROOT ${HUNTER_ROOT} LOCAL)

project(Optimet3D CXX C)

option(dompi "Enable mpi" off)
option(dotesting "Enable testing" on)
option(dobenchmarks "Enable Benchmarking" on)

# looks for all dependencies used by optimet
include(dependencies)

# adds flags for compiling with c++11
include(AddCPP11Flags)

include_directories(${PROJECT_SOURCE_DIR}/src ${PROJECT_SOURCE_DIR})
include_directories(SYSTEM
  ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS}
  ${F2C_INCLUDE_DIRS} ${Boost_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIR}
  ${Belos_INCLUDE_DIRS})

# list all object files
file(GLOB SRC src/*.c*  pugi/*.c*)
file(GLOB EXCLUDE .clang-format src/main.cpp)
list(REMOVE_ITEM SRC ${EXCLUDE})
if(dompi)
  set(OPTIMET_MPI TRUE)
  file(GLOB parallel_files src/scalapack/*.cpp src/mpi/*.cpp)
  list(APPEND SRC ${parallel_files})
  # Add fortran interface for blacs and scalapack
  include(DetectFortran)
  DetectFortran("${PROJECT_BINARY_DIR}/include/optimet/OptimetFC.h")
  include_directories(SYSTEM ${MPI_CXX_INCLUDE_PATH})
  # We need to check integer architecture
  include(DetectIntegerArchitecture)
  DetectIntegerArchitecture(OPTIMET)
else()
  file(GLOB EXCLUDE src/ScalapackSolver.cpp src/MatrixBelosSolver.cpp src/FMMBelosSolver.cpp)
  list(REMOVE_ITEM SRC ${EXCLUDE})
endif()

if(NOT OPTIMET_BELOS)
  file(GLOB EXCLUDE src/MatrixBelosSolver.cpp src/FMMBelosSolver.cpp)
  list(REMOVE_ITEM SRC ${EXCLUDE})
endif()

# configure a file with some build defaults
include_directories("${PROJECT_BINARY_DIR}/include/optimet")
configure_file(src/Types.in.h "${PROJECT_BINARY_DIR}/include/optimet/Types.h")

# add a library that can be linked against in unittests as well
add_library(optilib STATIC ${SRC})
if(TARGET Boost::boost)
  target_link_libraries(optilib Boost::boost)
endif()
set(library_dependencies
  ${GSL_LIBRARIES} ${BOOST_LIBRARIES} ${HDF5_C_LIBRARIES} ${F2C_LIBRARIES}
  ${Belos_LIBRARIES}
  )
if(dompi)
  list(APPEND library_dependencies ${MPI_LIBRARIES} ${SCALAPACK_LIBRARIES})
endif()
list(APPEND library_dependencies ${F2C_LIBRARIES})

# main executable
add_executable(Optimet3D src/main.cpp)
target_link_libraries(Optimet3D optilib ${library_dependencies})

# all testing stuff below
if(dotesting)
  enable_testing()

  include(AddCatchTest)
  add_subdirectory(unittests)
  if(EXISTS "${PROJECT_SOURCE_DIR}/test-data/OneParticle.h5")
    include(regression)
  else()
    message(WARNING "The git submodule test-data has not been downloaded.\n"
      "Do `git submodule update --init --recursive` to do so now.\n"
      "Until then, tests that need this data will be disabled")
  endif()
endif()

if(dobenchmarks)
  include(AddBenchmark)
  add_subdirectory(benchmarks)
endif()
