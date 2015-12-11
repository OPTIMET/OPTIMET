cmake_minimum_required(VERSION 3.0)

project(Optimet3D CXX C)

option(dompi "Enable mpi" on)
option(dotesting "Enable testing" on)

# Location of extra cmake includes for the project
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# Downloads and installs GreatCMakeCookOff
# It contains a number of cmake recipes
include(LookUp-GreatCMakeCookOff)

include(dependencies)


include(AddCPP11Flags)

include_directories(
  ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ${F2C_INCLUDE_DIRS})
if(EIGEN3_INCLUDE_DIR)
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()

file(GLOB SRC *.c* pugi/*.c*)
add_executable(Optimet3D ${SRC})
target_link_libraries(Optimet3D
  ${GSL_LIBRARIES} ${HDF5_LIBRARIES} ${ARMADILLO_LIBRARIES} ${F2C_LIBRARIES})


if(dotesting)
  enable_testing()
  if(EXISTS "${PROJECT_SOURCE_DIR}/test-data/OneParticle.h5")
    include(regression)
  else()
    message(WARNING "The git submodule test-data has not been downloaded.\n"
      "Do `git submodule update --init --recursive` to do so now.\n"
      "Until then, tests that need this data will be disabled")
  endif()
endif()
