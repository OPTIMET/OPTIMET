cmake_minimum_required(VERSION 3.0)

project(Optimet3D CXX C)

option(dompi "Enable mpi" on)
option(dotesting "Enable testing" on)

# Location of extra cmake includes for the project
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
# Downloads and installs GreatCMakeCookOff
# It contains a number of cmake recipes
include(LookUp-GreatCMakeCookOff)

# looks for all dependencies used by optimet
include(dependencies)

# adds flags for compiling with c++11
include(AddCPP11Flags)

include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${F2C_INCLUDE_DIRS})
if(EIGEN3_INCLUDE_DIR)
  include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
endif()

# list all object files
file(GLOB SRC *.c* pugi/*.c*)
list(REMOVE_ITEM SRC main.cpp)

# add a library that can be linked against in unittests as well
add_library(optilib STATIC ${SRC})
add_dependencies(optilib lookup_dependencies)
set(library_dependencies ${GSL_LIBRARIES} ${HDF5_LIBRARIES} ${F2C_LIBRARIES})

# main executable
add_executable(Optimet3D main.cpp)
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
