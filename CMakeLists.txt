cmake_minimum_required(VERSION 3.0)

# Location of extra cmake includes for the project
list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
include(HunterGate)
if(NOT HUNTER_ROOT AND "$ENV{HUNTER_ROOT}" STREQUAL "")
  set(HUNTER_ROOT "${CMAKE_CURRENT_BINARY_DIR}/hunter" CACHE PATH "Hunter root dir")
endif()
HunterGate(
  URL "https://github.com/ruslo/hunter/archive/v0.12.33.tar.gz"
  SHA1 "f64eda11326e12d84b6aa7f4aed80c09c182aa64"
)

project(Optimet3D CXX C)

option(dompi "Enable mpi" on)
option(dotesting "Enable testing" on)

# Downloads and installs GreatCMakeCookOff
# It contains a number of cmake recipes
include(LookUp-GreatCMakeCookOff)

# looks for all dependencies used by optimet
include(dependencies)

# adds flags for compiling with c++11
include(AddCPP11Flags)

include_directories(SYSTEM ${GSL_INCLUDE_DIRS} ${BOOST_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${F2C_INCLUDE_DIRS})
if(EIGEN3_INCLUDE_DIR)
  include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})
endif()

# list all object files
file(GLOB SRC *.c* pugi/*.c*)
list(REMOVE_ITEM SRC main.cpp)

# add a library that can be linked against in unittests as well
add_library(optilib STATIC ${SRC})
add_dependencies(optilib lookup_dependencies)
set(library_dependencies ${GSL_LIBRARIES} ${BOOST_LIBRARIES} ${HDF5_LIBRARIES} ${F2C_LIBRARIES})

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
