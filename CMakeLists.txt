cmake_minimum_required(VERSION 3.2)

project(Optimet3D)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

find_package(GSL REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(F2C REQUIRED)

include_directories(${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ${F2C_INCLUDE_DIRS})

file(GLOB SRC ${CMAKE_SOURCE_DIR}/*.c* ${CMAKE_SOURCE_DIR}/pugi/*.c*)

add_executable(Optimet3D ${SRC})
target_link_libraries(Optimet3D ${GSL_LIBRARIES} ${HDF5_LIBRARIES} ${ARMADILLO_LIBRARIES} ${F2C_LIBRARIES})

set_target_properties(Optimet3D PROPERTIES CXX_STANDARD 11)
if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # CXX_STANDARD does not work with intel compiler yet
  list(APPEND CMAKE_CXX_FLAGS -std=c++11)
endif()

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

if(EXISTS "${PROJECT_SOURCE_DIR}/test-data/OneParticle.h5")
  enable_testing()
  include(regression)
else()
  message(WARNING "The git submodule test-data has not been downloaded.\n"
    "Do `git submodule update --init --recursive` to do so now.\n"
    "Until then, tests that need this data will be disabled")
endif()
