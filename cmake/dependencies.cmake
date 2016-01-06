find_package(F2C REQUIRED)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)


macro(add_hunter_package package)
  find_package(${package} ${ARGN})
  if(NOT ${package}_FOUND)
    hunter_add_package(${package} ${ARGN})
    find_package(${package} REQUIRED ${ARGN})
  endif()
endmacro()

add_hunter_package(Boost)
find_package(Eigen3)
if(NOT Eigen3_FOUND)
  hunter_add_package(Eigen)
  find_package(Eigen3 REQUIRED)
endif()
add_hunter_package(GSL)
add_hunter_package(HDF5)
include(PackageLookup)
