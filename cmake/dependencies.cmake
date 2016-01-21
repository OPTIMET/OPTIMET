find_package(BOOST 1.35.0 REQUIRED COMPONENTS math)
find_package(GSL REQUIRED)
find_package(HDF5 REQUIRED)
find_package(F2C REQUIRED)

include(PackageLookup)
lookup_package(Eigen3 REQUIRED)

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)
