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


include(CMakeParseArguments)
macro(add_hunter_package package)
  cmake_parse_arguments(ah "CHECK_SYSTEM" "PACKAGE" "" ${ARGN})
  if(NOT ah_PACKAGE)
    set(ah_PACKAGE ${package})
  endif()
  string(TOUPPER "${ah_PACKAGE}" AH_PACKAGE)
  set(do_hunter true)
  if(ah_CHECK_SYSTEM AND NOT ${ah_PACKAGE}_ADDED_VIA_HUNTER)
    find_package(${ah_PACKAGE} ${ah_UNPARSED_ARGUMENTS})
    if(${ah_PACKAGE}_FOUND OR ${AH_PACKAGE}_FOUND)
      set(do_hunter false)
    endif()
  endif()
  if(do_hunter)
    hunter_add_package(${package} ${ah_UNPARSED_ARGUMENTS})
    find_package(${ah_PACKAGE} REQUIRED ${ah_UNPARSED_ARGUMENTS})
    set(${ah_PACKAGE}_ADDED_VIA_HUNTER TRUE CACHE "${ah_PACKAGE} added via hunter" FORCE)
    mark_as_advanced(${ah_PACKAGE}_ADDED_VIA_HUNTER)
  endif()
endmacro()

add_hunter_package(Boost CHECK_SYSTEM)
add_hunter_package(Eigen PACKAGE Eigen3 CHECK_SYSTEM)
add_hunter_package(hdf5 PACKAGE HDF5 CHECK_SYSTEM)
add_hunter_package(GSL CHECK_SYSTEM)
include(PackageLookup)
