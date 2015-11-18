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

find_package(Doxygen)
if(DOXYGEN_FOUND)
  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/doxygen.config
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config @ONLY)
  add_custom_target(doc ${DOXYGEN_EXECUTABLE}
    ${CMAKE_CURRENT_BINARY_DIR}/doxygen.config
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    COMMENT "Generating API documentation with Doxygen" VERBATIM)
endif(DOXYGEN_FOUND)

enable_testing()

configure_file(examples/OneParticle.xml examples/OneParticle.xml COPYONLY)
add_test(NAME OneParticle
  COMMAND ${CMAKE_COMMAND}
  -Dtest_cmd=${CMAKE_BINARY_DIR}/Optimet3D\;${CMAKE_BINARY_DIR}/examples/OneParticle
  -Ddiff_cmd=${HDF5_DIFF_EXECUTABLE}
  -Doutput_blessed=${CMAKE_SOURCE_DIR}/test-data/OneParticle.h5
  -Doutput_test=${CMAKE_BINARY_DIR}/examples/OneParticle.h5
  -P ${CMAKE_SOURCE_DIR}/cmake/run_test.cmake
)

configure_file(examples/OneParticleScan.xml examples/OneParticleScan.xml COPYONLY)
add_test(NAME OneParticleScan
  COMMAND ${CMAKE_COMMAND}
  -Dtest_cmd=${CMAKE_BINARY_DIR}/Optimet3D\;${CMAKE_BINARY_DIR}/examples/OneParticleScan
  -Doutput_blessed=${CMAKE_SOURCE_DIR}/test-data/OneParticleScan_AbsorptionCS.dat\;${CMAKE_SOURCE_DIR}/test-data/OneParticleScan_ExtinctionCS.dat
  -Doutput_test=${CMAKE_BINARY_DIR}/examples/OneParticleScan_AbsorptionCS.dat\;${CMAKE_BINARY_DIR}/examples/OneParticleScan_ExtinctionCS.dat
  -P ${CMAKE_SOURCE_DIR}/cmake/run_test.cmake
)

configure_file(examples/ThreeParticles.xml examples/ThreeParticles.xml COPYONLY)
add_test(NAME ThreeParticles
  COMMAND ${CMAKE_COMMAND}
  -Dtest_cmd=${CMAKE_BINARY_DIR}/Optimet3D\;${CMAKE_BINARY_DIR}/examples/ThreeParticles
  -Ddiff_cmd=${HDF5_DIFF_EXECUTABLE}
  -Doutput_blessed=${CMAKE_SOURCE_DIR}/test-data/ThreeParticles.h5
  -Doutput_test=${CMAKE_BINARY_DIR}/examples/ThreeParticles.h5
  -P ${CMAKE_SOURCE_DIR}/cmake/run_test.cmake
)

configure_file(examples/SpiralStructure.xml examples/SpiralStructure.xml COPYONLY)
add_test(NAME SpiralStructure
  COMMAND ${CMAKE_COMMAND}
  -Dtest_cmd=${CMAKE_BINARY_DIR}/Optimet3D\;${CMAKE_BINARY_DIR}/examples/SpiralStructure
  -Ddiff_cmd=${HDF5_DIFF_EXECUTABLE}
  -Doutput_blessed=${CMAKE_SOURCE_DIR}/test-data/SpiralStructure.h5
  -Doutput_test=${CMAKE_BINARY_DIR}/examples/SpiralStructure.h5
  -P ${CMAKE_SOURCE_DIR}/cmake/run_test.cmake
)
