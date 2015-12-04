configure_file(examples/OneParticle.xml examples/OneParticle.xml COPYONLY)
add_test(NAME OneParticle
  COMMAND ${CMAKE_COMMAND}
  -Dtest_cmd=${CMAKE_BINARY_DIR}/Optimet3D\;${CMAKE_BINARY_DIR}/examples/OneParticle
  -Ddiff_cmd=${HDF5_DIFF_EXECUTABLE}\;-d1.8e-12
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
  -Ddiff_cmd=${HDF5_DIFF_EXECUTABLE}\;-d6.4e-09
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

configure_file(examples/2TouchingSpheres.xml examples/2TouchingSpheres.xml COPYONLY)
