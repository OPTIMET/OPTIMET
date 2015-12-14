include(AddRegression)
add_regression_test(OneParticleScan
  BLESSED
    ${CMAKE_SOURCE_DIR}/test-data/OneParticleScan_AbsorptionCS.dat
    ${CMAKE_SOURCE_DIR}/test-data/OneParticleScan_ExtinctionCS.dat
  OUTPUTS
    ${CMAKE_BINARY_DIR}/examples/OneParticleScan_AbsorptionCS.dat
    ${CMAKE_BINARY_DIR}/examples/OneParticleScan_ExtinctionCS.dat
  LABELS "slow"
)

if(HDF5_DIFF_EXECUTABLE)
  add_regression_test(OneParticle
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/OneParticle.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/OneParticle.h5
    DIFF_CMD ${HDF5_DIFF_EXECUTABLE} -d1.8e-12
    LABELS "slow"
  )
  add_regression_test(ThreeParticles
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/ThreeParticles.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/ThreeParticles.h5
    DIFF_CMD ${HDF5_DIFF_EXECUTABLE} -d6.4e-9
    LABELS "slow"
  )
  add_regression_test(SpiralStructure
    DISABLED
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/SpiralStructure.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/SpiralStructure.h5
    DIFF_CMD ${HDF5_DIFF_EXECUTABLE}
    LABELS "slow"
  )
else()
  message(WARNING "The hdf5 diff executable was not found: some regression tests will not br run.")
endif()
