include(AddRegression)
add_regression_test(OneParticleScan
  DISABLE # comparison not done right
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
    HDF5_PRECISION 1.8e-12
    LABELS "slow"
  )
  add_regression_test(ThreeParticles
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/ThreeParticles.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/ThreeParticles.h5
    HDF5_PRECISION 6.4e-9
    LABELS "slow"
  )
  add_regression_test(SpiralStructure
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/SpiralStructure.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/SpiralStructure.h5
    HDF5_PRECISION 1.0e-8
    LABELS "medium"
  )

  add_regression_test(FastOneParticle
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/FastOneParticle.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/FastOneParticle.h5
    HDF5_PRECISION 1e-11
    LABELS "fast"
  )
  add_regression_test(FastThreeParticles
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/FastThreeParticles.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/FastThreeParticles.h5
    HDF5_PRECISION 1e-11
    LABELS "fast"
  )
  add_regression_test(FastSpiralStructure
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/FastSpiralStructure.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/FastSpiralStructure.h5
    HDF5_PRECISION 1e-11
    LABELS "fast"
  )
else()
  message(WARNING "The hdf5 diff executable was not found: some regression tests will not br run.")
endif()
