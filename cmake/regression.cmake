# (C) University College London 2017
# This file is part of Optimet, licensed under the terms of the GNU Public License
#
# Optimet is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Optimet is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Optimet. If not, see <http:#www.gnu.org/licenses/>.

include(AddRegression)
if(dompi)
  set(DOMPI DOMPI)
endif()
add_regression_test(OneParticleScan
  DISABLE # comparison not done right
  ${DOMPI}
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
    ${DOMPI}
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/OneParticle.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/OneParticle.h5
    HDF5_PRECISION 2.3e-12
    LABELS "slow"
  )
  add_regression_test(ThreeParticles
    ${DOMPI}
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/ThreeParticles.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/ThreeParticles.h5
    HDF5_PRECISION 6.4e-9
    LABELS "slow"
  )
  add_regression_test(SpiralStructure
    ${DOMPI}
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/SpiralStructure.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/SpiralStructure.h5
    HDF5_PRECISION 1.0e-8
    LABELS "medium"
  )

  add_regression_test(FastOneParticle
    ${DOMPI}
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/FastOneParticle.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/FastOneParticle.h5
    HDF5_PRECISION 1e-11
    LABELS "fast"
  )
  add_regression_test(FastThreeParticles
    ${DOMPI}
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/FastThreeParticles.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/FastThreeParticles.h5
    HDF5_PRECISION 1e-11
    LABELS "fast"
  )
  add_regression_test(FastSpiralStructure
    ${DOMPI}
    BLESSED ${CMAKE_SOURCE_DIR}/test-data/FastSpiralStructure.h5
    OUTPUTS ${CMAKE_BINARY_DIR}/examples/FastSpiralStructure.h5
    HDF5_PRECISION 1e-11
    LABELS "fast"
  )
else()
  message(WARNING "The hdf5 diff executable was not found: some regression tests will not be run.")
endif()
