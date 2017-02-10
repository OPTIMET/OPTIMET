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

include(CMakeParseArguments)
function(add_regression_test testname)
  cmake_parse_arguments(regr
    "DISABLE;DOMPI" "INPUTFILE;HDF5_PRECISION" "DIFF_CMD;BLESSED;OUTPUTS;LABELS" ${ARGN})
  if(NOT regr_BLESSED)
    message(FATAL_ERROR "No blessed outputs given for test ${testname}")
  endif()
  if(NOT regr_OUTPUTS)
    message(FATAL_ERROR "No actual outputs given for test ${testname}")
  endif()
  list(LENGTH regr_BLESSED l0)
  list(LENGTH regr_OUTPUTS l1)
  if(NOT l0 EQUAL l1)
    message(FATAL_ERROR "Number of blessed and outputs do not match: ${l0} vs ${l1}")
  endif()
  if(NOT regr_INPUTFILE)
    set(regr_INPUTFILE "examples/${testname}")
  endif()
  # Copies input file to examples directory
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/${regr_INPUTFILE}.xml")
    message(FATAL_ERROR "Input file ${PROJECT_SOURCE_DIR}/${regr_INPUTFILE}"
      " could not be found for test ${testname}")
  endif()
  configure_file(
    "${PROJECT_SOURCE_DIR}/${regr_INPUTFILE}.xml"
    "${PROJECT_BINARY_DIR}/${regr_INPUTFILE}.xml"
    COPYONLY
  )
  # Command is always optimet + input
  set(regr_CMD ${PROJECT_BINARY_DIR}/Optimet3D ${regr_INPUTFILE}.xml)
  if(regr_DOMPI)
    set(regr_CMD
      ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS}
      ${regr_CMD})
  endif()
  # Creates regression test file
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/regressions/")
  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/regression.in.cmake"
    "${CMAKE_BINARY_DIR}/regressions/${testname}.cmake" @ONLY
  )
  if(NOT regr_DISABLE)
    add_test(
      NAME ${testname}
      COMMAND ${CMAKE_COMMAND} -P ${CMAKE_BINARY_DIR}/regressions/${testname}.cmake
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
    list(APPEND regr_LABELS "regression")
    set_tests_properties(${testname} PROPERTIES LABELS "${regr_LABELS}")
  endif()
endfunction()
