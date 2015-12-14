include(CMakeParseArguments)
function(add_regression_test testname)
  cmake_parse_arguments(regr "DISABLE" "INPUTFILE" "DIFF_CMD;BLESSED;OUTPUTS;LABELS" ${ARGN})
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
  set(regr_CMD ${PROJECT_BINARY_DIR}/Optimet3D ${regr_INPUTFILE})
  # Creates regression test file
  file(MAKE_DIRECTORY "${CMAKE_BINARY_DIR}/regressions/")
  configure_file(
    "${PROJECT_SOURCE_DIR}/cmake/regression.in.cmake"
    "${CMAKE_BINARY_DIR}/regressions/${testname}.cmake" @ONLY
  )
  if(NOT regr_DISABLE)
    add_test(
      NAME ${testname}
      COMMAND ${CMAKE_COMMAND} -P ${PROJECT_SOURCE_DIR}/regression.in.cmake
      WORKING_DIRECTORY ${PROJECT_BINARY_DIR}
    )
    list(APPEND regr_LABELS "regression")
    set_tests_properties(${testname} PROPERTIES LABELS "${regr_LABELS}")
  endif()
endfunction()
