# some argument checking:
# test_cmd is the command to run with all its arguments
if( NOT test_cmd )
  message( FATAL_ERROR "Variable test_cmd not defined" )
endif( NOT test_cmd )

# output_blessed contains a list of the "blessed" output files
if( NOT output_blessed )
  message( FATAL_ERROR "Variable output_blessed not defined" )
endif( NOT output_blessed )

# output_test contains a list of the output files the test_cmd will produce
# (in the same order as output_blessed above)
if( NOT output_test )
  message( FATAL_ERROR "Variable output_test not defined" )
endif( NOT output_test )

# Sometimes HDF5 fails to manipulate a previously written file
# Might be because OPTIMET, might be because HDF5
foreach(output ${output_test})
  if(EXISTS "${output}")
    file(REMOVE "${output}")
  endif()
endforeach()

execute_process(
  COMMAND ${test_cmd}
  RESULT_VARIABLE test_not_successful
)

if( test_not_successful )
  message( SEND_ERROR "${test_cmd} returned ${test_not_successful}" )
endif( test_not_successful )

list(LENGTH output_blessed len1)
math(EXPR len2 "${len1} - 1")

# diff_cmd contains the executable that will be invoked to compare the output
# files
if( diff_cmd )
  foreach(index RANGE ${len2})
    list(GET output_blessed ${index} blessed)
    list(GET output_test ${index} test)
    execute_process(
      COMMAND ${diff_cmd} ${blessed} ${test}
      RESULT_VARIABLE test_not_successful
      OUTPUT_QUIET
      ERROR_QUIET
    )
  endforeach()
else( diff_cmd )
  foreach(index RANGE ${len2})
    list(GET output_blessed ${index} blessed)
    list(GET output_test ${index} test)
    execute_process(
      COMMAND ${CMAKE_COMMAND} -E compare_files ${blessed} ${test}
      RESULT_VARIABLE test_not_successful
      OUTPUT_QUIET
      ERROR_QUIET
    )
  endforeach()
endif( diff_cmd )

if( test_not_successful )
  message( SEND_ERROR "${output_test} does not match ${output_blessed}!")
endif( test_not_successful )
