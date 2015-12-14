#! @CMAKE_COMMAND@ -P
# Sometimes HDF5 fails to manipulate a previously written file
# Might be because OPTIMET, might be because HDF5
set(output_blessed @regr_BLESSED@)
set(test_outputs @regr_OUTPUTS@)
set(diff_cmd @regr_DIFF_CMD@)

if("${diff_cmd}" STREQUAL "")
  set(diff_cmd ${CMAKE_COMMAND} -E compare_files)
endif()

foreach(output ${test_blessed})
  if(EXISTS "${output}")
    file(REMOVE "${output}")
  endif()
endforeach()

execute_process(COMMAND @regr_CMD@ RESULT_VARIABLE test_not_successful)

if( test_not_successful )
  message( SEND_ERROR "@regr_CMD@ returned ${test_not_successful}" )
endif( test_not_successful )

list(LENGTH output_blessed len1)
math(EXPR len2 "${len1} - 1")

# diff_cmd contains the executable that will be invoked to compare the output
# files
foreach(index RANGE ${len2})
  list(GET output_blessed ${index} blessed)
  list(GET test_outputs ${index} test)
  execute_process(
    COMMAND ${diff_cmd};${blessed};${test}
    OUTPUT_VARIABLE test_output
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if(NOT test_output STREQUAL "")
    message(SEND_ERROR "${test} does not match ${blessed}! ${test_output} !")
  endif()
endforeach()
