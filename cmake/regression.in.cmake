#! @CMAKE_COMMAND@ -P
# Sometimes HDF5 fails to manipulate a previously written file
# Might be because OPTIMET, might be because HDF5
set(output_blessed @regr_BLESSED@)
set(test_outputs @regr_OUTPUTS@)
set(diff_cmd @regr_DIFF_CMD@)
unset(hdf5_args)
if(NOT "@regr_HDF5_PRECISION@" STREQUAL "")
  set(hdf5_args -d@regr_HDF5_PRECISION@)
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
  get_filename_component(extension "${output_blessed}" EXT)
  if(diff_cmd)
    set(cmd ${diff_cmd} ${blessed} ${test})
  elseif(extension STREQUAL ".h5")
    set(cmd @HDF5_DIFF_EXECUTABLE@ -v ${hdf5_args} ${blessed} ${test})
  else()
    set(cmd ${CMAKE_COMMAND};-E;compare_files;${blessed};${test})
  endif()
  execute_process(
    COMMAND ${cmd}
    ERROR_VARIABLE test_error
    OUTPUT_VARIABLE test_output
    RESULT_VARIABLE test_result
    OUTPUT_STRIP_TRAILING_WHITESPACE
  )

  if(NOT test_result EQUAL 0)
    message(STATUS "${test_output}")
    message(STATUS "${test_error}")
    message(SEND_ERROR "${test} does not match ${blessed}!")
  endif()
endforeach()
