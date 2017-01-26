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
