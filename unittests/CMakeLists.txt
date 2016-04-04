include_directories("${PROJECT_SOURCE_DIR}")
add_catch_test(single_particle LIBRARIES optilib ${library_dependencies})
add_catch_test(translation_addition LIBRARIES optilib ${library_dependencies})
add_catch_test(harmonics_iterator LIBRARIES optilib ${library_dependencies})

if(dompi)
  if(CATCH_INCLUDE_DIR)
    include_directories(${CATCH_INCLUDE_DIR})
  endif()
  if(MPIEXEC_MAX_NUMPROCS LESS 2)
    message(STATUS "MPI tests are disabled: MPIEXEC_MAX_NUMPROCS=${MPIEXEC_MAX_NUMPROCS}")
    return()
  elseif(MPIEXEC_MAX_NUMPROCS LESS 4)
    message(STATUS "Some MPI tests are disabled: MPIEXEC_MAX_NUMPROCS < 4")
  endif()

  add_library(common_mpi_catch_main_object OBJECT common_mpi_catch_main.cpp)
  function(add_mpi_test testname)
    add_catch_test(${testname}
      COMMON_MAIN common_mpi_catch_main_object
      NOTEST ${ARGN})
    unset(arguments)
    if(CATCH_JUNIT)
      set(arguments -r junit -o ${PROJECT_BINARY_DIR}/Testing/${testname}.xml)
    endif()
    add_test(
      NAME ${testname}
      COMMAND
      ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS}
        test_${testname} ${arguments}
    )
    set_tests_properties(${testname} PROPERTIES LABELS "catch;mpi")
  endfunction()

  add_catch_test(squarest LIBRARIES optilib ${library_dependencies})
  add_mpi_test(mpi_communicator LIBRARIES optilib ${library_dependencies})
  add_mpi_test(scalapack_context LIBRARIES optilib ${library_dependencies})
  add_mpi_test(scalapack_matrix LIBRARIES optilib ${library_dependencies})
  add_mpi_test(scalapack_linear_system LIBRARIES optilib ${library_dependencies})
  add_mpi_test(scalapack_collectives LIBRARIES optilib ${library_dependencies})
  add_mpi_test(serial_vs_parallel LIBRARIES optilib ${library_dependencies})
endif()
