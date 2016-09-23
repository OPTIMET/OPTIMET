include_directories("${PROJECT_SOURCE_DIR}")
add_library(common_catch_main_object OBJECT common_catch_main.cpp)
if(CATCH_INCLUDE_DIR)
  target_include_directories(common_catch_main_object SYSTEM PUBLIC ${CATCH_INCLUDE_DIR})
endif()

if(NOT dompi)
  add_catch_test(single_particle LIBRARIES optilib ${library_dependencies})
endif()
add_catch_test(translation_addition LIBRARIES optilib ${library_dependencies})
add_catch_test(harmonics_iterator LIBRARIES optilib ${library_dependencies})
add_catch_test(aux_coefficients LIBRARIES optilib ${library_dependencies})
add_catch_test(rotation_coaxial_decomposition LIBRARIES optilib ${library_dependencies})

add_catch_test(rotation_coefficients LIBRARIES optilib ${library_dependencies})

if(dompi)
  if(MPIEXEC_MAX_NUMPROCS LESS 2)
    message(STATUS "MPI tests are disabled: MPIEXEC_MAX_NUMPROCS=${MPIEXEC_MAX_NUMPROCS}")
    return()
  elseif(MPIEXEC_MAX_NUMPROCS LESS 4)
    message(STATUS "Some MPI tests are disabled: MPIEXEC_MAX_NUMPROCS < 4")
    return()
  endif()

  add_library(common_mpi_catch_main_object OBJECT common_mpi_catch_main.cpp)
  if(CATCH_INCLUDE_DIR)
    target_include_directories(common_mpi_catch_main_object SYSTEM PUBLIC ${CATCH_INCLUDE_DIR})
  endif()
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
  add_mpi_test(single_particle LIBRARIES optilib ${library_dependencies})
  add_mpi_test(mpi_communicator LIBRARIES optilib ${library_dependencies})
  add_mpi_test(scalapack_context LIBRARIES optilib ${library_dependencies})
  add_mpi_test(scalapack_matrix LIBRARIES optilib ${library_dependencies})
  add_mpi_test(scalapack_linear_system LIBRARIES optilib ${library_dependencies})
  add_mpi_test(scalapack_collectives LIBRARIES optilib ${library_dependencies})
  add_mpi_test(serial_vs_parallel LIBRARIES optilib ${library_dependencies})
  add_mpi_test(mpi_scattering_matrix LIBRARIES optilib ${library_dependencies})
  if(OPTIMET_BELOS)
    add_mpi_test(belos LIBRARIES optilib ${library_dependencies})
    add_mpi_test(scalapack_vs_belos LIBRARIES optilib ${library_dependencies})

    # Prints all available solvers and their parameters
    add_executable(print_belos_params print_belos_params.cpp)
    target_link_libraries(print_belos_params optilib ${library_dependencies})
  endif()
endif()
