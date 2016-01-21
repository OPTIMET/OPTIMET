include_directories("${PROJECT_SOURCE_DIR}")
add_catch_test(single_particle LIBRARIES optilib ${library_dependencies})
add_catch_test(translation_addition LIBRARIES optilib ${library_dependencies})

if(dompi)
  add_library(common_mpi_catch_main_object OBJECT common_mpi_catch_main.cpp)
  function(add_mpi_test testname)
    add_catch_test(${testname}
      COMMON_MAIN common_mpi_catch_main_object
      PRECOMMAND "${MPIEXEC} ${MPIEXEC_PREFLAGS} 4 ${MPIEXEC_PREFLAGS}" ${ARGN})
  endfunction()

  add_mpi_test(mpi_communicator LIBRARIES optilib ${library_dependencies})
endif()
