if(OPTIMET_BELOS)
  add_benchmark(benchmark benchmarks.cpp fcc.cpp cmdl.cpp LIBRARIES optilib ${library_dependencies})

endif()
if(OPTIMET_MPI)
  add_executable(parallel_fmm_multiplication fmm_multiplication.cpp)
  target_link_libraries(parallel_fmm_multiplication optilib ${library_dependencies})
endif()

add_executable(serial_fmm_multiplication fmm_multiplication.cpp)
target_link_libraries(serial_fmm_multiplication optilib ${library_dependencies})
target_compile_definitions(serial_fmm_multiplication PRIVATE OPTIMET_JUST_DO_SERIAL)
