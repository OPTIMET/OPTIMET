if(OPTIMET_BELOS)
  add_benchmark(benchmark benchmarks.cpp fcc.cpp cmdl.cpp LIBRARIES optilib ${library_dependencies})
endif()
