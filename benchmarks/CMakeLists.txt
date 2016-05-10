add_library(benchmark_helpers STATIC cmdl.cpp)
add_benchmark(fcc LIBRARIES ${library_dependencies} benchmark_helpers optilib)
