cmake_minimum_required(VERSION 3.2)
project(Optimet3DP)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

find_package(MKL REQUIRED)
find_package(HDF5 REQUIRED)
find_package(GSL REQUIRED)
find_package(MPI REQUIRED)
find_package(F2C REQUIRED)

set(MKL_LIBRARIES iomp5 mkl_scalapack_lp64 mkl_intel_lp64 mkl_core
    mkl_intel_thread mkl_blacs_intelmpi_lp64 pthread m)

include_directories(${MKL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${GSL_INCLUDE_DIRS}
                    ${MPI_INCLUDE_DIRS} ${F2C_INCLUDE_DIRS})

file(GLOB SRC ${CMAKE_SOURCE_DIR}/*.c*)

add_executable(Optimet3DP ${SRC})
target_link_libraries(Optimet3DP ${MKL_LIBRARIES} ${HDF5_LIBRARIES}
                      ${GSL_LIBRARIES} ${MPI_LIBRARIES} ${F2C_LIBRARIES})
