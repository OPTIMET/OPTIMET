cmake_minimum_required(VERSION 2.6)

project(Optimet3D)

find_package(GSL REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Armadillo REQUIRED)

include_directories(${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS})

file(GLOB SRC ${CMAKE_SOURCE_DIR}/*.c* ${CMAKE_SOURCE_DIR}/pugi/*.c*)

add_executable(Optimet3D ${SRC})
target_link_libraries(Optimet3D ${GSL_LIBRARIES} ${HDF5_LIBRARIES} ${ARMADILLO_LIBRARIES} libf2c.a)
