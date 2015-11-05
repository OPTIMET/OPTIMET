cmake_minimum_required(VERSION 3.2)

project(Optimet3D)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

find_package(GSL REQUIRED)
find_package(HDF5 REQUIRED)
find_package(Armadillo REQUIRED)
find_package(F2C REQUIRED)

include_directories(${GSL_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS} ${F2C_INCLUDE_DIRS})

file(GLOB SRC ${CMAKE_SOURCE_DIR}/*.c* ${CMAKE_SOURCE_DIR}/pugi/*.c*)

add_executable(Optimet3D ${SRC})
target_link_libraries(Optimet3D ${GSL_LIBRARIES} ${HDF5_LIBRARIES} ${ARMADILLO_LIBRARIES} ${F2C_LIBRARIES})
