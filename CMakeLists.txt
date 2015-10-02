cmake_minimum_required(VERSION 2.6)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/cmake)

set(PROJECT "Optimet3D")

file(GLOB SRC ${CMAKE_SOURCE_DIR}/*.c* ${CMAKE_SOURCE_DIR}/pugi/*.c*)
file(GLOB HPP ${CMAKE_SOURCE_DIR}/*.h* ${CMAKE_SOURCE_DIR}/*.h*)

#Libraries and headers NOT included as part ofa package
set(STATIC_LIBS ${CMAKE_SOURCE_DIR}/libf2c.a)
set(DYNAMIC_LIBS "")
set(INCLUDES "")

project(${PROJECT})

include_directories(${INCLUDES})
add_executable(${PROJECT} ${SRC})
target_link_libraries(${PROJECT} ${STATIC_LIBS})
target_link_libraries(${PROJECT} ${DYNAMIC_LIBS})

#Find extra packages (Armadillo, HDF5, GSL)

find_package(HDF5 REQUIRED)
if(HDF5_FOUND)
	include_directories(${HDF5_INCLUDE_DIRS})
	target_link_libraries(${PROJECT} ${HDF5_LIBRARIES})
endif(HDF5_FOUND)

find_package(ARMA REQUIRED)
if(ARMA_FOUND)
	include_directories(${ARMA_INCLUDE_DIRS})
	target_link_libraries(${PROJECT} ${ARMA_LIBRARIES})
endif(ARMA_FOUND)

find_package(GSL REQUIRED)
if(GSL_FOUND)
	include_directories(${GSL_INCLUDE_DIRS})
	target_link_libraries(${PROJECT} ${GSL_LIBRARIES})
endif(GSL_FOUND)
