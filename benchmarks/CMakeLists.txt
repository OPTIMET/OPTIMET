# (C) University College London 2017
# This file is part of Optimet, licensed under the terms of the GNU Public License
#
# Optimet is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# Optimet is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Optimet. If not, see <http:#www.gnu.org/licenses/>.

if(OPTIMET_BELOS)
  add_benchmark(benchmark benchmarks.cpp fcc.cpp cmdl.cpp LIBRARIES optilib ${library_dependencies})

endif()
if(OPTIMET_MPI)
  add_executable(parallel_fmm_multiplication fmm_multiplication.cpp)
  target_link_libraries(parallel_fmm_multiplication optilib ${library_dependencies})

	add_executable(parallel_scalapack_multiplication scalapack_multiplication.cpp)
	target_link_libraries(parallel_scalapack_multiplication optilib ${library_dependencies})

	add_executable(parallel_solver solver.cpp cmdl.cpp)
	target_link_libraries(parallel_solver optilib ${library_dependencies})
endif()

add_executable(serial_fmm_multiplication fmm_multiplication.cpp)
target_link_libraries(serial_fmm_multiplication optilib ${library_dependencies})
target_compile_definitions(serial_fmm_multiplication PRIVATE OPTIMET_JUST_DO_SERIAL)
