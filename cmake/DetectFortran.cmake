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

# Detects fortran using the scalapack library
# Since we know what symbols to detect, we do not have to include the fortran compiler into the
# project, unlike the standard cmake approach.
# We do not try and dectect module mangling, since Scalapack does not use it.
function(DetectFortran OUTFILE)
  if(FORTRAN_MANGLING AND FORTRAN_MANGLING_)
    return()
  endif()
  if(NOT scalapack_FOUND AND "$ENV{CRAYOS_VERSION}" STREQUAL "")
    message(FATAL_ERROR "Please find package scalapack before calling this function")
  endif()
  file(MAKE_DIRECTORY "${PROJECT_BINARY_DIR}/DetectFortran")
  file(WRITE "${PROJECT_BINARY_DIR}/DetectFortran/main.cpp"
    "#include \"fcmangling.h\"\n"
    "extern \"C\" {\n"
    "void OPTIMET_FC_GLOBAL(pdgemr2d, PDGEMR2D)(int *m, int *n, double *A,\n"
    "            int *IA, int *JA, int *descA, double *B, int *IB, int *JB,\n"
    "            int *descB, int *gcontext);\n"
    "}\n"
    "int main(int nargs, char **argv) {\n"
    "  OPTIMET_FC_GLOBAL(pdgemr2d, PDGEMR2D)(nullptr, nullptr, nullptr,\n"
    "            nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr, nullptr);\n"
    "  return 0;\n"
    "}\n"
  )
  file(WRITE "${PROJECT_BINARY_DIR}/DetectFortran/main_underscore.cpp"
    "#include \"fcmangling.h\"\n"
    "extern \"C\" {\n"
    "void OPTIMET_FC_GLOBAL_(blacs_setup, BLACS_SETUP)(int *nprow, int *npcol);\n"
    "}\n"
    "int main(int nargs, char **argv) {\n"
    "  OPTIMET_FC_GLOBAL_(blacs_setup, BLACS_SETUP)(nullptr, nullptr);\n"
    "  return 0;\n"
    "}\n"
  )

  set(manglings "NAME" "name" "name ## _" "NAME ## _")
  foreach(mangling ${manglings})
      file(WRITE "${PROJECT_BINARY_DIR}/DetectFortran/fcmangling.h"
        "#ifndef OPTIMET_FC_MANGLING_H_\n"
        "#define OPTIMET_FC_MANGLING_H_\n"
        "#define OPTIMET_FC_GLOBAL(name, NAME) ${mangling}\n"
        "#endif\n"
      )
      try_compile(result
        "${PROJECT_BINARY_DIR}/DetectFortran"
        "${PROJECT_BINARY_DIR}/DetectFortran/main.cpp"
        LINK_LIBRARIES ${SCALAPACK_LIBRARIES} ${MPI_LIBRARIES}
        OUTPUT_VARIABLE output
      )
      if(result)
        set(found_mangling "${mangling}")
        message(STATUS "Dected fortran mangling without underscore: ${mangling}")
        break()
      endif()
  endforeach()
  if(NOT found_mangling)
    message(FATAL_ERROR "Could not detect fortran mangling")
  endif()

  foreach(underscore "NAME ## __" "name ## __" ${manglings})
      file(WRITE "${PROJECT_BINARY_DIR}/DetectFortran/fcmangling.h"
        "#ifndef OPTIMET_FC_MANGLING_H_\n"
        "#define OPTIMET_FC_MANGLING_H_\n"
        "#define OPTIMET_FC_GLOBAL_(name, NAME) ${underscore}\n"
        "#endif\n"
      )
      try_compile(result
        "${PROJECT_BINARY_DIR}/DetectFortran"
        "${PROJECT_BINARY_DIR}/DetectFortran/main_underscore.cpp"
        LINK_LIBRARIES ${SCALAPACK_LIBRARIES} ${MPI_LIBRARIES}
        OUTPUT_VARIABLE output
      )
      if(result)
        set(found_underscore_mangling "${underscore}")
        message(STATUS "Dected fortran mangling with underscore: ${underscore}")
        break()
      endif()
  endforeach()

  if(found_underscore_mangling AND found_mangling)
    file(WRITE "${OUTFILE}"
        "#ifndef OPTIMET_FC_MANGLING_H_\n"
        "#define OPTIMET_FC_MANGLING_H_\n"
        "#define OPTIMET_FC_GLOBAL(name, NAME) ${found_mangling}\n"
        "#define OPTIMET_FC_GLOBAL_(name, NAME) ${found_underscore_mangling}\n"
        "#endif\n"
    )
  set(FORTRAN_MANGLING "${found_mangling}" CACHE INTERNAL "Fortran mangling without underscore")
  set(FORTRAN_MANGLING_
    "${found_underscore_mangling}" CACHE INTERNAL "Fortran mangling with underscore")
  else()
    message(FATAL_ERROR "Could not detect fortran mangling")
  endif()
endfunction()
