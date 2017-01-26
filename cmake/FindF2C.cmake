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

find_path(F2C_INCLUDE_DIR NAMES f2c.h)
find_library(F2C_LIBRARY NAMES libf2c.a f2c)

set(F2C_INCLUDE_DIRS ${F2C_INCLUDE_DIR})
set(F2C_LIBRARIES ${F2C_LIBRARY})

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(F2C REQUIRED_VARS F2C_INCLUDE_DIR F2C_LIBRARY)
