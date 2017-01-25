// (C) University College London 2017
// This file is part of Optimet, licensed under the terms of the GNU Public License
//
// Optimet is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Optimet is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Optimet. If not, see <http://www.gnu.org/licenses/>.

#ifndef READER_H_
#define READER_H_

#include "Run.h"
#include "pugi/pugixml.hpp"
#include <string>

#ifdef OPTIMET_BELOS
#include <Teuchos_ParameterList.hpp>
#endif

namespace optimet {
//! Reads simulation configuration from input
Run simulation_input(std::string const &filename);
//! Reads simulation configuration from string buffer
Run simulation_input(std::istream &buffer);
}

#endif /* READER_H_ */
