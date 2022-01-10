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

#ifndef OPTIMET_ATTACH_DEBUG_H
#define OPTIMET_ATTACH_DEBUG_H

#include <unistd.h>
#include <iostream>

namespace optimet {
void wait_for_debugger()
{
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    std::cout << "PID " << getpid() << " on " << hostname << " ready for attach" << std::endl;
    while (0 == i)
        sleep(5);
}
}
#endif
