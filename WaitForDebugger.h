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
