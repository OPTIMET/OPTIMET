#ifndef OPTIMET_BENCHMARK_CMDL_H
#define OPTIMET_BENCHMARK_CMDL_H

#include "Types.h"

#ifdef OPTIMET_BELOS
#include <Teuchos_ParameterList.hpp>

namespace optimet {
Teuchos::RCP<Teuchos::ParameterList> parse_cmdl(int argc, char *argv[]);
}
#endif
#endif
