if(ARMA_LIBRARIES)
    set(ARMA_FIND_QUIETLY TRUE)
endif(ARMA_LIBRARIES)

find_library(ARMA_LIB
    armadillo
    HINTS
	/usr/lib
    )

if (ARMA_LIB)
    set(ARMA_LIBRARIES ${ARMA_LIB} m)
endif (ARMA_LIB)

find_path(ARMA_INCLUDE_DIRS
    NAMES armadillo
    HINTS
	/usr/include
    )

if(ARMA_LIBRARIES AND ARMA_INCLUDE_DIRS)
    set(ARMA_FOUND TRUE)
endif(ARMA_LIBRARIES AND ARMA_INCLUDE_DIRS)
