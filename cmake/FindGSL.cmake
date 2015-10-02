if (GSL_LIBRARIES)
    set(GSL_FIND_QUIETLY TRUE)
endif (GSL_LIBRARIES)

find_library(GSL_LIB
    gsl
    HINTS
	/usr/lib
    )
    
find_library(GSL_CBLAS_LIB
	gslcblas
	HINTS
	/usr/lib
	)

if(GSL_LIB)
    set(GSL_LIBRARIES ${GSL_LIB} ${GSL_CBLAS_LIB})
endif(GSL_LIB)

find_path(GSL_INCLUDE_DIRS
    NAMES gsl/gsl_sf.h
    HINTS
	/usr/include
    )

if(GSL_LIBRARIES AND GSL_INCLUDE_DIRS)
    set(GSL_FOUND TRUE)
endif(GSL_LIBRARIES AND GSL_INCLUDE_DIRS)
