if(HDF5_LIBRARIES)
    set(HDF5_FIND_QUIETLY TRUE)
endif(HDF5_LIBRARIES)

find_library(HDF5_LIB
    hdf5
    HINTS
	/usr/lib/i386-linux-gnu
    )

if (HDF5_LIB)
    set(HDF5_LIBRARIES ${HDF5_LIB})
endif (HDF5_LIB)

find_path(HDF5_INCLUDE_DIRS
    NAMES hdf5.h
    HINTS
	/usr/include
    )

if(HDF5_LIBRARIES AND HDF5_INCLUDE_DIRS)
    set(HDF5_FOUND TRUE)
endif(HDF5_LIBRARIES AND HDF5_INCLUDE_DIRS)
