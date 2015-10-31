find_path(Hdf5_INCLUDE_DIR
  hdf5.h
  HINTS $ENV{HDF5_ROOT}
  PATH_SUFFIXES hdf5 hdf5/openmpi hdf5/serial
)
find_library(Hdf5_LIBRARY
  hdf5
  HINTS $ENV{HDF5_ROOT}
  PATH_SUFFIXES hdf5 hdf5/openmpi hdf5/serial
)
