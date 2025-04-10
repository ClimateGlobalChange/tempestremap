# Copyright (c) 2025 Hongyu Chen

CXX=               CC
F90=               ftn
MPICXX=            CC
MPIF90=            ftn

# NetCDF
NETCDF_ROOT=       $(NETCDF_DIR)
NETCDF_CXXFLAGS=   -I$(NETCDF_ROOT)/include
NETCDF_LIBRARIES=  -lnetcdf
NETCDF_LDFLAGS=    -L$(NETCDF_ROOT)/lib

# HDF5
HDF5_ROOT=         $(HDF5_DIR)
CPPFLAGS +=        -I$(HDF5_ROOT)/include
LDFLAGS  +=        -L$(HDF5_ROOT)/lib
LDLIBS   +=        -lhdf5_hl -lhdf5

# LAPACK (optional: most likely not used)
LAPACK_INTERFACE=  FORTRAN
LAPACK_CXXFLAGS=
LAPACK_LIBRARIES=
LAPACK_LDFLAGS=

# DO NOT DELETE
