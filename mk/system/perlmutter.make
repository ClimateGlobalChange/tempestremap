# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# NERSC Babbage Testbed

CXX=               CC
F90=               ftn
MPICXX=            CC
MPIF90=            ftn

# NetCDF
NETCDF_ROOT=       $(NETCDF_DIR)
NETCDF_CXXFLAGS=   -I$(NETCDF_ROOT)/include
NETCDF_LIBRARIES=  -lnetcdf
NETCDF_LDFLAGS=    -L$(NETCDF_ROOT)/lib

# LAPACK (Intel MKL)
LAPACK_INTERFACE=  FORTRAN
LAPACK_CXXFLAGS=
LAPACK_LIBRARIES=
LAPACK_LDFLAGS=   

# DO NOT DELETE
