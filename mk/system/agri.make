# Copyright (c) 2016      Bryce Adelstein-Lelbach aka wash
# Copyright (c) 2000-2016 Paul Ullrich 
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying 
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

# UC Davis Agri 

CXX=               g++
F90=               gfortran
MPICXX=            mpiCC
MPIF90=            mpif90

CXXFLAGS+=         -fPIC -Wno-literal-suffix
F90FLAGS+=         -fPIC

F90_RUNTIME=       -lgfortran

# NetCDF
NETCDF_ROOT=       /opt/local
NETCDF_CXXFLAGS=   -I$(NETCDF_ROOT)/include
NETCDF_LIBRARIES=  -lnetcdf -lnetcdf_c++
NETCDF_LDFLAGS=    -L$(NETCDF_ROOT)/lib

# PetSc
PETSC_ROOT=        /opt/local/lib/petsc
X11_ROOT=          /opt/X11
PETSC_CXXFLAGS=    -I$(PETSC_ROOT)/include
PETSC_LIBRARIES=   -lpetsc -lX11
PETSC_LDFLAGS=     -L$(PETSC_ROOT)/lib -L$(X11_ROOT)/lib

# LAPACK (Mac OS X Accelerate Framework)
LAPACK_INTERFACE=  FORTRAN
LAPACK_CXXFLAGS=
LAPACK_LIBRARIES=  -llapack -lblas
LAPACK_LDFLAGS=    

# DO NOT DELETE
