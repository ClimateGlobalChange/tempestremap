#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Usage: autoreconf -fi
#
AC_PREREQ([2.68])
AC_INIT([TempestRemap], [1.0-beta], [paullrich@ucdavis.edu])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
#AC_CONFIG_SRCDIR([test/optimtest.m])
AC_CONFIG_HEADERS([src/TempestConfig.h])

AM_INIT_AUTOMAKE([foreign subdir-objects])
AM_PROG_LIBTOOL
AM_SILENT_RULES([yes])

LT_INIT([disable-shared])
LT_LANG([C++])

# Checks for programs.
AC_PROG_CXX([icpc g++])
AC_PROG_CC([icc gcc])
AC_PROG_MAKE_SET

DISTCHECK_CONFIGURE_FLAGS="$DISTCHECK_CONFIGURE_FLAGS CC=$CC CXX=$CXX FC=$FC"

# Check how to add C++0x flags
# TODO: Should we check for partial C++11 support also ?
AX_CXX_COMPILE_STDCXX_0X

# Checks for FORTRAN.
AC_PROG_F77([ifort gfortran g77 f77])    
AC_F77_LIBRARY_LDFLAGS
LIBS="$LIBS $FCLIBS -lm"

# Checks for BLAS/LAPACK libraries:
AX_BLAS([], [AC_MSG_ERROR([BLAS library not found])])
AX_LAPACK([], [AC_MSG_ERROR([LAPACK library not found])])

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lnetcdf':
# checking for netCDF
AC_LANG_PUSH([C++])
ACX_NETCDF([],[AC_MSG_ERROR(cannot find netCDF fortran library)])
AC_LANG_POP([C++])

# Checks for header files.
AC_CHECK_HEADERS([malloc.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_SIZE_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_ALLOCA
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([memset pow sqrt])

AC_SUBST(DISTCHECK_CONFIGURE_FLAGS)
AC_CONFIG_FILES([Makefile])
AC_OUTPUT

