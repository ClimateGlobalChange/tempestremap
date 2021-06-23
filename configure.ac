#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Usage: autoreconf -fi
#
AC_PREREQ([2.68])
AC_INIT([TempestRemap], [2.1.0], [paullrich@ucdavis.edu])
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
AX_CXX_COMPILE_STDCXX([11], [], [mandatory])

# Checks for some standard libraries.
AC_CHECK_LIB(dl, dlopen, LIBS="$LIBS -ldl")
AC_CHECK_LIB(m, pow, LIBS="$LIBS -lm")

# Checks for BLAS/LAPACK libraries:
AX_BLAS([], [AC_MSG_ERROR([BLAS library not found])])
AX_LAPACK([], [AC_MSG_ERROR([LAPACK library not found])])

# checking for netCDF
AC_LANG_PUSH([C++])
ACX_NETCDF([],[AC_MSG_ERROR(cannot find NetCDF library. Please specify the directory using --with-netcdf=DIR configure option.)])
AC_LANG_POP([C++])

# Checks for header files.
AC_CHECK_HEADERS([malloc.h stdlib.h string.h unistd.h])
AC_CHECK_HEADER(unordered_map, [AC_DEFINE([OVERLAPMESH_USE_UNSORTED_MAP], [1], [Use unordered maps in Overlap mesh],[])])

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

