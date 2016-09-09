#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.
# Usage: autoreconf -fi
#
AC_PREREQ([2.69])
AC_INIT([TempestRemap], [1.0-beta], [paullrich@ucdavis.edu])
AC_CONFIG_AUX_DIR([config])
AC_CONFIG_MACRO_DIR([config])
#AC_CONFIG_SRCDIR([test/optimtest.m])
AC_CONFIG_HEADERS([src/TempestConfig.h])

AM_INIT_AUTOMAKE([subdir-objects])
LT_INIT

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_MAKE_SET

# Check how to add C++0x flags
# TODO: Should we check for partial C++11 support also ?
AX_CXX_COMPILE_STDCXX_0X

# Checks for some standard libraries
AC_CHECK_LIB([m],[cos])
AC_CHECK_LIB([dl], [dlopen])

# Checks for BLAS/LAPACK libraries:
AX_BLAS([], [AC_MSG_ERROR([BLAS library not found])])
AX_LAPACK([], [AC_MSG_ERROR([LAPACK library not found])])

# Checks for libraries.
# FIXME: Replace `main' with a function in `-lnetcdf':
# checking for netCDF
AC_LANG_PUSH(C++)
ACX_NETCDF([],[AC_MSG_ERROR(cannot find NetCDF C++ library)])
AC_LANG_POP(C++)

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

AC_CONFIG_FILES([Makefile])
AC_OUTPUT

