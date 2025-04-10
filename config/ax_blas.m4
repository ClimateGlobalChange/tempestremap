# ===========================================================================
#         https://www.gnu.org/software/autoconf-archive/ax_blas.html
# ===========================================================================
#
# SYNOPSIS
#
#   AX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a library that implements the BLAS linear-algebra
#   interface (see http://www.netlib.org/blas/). On success, it sets the
#   BLAS_LIBS output variable to hold the requisite library linkages.
#
#   To link with BLAS, you should link with:
#
#     $BLAS_LIBS $LIBS $FLIBS
#
#   in that order. FLIBS is the output variable of the
#   AC_F77_LIBRARY_LDFLAGS macro (called if necessary by AX_BLAS), and is
#   sometimes necessary in order to link with F77 libraries. Users will also
#   need to use AC_F77_DUMMY_MAIN (see the autoconf manual), for the same
#   reason.
#
#   Many libraries are searched for, from ATLAS to CXML to ESSL. The user
#   may also use --with-blas=<lib> in order to use some specific BLAS
#   library <lib>. In order to link successfully, however, be aware that you
#   will probably need to use the same Fortran compiler (which can be set
#   via the F77 env. var.) as was used to compile the BLAS library.
#
#   ACTION-IF-FOUND is a list of shell commands to run if a BLAS library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_BLAS.
#
# LICENSE
#
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#   Copyright (c) 2019 Geoffrey M. Oxberry <goxberry@gmail.com>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <https://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Archive. When you make and distribute a
#   modified version of the Autoconf Macro, you may extend this special
#   exception to the GPL to apply to your modified version as well.

#serial 21

AU_ALIAS([ACX_BLAS], [AX_BLAS])
AC_DEFUN([AX_BLAS], [
AC_PREREQ([2.55])
AC_REQUIRE([AC_CANONICAL_HOST])
ax_blas_ok=no

AC_ARG_WITH(blas,
	[AS_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) ax_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.dylib | *.dylib.* | *.o)
		BLAS_LIBS="$with_blas"
	;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

# Get fortran linker names of BLAS functions to check for.
#AC_F77_FUNC(daxpy)
#AC_F77_FUNC(dgemm)
#daxpy="daxpy DAXPY daxpy_ daxpy__ _daxpy DAXPY_ _DAXPY"
#dgemm="dgemm DGEMM dgemm_ dgemm__ _dgemm DGEMM_ _DGEMM"
#daxpy="daxpy_"
#dgemm="dgemm_"

ax_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $ax_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $daxpy in $BLAS_LIBS])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$daxpy])], [ax_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($ax_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $ax_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_MSG_CHECKING([if $daxpy is being linked in already])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$daxpy])], [ax_blas_ok=yes])
	AC_MSG_RESULT($ax_blas_ok)
	LIBS="$save_LIBS"
fi

# BLAS linked to by flexiblas? (Default on FC33+ and RHEL9+)

if test $ax_blas_ok = no; then
	AC_CHECK_LIB(flexiblas, $daxpy, [ax_blas_ok=yes
			                BLAS_LIBS="-lflexiblas"])
fi

# BLAS in OpenBLAS library? (http://xianyi.github.com/OpenBLAS/)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(openblas, $daxpy, [ax_blas_ok=yes
			                BLAS_LIBS="-lopenblas"])
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $daxpy,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[ax_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $daxpy,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(daxpy, $daxpy,
			[ax_blas_ok=yes; BLAS_LIBS="-ldaxpy -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Intel MKL library?
if test $ax_blas_ok = no; then
	# MKL for gfortran
	if test x"$ac_cv_cxx_compiler_gnu" = xyes; then
		# 64 bit
		if test $host_cpu = x86_64; then
			AC_CHECK_LIB(mkl_gf_lp64, $daxpy,
			[ax_blas_ok=yes;BLAS_LIBS="-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread"],,
			[-lmkl_gf_lp64 -lmkl_sequential -lmkl_core -lpthread])
		# 32 bit
		elif test $host_cpu = i686; then
			AC_CHECK_LIB(mkl_gf, $daxpy,
				[ax_blas_ok=yes;BLAS_LIBS="-lmkl_gf -lmkl_sequential -lmkl_core -lpthread"],,
				[-lmkl_gf -lmkl_sequential -lmkl_core -lpthread])
		fi
	# MKL for other compilers (Intel, PGI, ...?)
	else
		# 64-bit
		if test $host_cpu = x86_64; then
			AC_CHECK_LIB(mkl_intel_lp64, $daxpy,
				[ax_blas_ok=yes;BLAS_LIBS="-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -limf -lpthread"],,
				[-lmkl_intel_lp64 -lmkl_sequential -lmkl_core -limf -lpthread])
		# 32-bit
		elif test $host_cpu = i686; then
			AC_CHECK_LIB(mkl_intel, $daxpy,
				[ax_blas_ok=yes;BLAS_LIBS="-lmkl_intel -lmkl_sequential -lmkl_core -limf -lpthread"],,
				[-lmkl_intel -lmkl_sequential -lmkl_core -limf -lpthread])
		fi
	fi
fi
# Old versions of MKL
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(mkl, $daxpy, [ax_blas_ok=yes;BLAS_LIBS="-lmkl -lguide -lpthread"],,[-lguide -lpthread])
fi

# BLAS in Apple vecLib library?
if test $ax_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="-framework Accelerate $LIBS"
	AC_MSG_CHECKING([for $daxpy in -framework Accelerate])
	AC_LINK_IFELSE([AC_LANG_CALL([], [$daxpy])], [ax_blas_ok=yes;BLAS_LIBS="-framework Accelerate"])
	AC_MSG_RESULT($ax_blas_ok)
	LIBS="$save_LIBS"
fi

# BLAS in Alpha CXML library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(cxml, $daxpy, [ax_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(dxml, $daxpy, [ax_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $ax_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $daxpy,
				[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 ax_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(scs, $daxpy, [ax_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $daxpy,
		     [ax_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $daxpy,
		[AC_CHECK_LIB(essl, $daxpy,
			[ax_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $daxpy, [ax_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

AC_SUBST(BLAS_LIBS)

LIBS="$ax_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_blas_ok" = xyes; then
        ifelse([$1],,AC_DEFINE(HAVE_BLAS,1,[Define if you have a BLAS library.]),[$1])
        :
else
        ax_blas_ok=no
        $2
fi
])dnl AX_BLAS
