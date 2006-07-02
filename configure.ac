#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(latte, 1.2-mk-0.3+cvs)
AC_CONFIG_HEADERS([code/config.h:code/config.hin])
AM_INIT_AUTOMAKE(gnu)
AM_MAINTAINER_MODE
AC_CONFIG_FILES([Makefile EXAMPLES/Makefile code/Makefile code/latte/Makefile code/test-suite/Makefile])
AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_CXX

AC_LANG(C++)

# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h stdlib.h string.h sys/param.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([pow sqrt])

AC_CHECK_HEADERS([sstream])

LB_MISC
LB_CHECK_GMP(,,[AC_MSG_ERROR([GMP not found!
GMP version 3.1.1 or greater is required for LattE to compile. Please
make sure GMP is installed and specify its location with the option
--with-gmp=<prefix> when running configure.])
])
LB_CHECK_NTL(,,[AC_MSG_ERROR([NTL not found!
NTL version 5.4 or greater is required for LattE to compile. Please
make sure NTL is installed and specify its location with the option
--with-ntl=<prefix> when running configure.])
])
LB_CHECK_LIDIA(,,
  AC_MSG_WARN([LiDIA not found.  We need LiDIA for computing Smith normal forms!]))

# From CDD+
AC_PATH_PROGS(CDD, [cddr+_gmp cddr+])
AC_DEFINE_UNQUOTED(CDD_PATH, "${CDD}", [The path to the CDD+ binary.])
if test -z "$CDD"; then
  AC_MSG_ERROR([The `cddr+' binary from CDD+ is required for LattE.])			
fi

# From cddlib
AC_PATH_PROGS(REDCHECK, [redcheck_gmp])
AC_DEFINE_UNQUOTED(REDCHECK_PATH, "${REDCHECK}", [The path to the cddlib redcheck binary.])
if test -z "$REDCHECK"; then
  AC_MSG_ERROR([The `redcheck' binary from cddlib is required for LattE.])			
fi
	  
# From LRS
AC_PATH_PROGS(LRS, [lrs1])
AC_DEFINE_UNQUOTED(LRS_PATH, "${LRS}", [The path to the lrs binary.])
if test -z "$LRS"; then
  AC_MSG_WARN([The `lrs1' binary from LRS is missing; LattE can optionally use LRS instead of CDD.])
fi

# We need Maple to generate the code for evaluation of the Todd polynomial. 
AC_PATH_PROG(MAPLEBIN, maple)
AM_CONDITIONAL(HAVE_MAPLE, test x$MAPLEBIN != x )

AC_ARG_ENABLE(experiments, AS_HELP_STRING([--enable-experiments],
	[use experimental code]))
AM_CONDITIONAL(HAVE_EXPERIMENTS, test x$enable_experiments = xyes)
if test x$enable_experiments = xyes ; then
   AC_DEFINE(HAVE_EXPERIMENTS, 1, [Define if experimental code is enabled])
fi

AC_OUTPUT
