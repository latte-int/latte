AC_INIT(latte, 1.2-mk-0.1)
AC_CONFIG_HEADERS([code/config.h:code/config.hin])
AM_INIT_AUTOMAKE(foreign)
AM_MAINTAINER_MODE
AC_CONFIG_FILES([Makefile code/Makefile code/latte/Makefile])
AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_CXX
AC_LANG(C++)

LB_MISC
LB_CHECK_GMP(,,[AC_MSG_ERROR([GMP not found!
GMP version 3.1.1 or greater is required for LattE to compile. Please
make sure GMP is installed and specify its location with the option
--with-gmp=<prefix> when running configure.])
])
LB_CHECK_NTL(,,[AC_MSG_ERROR([NTL not found!
NTL version 5.0 or greater is required for LattE to compile. Please
make sure NTL is installed and specify its location with the option
--with-ntl=<prefix> when running configure.])
])
LB_CHECK_LIDIA(,,
  AC_MSG_WARN([LiDIA not found.  We need LiDIA for computing Smith normal forms!]))

# From CDD+
AC_PATH_PROGS(CDD, [cddr+])
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

AC_OUTPUT
