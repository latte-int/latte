AC_INIT(latte, 1.2)
AC_CONFIG_HEADERS([code/config.h:code/config.hin])
AM_INIT_AUTOMAKE(foreign)
AM_MAINTAINER_MODE
AC_CONFIG_FILES([Makefile code/Makefile code/latte/Makefile])
AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_CXX
AC_LANG(C++)

LB_CHECK_NTL

# From CDD+
AC_PATH_PROGS(CDD, [cddr+])
AC_DEFINE_UNQUOTED(CDD_PATH, "${CDD}", [The path to the CDD+ binary.])

# From cddlib
AC_PATH_PROGS(REDCHECK, [redcheck_gmp])
AC_DEFINE_UNQUOTED(REDCHECK_PATH, "${REDCHECK}", [The path to the cddlib redcheck binary.])

# From LRS
AC_PATH_PROGS(LRS, [lrs1])
AC_DEFINE_UNQUOTED(LRS_PATH, "${LRS}", [The path to the lrs binary.])

AC_OUTPUT
