dnl LB_CHECK_CDDLIB()
dnl
dnl Test for CDDLIB Library and define CDDLIB_CFLAGS and CDDLIB_LIBS

AC_DEFUN([LB_CHECK_CDDLIB],
[

AC_ARG_WITH([cddlib],
	[AS_HELP_STRING([--with-cddlib=auto|yes|no], [Make use of cddlib])],
	[want_cddlib="$withval"], [want_cddlib=auto])

have_cddlib=0
AS_IF([test "$want_cddlib" != no], [
	PKG_CHECK_MODULES([CDDLIB], [cddlib >= 0.94l], [
		have_cddlib=1
		AH_TEMPLATE([HAVE_CDDLIB], [Whether cddlib is available])
		AC_DEFINE([HAVE_CDDLIB], [1])
        ], [
		AS_IF([test "$want_cddlib" = yes], [AC_MSG_ERROR([$libcrypto_PKG_ERRORS])])
	])
])
AM_CONDITIONAL([HAVE_CDDLIB], [test "$have_cddlib" = 1])

])
