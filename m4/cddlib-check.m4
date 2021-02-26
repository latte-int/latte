# Check for cddlib, derived from lidia-check.m4

dnl LB_CDDLIB_TEST_HEADERS  ([header_prefix, path])

AC_DEFUN([LB_CDDLIB_TEST_HEADERS], [
cdd_header_found=no
AS_IF([test "x$2" = "x"],
[
    AC_DEFINE_UNQUOTED(CDDLIB_SETOPER_H, [<$1setoper.h>], [header setoper.h])
    AC_DEFINE_UNQUOTED(CDDLIB_CDDMP_H, [<$1cddmp.h>], [header cddmp.h])
    AC_DEFINE_UNQUOTED(CDDLIB_CDD_H, [<$1cdd.h>], [header cdd.h])
],[
    # Use the specified path and prevent other installs from leaking in.
    AC_DEFINE_UNQUOTED(CDDLIB_SETOPER_H, ["$2/include/$1setoper.h"], [header setoper.h])
    AC_DEFINE_UNQUOTED(CDDLIB_CDDMP_H, ["$2/include/$1cddmp.h"], [header cddmp.h])
    AC_DEFINE_UNQUOTED(CDDLIB_CDD_H, ["$2/include/$1cdd.h"], [header cdd.h])
])
AC_TRY_LINK([
#define GMPRATIONAL
#include CDDLIB_SETOPER_H
#include CDDLIB_CDDMP_H
#include CDDLIB_CDD_H
],
[ mytype a;
  dd_init(a);
  dd_abs(a, a);
],
[       cdd_header_found=yes
]
)
])

dnl LB_CHECK_CDDLIB ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for CDDLIB Library and define CDDLIB_CFLAGS and CDDLIB_LIBS

AC_DEFUN([LB_CHECK_CDDLIB],
[


want_cddlib=auto

AC_ARG_WITH(cddlib,
            [  --with-cddlib=<path>|yes|no
                                           Use cddlib.
                                           If argument is no, you do not have the library installed on your machine (set as default).
                                           If argument is yes or <empty> that means the library is reachable with the standard
                                           search path (/usr or /usr/local).
                                           Otherwise you give the <path> to the directory which contain the library.
            ],
            [if test "$withval" = no ; then
                want_cddlib=no
            elif test "$withval" = yes ; then
                want_cddlib=yes
            else
                CDDLIB_HOME="$withval"
            fi],
            [])

min_cddlib_version=ifelse([$1], , 0.93c,$1)

cddlib_found=no
AS_IF([test "x$want_cddlib" != "xno"],
[
    AC_MSG_CHECKING(for CDDLIB >= $min_cddlib_version)
    AS_IF([test -n "$CDDLIB_HOME"],
    [
        cddlib_found=yes
        CDDLIB_CFLAGS="-I${CDDLIB_HOME}/include"
        CDDLIB_LIBS="-L${CDDLIB_HOME}/lib -lcddgmp"
    ], [
        PKG_CHECK_MODULES([CDDLIB], [cddlib >= "$min_cddlib_version"],
            [cddlib_found=yes],
            [AS_IF([test "$want_cddlib" = yes], [AC_MSG_ERROR([$libcrypto_PKG_ERRORS])])])
    ])
])

dnl Check for the correct header files.
dnl E.g. <cdd.h> or <cdd/cdd.h> or <cddlib/cdd.h>.

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_CFLAGS=${CFLAGS}
BACKUP_LIBS=${LIBS}

CXXFLAGS="${BACKUP_CXXFLAGS} ${CDDLIB_CFLAGS} ${GMP_CFLAGS}"
CFLAGS="${BACKUP_CFLAGS} ${CDDLIB_CFLAGS} ${GMP_CFLAGS}"
LIBS="${BACKUP_LIBS} ${CDDLIB_LIBS} ${GMP_LIBS}"

AS_IF([test "x$cddlib_found" = "xyes"],
[
    AC_MSG_NOTICE([found cdd, trying to find the headers])
    cddlib_found=no
    for cdd_h_path in "" "cdd/" "cddlib/"
        do
            AS_IF([test -z $cdd_h_path], [cdd_h_path_print="inclusion directory"], [cdd_h_path_print="${cdd_h_path}"])
            LB_CDDLIB_TEST_HEADERS([${cdd_h_path}], [${CDDLIB_HOME}])
            AS_IF([test "x$cdd_header_found" = "xyes"],
            [
                AC_MSG_NOTICE([headers are in ${cdd_h_path_print}])
                cddlib_found=yes
                break
            ],
            [
                AC_MSG_NOTICE([headers are not in ${cdd_h_path_print}])
            ])
    done
])

if test "x$cddlib_found" = "xyes" ; then
        AC_SUBST(CDDLIB_CFLAGS)
        AC_SUBST(CDDLIB_LIBS)
        AC_DEFINE(HAVE_CDDLIB,1,[Define if CDDLIB is installed])
        HAVE_CDDLIB=yes
        AC_MSG_RESULT(found)
else
        AC_MSG_RESULT(not found)
        ifelse([$3], , :, [$3])
fi

AM_CONDITIONAL(HAVE_CDDLIB, test "x$HAVE_CDDLIB" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
CFLAGS=${BACKUP_CFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
