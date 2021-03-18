# Check for cddlib, derived from lidia-check.m4

dnl LB_CDDLIB_TEST_HEADERS  ([header_prefix, path])

AC_DEFUN([LB_CDDLIB_TEST_HEADERS], [
cdd_header_found=no
    AC_DEFINE_UNQUOTED(CDDLIB_SETOPER_H, [<$1setoper.h>], [header setoper.h])
    AC_DEFINE_UNQUOTED(CDDLIB_CDDMP_H, [<$1cddmp.h>], [header cddmp.h])
    AC_DEFINE_UNQUOTED(CDDLIB_CDD_H, [<$1cdd.h>], [header cdd.h])
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

AC_DEFUN([LB_CDDLIB_TEST_HEADERS_ALL], [
    cddlib_found=no
    dnl cddlib is the prefix used by modern cddlib (0.94m).
    for cdd_h_path in "cddlib/" "" "cdd/"
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

dnl LB_CHECK_CDDLIB ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for CDDLIB Library and define CDDLIB_CFLAGS and CDDLIB_LIBS

AC_DEFUN([LB_CHECK_CDDLIB],
[

CDDLIB_HOME_PATH="${DEFAULT_CHECKING_PATH}"

AC_ARG_WITH(cddlib,
	    [  --with-cddlib=<path>|yes|no 
					   Use cddlib. 
					   If argument is no, you do not have the library installed on your machine (set as default).
					   If argument is yes or <empty> that means the library is reachable with the standard
					   search path.
	 				   Otherwise you give the <path> to the directory which contain the library. 
	     ],
	     [if test "$withval" = yes ; then
			CDDLIB_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	      elif test "$withval" != no ; then
			CDDLIB_HOME_PATH="$withval"
	     fi],
	     [])

min_cddlib_version=ifelse([$1], , 093c,$1)


dnl Check for existence

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_CFLAGS=${CFLAGS}
BACKUP_LIBS=${LIBS}

if test -n "$CDDLIB_HOME_PATH" ; then
AC_MSG_CHECKING(for CDDLIB >= $min_cddlib_version)
fi

for CDDLIB_HOME in ${CDDLIB_HOME_PATH} 
 do	
	if test "$CDDLIB_HOME" != "DEFAULTS" ; then
		CDDLIB_CFLAGS="-I${CDDLIB_HOME}/include"
		CDDLIB_LIBS="-L${CDDLIB_HOME}/lib -lcddgmp"
	else
		CDDLIB_CFLAGS=
		CDDLIB_LIBS="-lcddgmp"		
	fi	
	CXXFLAGS="${CDDLIB_CFLAGS} ${GMP_CFLAGS} ${BACKUP_CXXFLAGS} "
	CFLAGS="${CDDLIB_CFLAGS} ${BACKUP_CFLAGS} {GMP_CFLAGS}"
	LIBS="${CDDLIB_LIBS} ${GMP_LIBS} ${BACKUP_LIBS}"

	LB_CDDLIB_TEST_HEADERS_ALL
	if test "x$cddlib_found" = "xyes" ; then
	    break
	fi
done

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
