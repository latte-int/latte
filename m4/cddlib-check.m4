# Check for cddlib, derived from lidia-check.m4

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
					   search path (/usr or /usr/local).
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
BACKUP_LIBS=${LIBS}

if test -n "$CDDLIB_HOME_PATH" ; then
AC_MSG_CHECKING(for CDDLIB >= $min_cddlib_version)
fi

for CDDLIB_HOME in ${CDDLIB_HOME_PATH} 
 do	
if test -r "$CDDLIB_HOME/include/cdd.h"; then
	if test "x$CDDLIB_HOME" != "x/usr" -a "x$CDDLIB_HOME" != "x/usr/local"; then
		CDDLIB_CFLAGS="-I${CDDLIB_HOME}/include"
		CDDLIB_LIBS="-L${CDDLIB_HOME}/lib -lcddgmp"
	else
		CDDLIB_CFLAGS=
		CDDLIB_LIBS="-lcddgmp"		
	fi	
	CXXFLAGS="${BACKUP_CXXFLAGS} ${CDDLIB_CFLAGS} ${GMP_CFLAGS}" 
	LIBS="${BACKUP_LIBS} ${CDDLIB_LIBS} ${GMP_LIBS}"

	AC_TRY_LINK([
#define GMPRATIONAL
#include <setoper.h>
#include <cddmp.h>
],
[ mytype a;
  dd_init(a);
],
[	cddlib_found="yes"
	break
]
)
else
	cddlib_found="no"
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
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
