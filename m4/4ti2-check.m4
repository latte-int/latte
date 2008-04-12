# Check for 4ti2 library, derived from cddlib-check.m4

dnl LB_CHECK_FORTYTWO ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for FORTYTWO Library and define FORTYTWO_CXXFLAGS and FORTYTWO_LIBS

AC_DEFUN([LB_CHECK_FORTYTWO],
[

FORTYTWO_HOME_PATH="${DEFAULT_CHECKING_PATH}"

AC_ARG_WITH(4ti2,
	    [  --with-4ti2=<path>|yes|no 
					   Use 4ti2. 
					   If argument is no, you do not have the library installed on your machine (set as default).
					   If argument is yes or <empty> that means the library is reachable with the standard
					   search path (/usr or /usr/local).
	 				   Otherwise you give the <path> to the directory which contain the library. 
	     ],
	     [if test "$withval" = yes ; then
			FORTYTWO_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	      elif test "$withval" != no ; then
			FORTYTWO_HOME_PATH="$withval"
	     fi],
	     [])

dnl Check for existence

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

if test -n "$FORTYTWO_HOME_PATH" ; then
AC_MSG_CHECKING(for 4ti2 library)
fi

for FORTYTWO_HOME in ${FORTYTWO_HOME_PATH} 
 do	
    if test -r "$FORTYTWO_HOME/include/groebner/RayAlgorithm.h"; then
	if test "x$FORTYTWO_HOME" != "x/usr" -a "x$FORTYTWO_HOME" != "x/usr/local"; then
		FORTYTWO_CXXFLAGS="-I${FORTYTWO_HOME}/include -D__STDC_LIMIT_MACROS -D_4ti2_GMP_"
		FORTYTWO_LIBS="-L${FORTYTWO_HOME}/lib -l4ti2gmp -lzsolve"
	else
		FORTYTWO_CXXFLAGS="-D__STDC_LIMIT_MACROS -D_4ti2_GMP_"
		FORTYTWO_LIBS="-l4ti2gmp -lzsolve"
	fi	
	CXXFLAGS="${BACKUP_CXXFLAGS} ${FORTYTWO_CXXFLAGS} ${GMP_CFLAGS}" 
	LIBS="${BACKUP_LIBS} ${FORTYTWO_LIBS} ${GMP_LIBS}"

	AC_TRY_LINK([
#include "groebner/RayAlgorithm.h"
#include "zsolve/VectorArray.hpp"
],
[ _4ti2_::RayAlgorithm algorithm;
  _4ti2_zsolve::VectorArray<int> array; 
],
[	FORTYTWO_found="yes"
	break
])
fi
done

if test "x$FORTYTWO_found" = "xyes" ; then		
	AC_SUBST(FORTYTWO_CXXFLAGS)
	AC_SUBST(FORTYTWO_LIBS)
	AC_DEFINE(HAVE_FORTYTWO_LIB,1,[Define if the 4ti2 library is installed])
	HAVE_FORTYTWO_LIB=yes
	AC_MSG_RESULT(found)
else
	AC_MSG_RESULT(not found)
	ifelse([$3], , :, [$3])
fi	

AM_CONDITIONAL(HAVE_FORTYTWO_LIB, test "x$HAVE_FORTYTWO_LIB" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

])
