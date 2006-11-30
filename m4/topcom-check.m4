# Check for TOPCOM library, derived from cddlib-check.m4

dnl LB_CHECK_TOPCOM ([MINIMUM-VERSION [, ACTION-IF-FOUND [, ACTION-IF-NOT-FOUND]]])
dnl
dnl Test for TOPCOM Library and define TOPCOM_CXXFLAGS and TOPCOM_LIBS

AC_DEFUN([LB_CHECK_TOPCOM],
[

TOPCOM_HOME_PATH="${DEFAULT_CHECKING_PATH}"

AC_ARG_WITH(topcom,
	    [  --with-topcom=<path>|yes|no 
					   Use TOPCOM. 
					   If argument is no, you do not have the library installed on your machine (set as default).
					   If argument is yes or <empty> that means the library is reachable with the standard
					   search path (/usr or /usr/local).
	 				   Otherwise you give the <path> to the directory which contain the library. 
	     ],
	     [if test "$withval" = yes ; then
			TOPCOM_HOME_PATH="${DEFAULT_CHECKING_PATH}"
	      elif test "$withval" != no ; then
			TOPCOM_HOME_PATH="$withval"
	     fi],
	     [])

dnl Check for existence

BACKUP_CXXFLAGS=${CXXFLAGS}
BACKUP_LIBS=${LIBS}

if test -n "$TOPCOM_HOME_PATH" ; then
AC_MSG_CHECKING(for TOPCOM library)
fi

for TOPCOM_HOME in ${TOPCOM_HOME_PATH} 
 do	
if test -r "$TOPCOM_HOME/include/PointConfiguration.hh"; then
	if test "x$TOPCOM_HOME" != "x/usr" -a "x$TOPCOM_HOME" != "x/usr/local"; then
		TOPCOM_CXXFLAGS="-I${TOPCOM_HOME}/include"
		TOPCOM_LIBS="-L${TOPCOM_HOME}/lib -lTOPCOM -lCHECKREG -lwrapgmp-gcc4"
	else
		TOPCOM_CXXFLAGS=
		TOPCOM_LIBS="-lTOPCOM -lCHECKREG -lwrapgmp-gcc4"
	fi	
	CXXFLAGS="${BACKUP_CXXFLAGS} ${TOPCOM_CXXFLAGS} ${GMP_CFLAGS}" 
	LIBS="${BACKUP_LIBS} ${TOPCOM_LIBS} ${GMP_LIBS}"

	AC_TRY_LINK([
#include "PointConfiguration.hh"
],
[ PointConfiguration pc;
],
[	topcom_found="yes"
	break
]
)
else
	topcom_found="no"
fi
done

if test "x$topcom_found" = "xyes" ; then		
	AC_SUBST(TOPCOM_CXXFLAGS)
	AC_SUBST(TOPCOM_LIBS)
	AC_DEFINE(HAVE_TOPCOM_LIB,1,[Define if the TOPCOM library is installed])
	HAVE_TOPCOM_LIB=yes
	AC_MSG_RESULT(found)
else
	AC_MSG_RESULT(not found)
	ifelse([$3], , :, [$3])
fi	

AM_CONDITIONAL(HAVE_TOPCOM_LIB, test "x$HAVE_TOPCOM_LIB" = "xyes")

CXXFLAGS=${BACKUP_CXXFLAGS}
LIBS=${BACKUP_LIBS}
#unset LD_LIBRARY_PATH

BACKUP_PATH=${PATH}
for TOPCOM_HOME in ${TOPCOM_HOME_PATH} 
do	
	PATH=${TOPCOM_HOME}/bin:$PATH
done
AC_PATH_PROG(TOPCOM_POINTS2TRIANGS, points2triangs)
AC_PATH_PROG(TOPCOM_POINTS2PLACINGTRIANG, points2placingtriang)
PATH=${BACKUP_PATH}

AC_DEFINE_UNQUOTED(TOPCOM_POINTS2TRIANGS, "${TOPCOM_POINTS2TRIANGS}", 
		   [The path to the TOPCOM program points2triangs.])
AC_DEFINE_UNQUOTED(TOPCOM_POINTS2PLACINGTRIANG, "${TOPCOM_POINTS2PLACINGTRIANG}", 
                   [The path to the TOPCOM program points2placingtriang.])
AM_CONDITIONAL(HAVE_TOPCOM_BIN, test x$TOPCOM_POINTS2TRIANGS != x -a x$TOPCOM_POINTS2PLACINGTRIANG != x )
if test x$TOPCOM_POINTS2TRIANGS != x -a x$TOPCOM_POINTS2PLACINGTRIANG != x; then
	AC_DEFINE(HAVE_TOPCOM_BIN,1,[Define if the TOPCOM binaries are installed])
fi

])
