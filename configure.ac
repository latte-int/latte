#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.59)
AC_INIT(latte, 1.2-mk-0.9+cvs)
AC_CONFIG_HEADERS([code/config.h:code/config.hin])
AM_INIT_AUTOMAKE(gnu)
AM_MAINTAINER_MODE
AC_CONFIG_FILES([Makefile EXAMPLES/Makefile code/Makefile code/latte/Makefile code/latte/normalize/Makefile code/test-suite/Makefile code/latte/adjacency/Makefile code/lidia-include.h:code/lidia-include.hin])
AC_CONFIG_FILES([code/test-suite/test.pl], [chmod +x code/test-suite/test.pl])
AC_PROG_RANLIB
AC_PROG_CC
AC_PROG_CXX

AC_LANG(C++)
AX_CXX_MAXOPT

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

# Compute installed location of auxiliary binaries.
computeadjacency_path=${libexecdir}/`echo ComputeAdjacency | sed "${program_transform_name}"`
AC_DEFINE_DIR([COMPUTEADJACENCY_PATH], [computeadjacency_path], [The path to the ComputeAdjacency binary.])

maple_script_dir=${datadir}/latte
AC_DEFINE_DIR([MAPLE_SCRIPT_DIR], [maple_script_dir], [The directory where LattE's Maple scripts are installed.])

LB_MISC
LB_CHECK_GMP(,,[AC_MSG_ERROR([GMP not found!
GMP version 3.1.1 or greater is required for LattE to compile. Please
make sure GMP is installed and specify its location with the option
--with-gmp=<prefix> when running configure.])
])
if test -z "GMP_HAVE_CXX" ; then 
  AC_MSG_ERROR([Unable to use the GMP C++ library.  
Please make sure that GMP was configured with "--enable-cxx".])
fi

AC_CACHE_SAVE
LB_CHECK_NTL(,,[AC_MSG_ERROR([NTL not found!
NTL version 5.4 or greater is required for LattE to compile. Please
make sure NTL is installed and specify its location with the option
--with-ntl=<prefix> when running configure.])
])
AC_CACHE_SAVE
LB_CHECK_LIDIA(,,
  AC_MSG_WARN([LiDIA not found.  We need LiDIA for computing Smith normal forms!]))
AC_CACHE_SAVE

# From CDD+
AC_PATH_PROGS(CDD, [cddr+_gmp cddr+])
AC_DEFINE_UNQUOTED(CDD_PATH, "${CDD}", [The path to the CDD+ binary.])
if test -z "$CDD"; then
  AC_MSG_ERROR([The `cddr+' binary from CDD+ is required for LattE.])			
fi
AC_CACHE_SAVE

# From cddlib
AC_PATH_PROGS(REDCHECK, [redcheck_gmp])
AC_DEFINE_UNQUOTED(REDCHECK_PATH, "${REDCHECK}", [The path to the cddlib redcheck binary.])
if test -z "$REDCHECK"; then
  AC_MSG_ERROR([The `redcheck' binary from cddlib is required for LattE.])			
fi

# Cddlib library
LB_CHECK_CDDLIB
if test "x$cddlib_found" = "xno" ; then
   # Older versions of cddlib are compiled with the C compiler;
   # also some people compile newer versions with the C compiler instead of the C++ compiler.
   AC_LANG(C)
   LB_CHECK_CDDLIB
   if test "x$cddlib_found" = "xyes" ; then
      AC_DEFINE(HAVE_CDDLIB_C_LINKAGE, 1, [Define if the CDDLIB library has C linkage])
   fi
   AC_LANG(C++)
fi   
if test "x$cddlib_found" = "xno" ; then
   AC_MSG_WARN([cddlib not found!
LattE needs it for building the CheckAdjacency program, which is
needed for the "vrep" mode.])
fi
AC_CACHE_SAVE
	  
# From LRS
AC_PATH_PROGS(LRS, [lrs1])
AC_DEFINE_UNQUOTED(LRS_PATH, "${LRS}", [The path to the lrs binary.])
if test -z "$LRS"; then
  AC_MSG_WARN([The `lrs1' binary from LRS is missing; LattE can optionally use LRS instead of CDD.])
fi
AC_CACHE_SAVE

# We need Maple to simplify some rational functions. 
AC_PATH_PROG(MAPLEBIN, maple)
AM_CONDITIONAL(HAVE_MAPLE, test x$MAPLEBIN != x )
AC_DEFINE_UNQUOTED(MAPLE_PATH, "${MAPLEBIN}", [The path to the Maple program.])

# TOPCOM for triangulations (optional)
LB_CHECK_TOPCOM

# 4ti2 for triangulations (optional)
LB_CHECK_FORTYTWO

# CPLEX for experimental code (optional)
AC_ARG_WITH(cplex-sharedlib,
    AS_HELP_STRING([--with-cplex-sharedlib=/PATH/TO/FILE.so],
                   [Specify name of CPLEX shared library]),
    [  CPLEX_SHAREDLIB=$withval ])

AC_ARG_WITH(cplex-include, 
    AS_HELP_STRING([--with-cplex-include=DIR], [Specify location of CPLEX include files]),
    [  CPLEX_CFLAGS="-I$withval" ])
AC_SUBST(CPLEX_CFLAGS)

AM_CONDITIONAL(HAVE_CPLEX, test x$CPLEX_SHAREDLIB != x -a x$CPLEX_CFLAGS != x)
if test x$CPLEX_SHAREDLIB != x -a x$CPLEX_CFLAGS != x ; then
   AC_DEFINE(HAVE_CPLEX, 1, [Define if the CPLEX callable library is available])
  CPLEX_LDFLAGS="-L`dirname ${CPLEX_SHAREDLIB}` -R`dirname ${CPLEX_SHAREDLIB}` -l`basename ${CPLEX_SHAREDLIB} .so | sed 's/^lib//' ` -lm -lnsl -lsocket -lpthread "
fi
AC_SUBST(CPLEX_LDFLAGS)

# Experimental, not release-ready code.
AC_ARG_ENABLE(experiments, AS_HELP_STRING([--enable-experiments],
	[use experimental code]))
AM_CONDITIONAL(HAVE_EXPERIMENTS, test x$enable_experiments = xyes)
if test x$enable_experiments = xyes ; then
   AC_DEFINE(HAVE_EXPERIMENTS, 1, [Define if experimental code is enabled])
fi

AC_OUTPUT
