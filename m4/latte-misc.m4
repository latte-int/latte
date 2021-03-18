# Stolen from LinBox

AC_DEFUN([LB_MISC],
[

AC_ARG_WITH(default,
	    [ --with-default=<path>
					Add <path> to the default path for external package checking
	    ],
	    [if test "$withval" = yes ; then
			DEFAULT_CHECKING_PATH="DEFAULTS"
	      else
			DEFAULT_CHECKING_PATH="$withval DEFAULTS"
	     fi
	     ],
	     [
		DEFAULT_CHECKING_PATH="DEFAULTS"
             ])

])
