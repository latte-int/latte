# Stolen from LinBox

AC_DEFUN([LB_MISC],
[

AC_ARG_WITH(default,
	    [ --with-default=<path>
					Add <path> to the default path for external package checking
					Set as default with /usr and /usr/local
	    ],
	    [if test "$withval" = yes ; then
			echo "Default path = /usr /usr/local"
			DEFAULT_CHECKING_PATH="/usr /usr/local"
	      else
			echo "Default path = $withval /usr /usr/local"
			DEFAULT_CHECKING_PATH="$withval /usr /usr/local"
	     fi
	     ],
	     [
		echo "Default path = /usr /usr/local"
		DEFAULT_CHECKING_PATH="/usr /usr/local"
             ])

])