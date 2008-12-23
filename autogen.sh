#!/bin/sh
#gnulib-tool --m4-base=m4/gnulib --update
aclocal -I m4
autoheader
autoconf
automake --add-missing
