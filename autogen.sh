#!/bin/sh
#gnulib-tool --m4-base=m4/gnulib --update
aclocal -I m4 -I m4/gnulib
autoheader
autoconf
libtoolize --force || glibtoolize --force
automake --add-missing
