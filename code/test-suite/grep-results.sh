#! /bin/sh
if test $# -ne 1 ; then
    echo "usage: grep-results.sh INSTANCE-REGEXP"
    exit 1
fi
INSTANCE_REGEXP=$1
cd /home/mkoeppe/w/latte/results
find . -name summary -print0 | xargs -0 grep $INSTANCE_REGEXP | grep -v Skipped | grep -v ERROR | sort 
