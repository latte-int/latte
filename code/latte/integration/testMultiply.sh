#!/bin/sh
echo > integration/randomPolys.txt
echo > integration/results.txt
maple -q integration/multiplyTest.mpl
if [ "$?" -ne 0 ]; then echo "Multiplication failed"; exit 1; fi 