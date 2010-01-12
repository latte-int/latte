#!/bin/sh
cat /dev/null > randomPolys.txt
cat /dev/null > forms.txt
errors=`maple -q integration/decomposeTest.mpl 2>&1| egrep "tests failed" | tr -d '"' | awk {'print $1'}`
echo "$errors tests failed."
if [ $errors -gt 0 ]; then
exit 1
fi;