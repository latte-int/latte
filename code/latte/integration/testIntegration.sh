#!/bin/sh
cat /dev/null > integration/status.txt

maple -q integration/decomposeTest.mpl #| tee 2>&1 | tee -a integration/status.txt | egrep "total errors" | tr -d '"' | awk {'exit $1'}