#!/bin/sh
if [ x/usr/local/bin/maple = x ]; then
  echo "Skipping integrating Hyper-rectangle test (needs Maple, which is not installed)"
  exit 0
fi
cat /dev/null > valuation/integrateHyperrectangleTest.status.txt
/usr/local/bin/maple -q valuation/integrateHyperrectangle.mpl #| tee 2>&1 | tee -a valuation/integrateHyperrectangleTest.status.txt | egrep "total errors" | tr -d '"' | awk {'exit $1'}
