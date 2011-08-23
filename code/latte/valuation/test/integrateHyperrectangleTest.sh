#!/bin/sh
if [ x/usr/local/bin/maple = x ]; then
  echo "Skipping integrating Hyper-rectangle test (needs Maple, which is not installed)"
  exit 0
fi
cat /dev/null > valuation/test/integrateHyperrectangleTest.status.txt
/usr/local/bin/maple -q valuation/test/integrateHyperrectangle.mpl 
