#! /bin/sh
PARAMS="--all-primal --exp --maxdet=1000"
echo "#################################"
echo "Checking count $PARAMS"
echo "#################################"
./test.pl "$PARAMS"
