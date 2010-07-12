#! /bin/sh
PARAMS="--triangulation=cdd --dualization=cdd --compute-vertex-cones=cdd"
echo "#################################"
echo "Checking count $PARAMS"
echo "#################################"
./test.pl "$PARAMS"
