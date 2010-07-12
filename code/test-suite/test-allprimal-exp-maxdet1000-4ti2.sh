#! /bin/sh
PARAMS="--all-primal --exp --maxdet=1000 --triangulation=4ti2 --dualization=4ti2 --compute-vertex-cones=4ti2"
echo "#################################"
echo "Checking count $PARAMS"
echo "#################################"
./test.pl "$PARAMS"
