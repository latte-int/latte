#! /bin/sh
./test.pl "" > /dev/null &
./test.pl "--exp" > /dev/null &
./test.pl "--exp --maxdet=10" > /dev/null  &
./test.pl "--exp --maxdet=100" > /dev/null &
./test.pl "--exp --maxdet=1000" > /dev/null &
./test.pl "homog" > /dev/null &
./test.pl "--maxdet=10" > /dev/null  &
./test.pl "--maxdet=100" > /dev/null &
./test.pl "--maxdet=1000" > /dev/null &
./test.pl "--irr" > /dev/null &
./test.pl "--irr --maxdet=10" > /dev/null &
./test.pl "--irr --maxdet=100" > /dev/null &
./test.pl "--irr --maxdet=1000" > /dev/null &
./test.pl "--exp --irr" > /dev/null &
./test.pl "--exp --irr --maxdet=10" > /dev/null &
./test.pl "--exp --irr --maxdet=100" > /dev/null &
./test.pl "--exp --irr --maxdet=1000" > /dev/null &
