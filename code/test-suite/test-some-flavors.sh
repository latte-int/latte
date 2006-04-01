#! /bin/sh
./test.pl "" > /dev/null &
./test.pl "--exp" > /dev/null &
./test.pl "--exp --maxdet=10" > /dev/null  &
./test.pl "--exp --maxdet=100" > /dev/null &
./test.pl "--exp --maxdet=1000" > /dev/null &
./test.pl "homog" > /dev/null &
