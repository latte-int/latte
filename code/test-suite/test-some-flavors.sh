#! /bin/sh
on omega 'nohup ./test.pl "--exp --irr" > /dev/null &' >/dev/null &
#./test.pl "--exp --irr --maxdet=10" > /dev/null &
#./test.pl "--exp --irr --maxdet=100" > /dev/null &
on omega 'nohup ./test.pl "--exp --irr --maxdet=500" > /dev/null &'  >/dev/null &
#./test.pl "--exp --irr --maxdet=1000" > /dev/null &
#./test.pl "--exp --irr --maxdet=10000" > /dev/null &
#./test.pl "--irr" > /dev/null &
#./test.pl "--irr --maxdet=10" > /dev/null &
#./test.pl "--irr --maxdet=100" > /dev/null &
#./test.pl "--irr --maxdet=1000" > /dev/null &
#./test.pl "--irr --maxdet=10000" > /dev/null &
on epsilon 'nohup ./test.pl "--all-primal --exp" > /dev/null &' >/dev/null &
#./test.pl "--all-primal --maxdet=10" > /dev/null &
#./test.pl "--all-primal --maxdet=100" > /dev/null &
on epsilon 'nohup ./test.pl "--all-primal --maxdet=500 --exp" > /dev/null &' >/dev/null &
#./test.pl "--all-primal --maxdet=1000" > /dev/null &
#./test.pl "--all-primal --maxdet=10000" > /dev/null &
on epsilon 'nohup ./test.pl "" > /dev/null &' >/dev/null &
#./test.pl "--exp" > /dev/null &
#./test.pl "--exp --maxdet=10" > /dev/null  &
#./test.pl "--exp --maxdet=100" > /dev/null &
#./test.pl "--exp --maxdet=1000" > /dev/null &
#./test.pl "--exp --maxdet=10000" > /dev/null &
on epsilon 'nohup ./test.pl "homog" > /dev/null &' >/dev/null &
#./test.pl "--maxdet=10" > /dev/null  &
#./test.pl "--maxdet=100" > /dev/null &
#./test.pl "--maxdet=1000" > /dev/null &
#./test.pl "--maxdet=10000" > /dev/null &
