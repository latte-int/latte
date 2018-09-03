#! /bin/sh
# First test-suite summary logs
LOGS=`find . -name "testingSLandTopEhrhart" -prune -o -name "config.log" -prune -o -name "*.log" -print`
# Then detailed results logs
LOGS="$LOGS `find . -name "log" -print`"
# Finally the long config.log
LOGS="$LOGS `find . -name config.log`"
for a in $LOGS; do
    section=`echo $a | sed 's/[^[:alnum:]]/_/g;'`
    echo "###########################################################################"
    echo "travis_fold:start:${section}"
    echo "###LOG### $a ##############################"
    echo "###########################################################################"
    cat "$a"
    echo "travis_fold:end:${section}"
done
