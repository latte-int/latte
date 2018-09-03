#! /bin/sh
for a in `find . -name "testingSLandTopEhrhart" -prune -o -name "*.log" -print -o -name "log" -print`; do
    section=`echo $a | sed 's/[^[:alnum:]]/_/g;'`
    echo "###########################################################################"
    echo "travis_fold:start:${section}"
    echo "###LOG### $a ##############################"
    echo "###########################################################################"
    cat "$a"
    echo "travis_fold:end:${section}"
done
