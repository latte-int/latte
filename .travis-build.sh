#! /bin/bash
autoreconf -fi || exit 1
./configure
if [ $? != 0 ]; then
    echo "==================== config.log ===================="
    cat config.log
    exit 1
fi
make -j2 distcheck || exit 1
