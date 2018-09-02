#! /bin/bash
autoreconf -fi || exit 1
./configure || exit 1
make -j2 distcheck
