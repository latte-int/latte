#!/bin/sh
cat /dev/null > integration/status.txt
echo "read(\"integration/benchmark.mpl\"):" > integration/test.mpl
echo "global filename:" >> integration/test.mpl
echo "filename:= \"$2\":" >> integration/test.mpl
echo "draw_table$1():" >> integration/test.mpl
maple -q integration/test.mpl
#rm -f integration/test.mpl
#maple -q integration/benchmark.mpl #| tee 2>&1 | tee -a integration/status.txt | egrep "total errors" | tr -d '"' | awk {'exit $1'}