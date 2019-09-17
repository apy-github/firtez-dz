#!/usr/bin/env bash
#
files=`ls -d ./test_* | grep -v 'inv'`
#
cwd="`pwd`"
#
for it in ${files}
do
  echo "Test: ", ${it}
  cd ${it}
  ./exec_test.sh > _out_test.txt
  echo "."
  cd "${cwd}"
done
#
# Inversions (they depend on previous synthesis)
#
files=`ls -d ./test_* | grep 'inv'`
#
for it in ${files}
do
  echo "Test: ", ${it}
  cd ${it}
  ./exec_test.sh > _out_test.txt
  echo "."
  cd "${cwd}"
done
#
