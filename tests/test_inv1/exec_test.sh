#!/usr/bin/env bash
#
echo "Test Inv1:"
#
if [ ! -d out_dir ]
then
  mkdir out_dir
fi
#
echo "\tlink executable:"
ln -vfs ../../src/firtez-dz.x .
echo "\trun firtez-dz:"
mpirun -n 2 ./firtez-dz.x control_file.dat
#

