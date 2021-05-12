#!/usr/bin/env bash
#
execfile=FIRTEZ-dz.x
echo "Test Syn1:"
#
if [ ! -d out_dir ]
then
  mkdir out_dir
fi
#
echo "\tlink executable:"
ln -vfs ../../src/${execfile} .
echo "\trun firtez-dz:"
mpirun -n 2 ./${execfile} control_file.dat
#

