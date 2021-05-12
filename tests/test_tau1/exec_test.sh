#!/usr/bin/env bash
#
execfile=FIRTEZ-dz.x
ln -vfs ../../src/${execfile} .
echo "Test tau1:"
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

