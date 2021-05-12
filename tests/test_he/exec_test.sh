#!/usr/bin/env bash
#
echo "Test hydrostatic equilibrium:"
execfile=FIRTEZ-dz.x
#
if [ ! -d out_dir ]
then
  mkdir out_dir
fi
#
echo "\tlink executable:"
ln -vfs ../../src/${execfile} .
echo "\trun firtez-dz:"
mpirun -n 2 ${execfile} control_file.dat
#

