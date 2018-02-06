#!/bin/bash

for i in `seq 1 10`
do

  rm -rfv interpol_${i}.pdb
  cp -v min${i}.pdb interpol_${i}.pdb

done

f95 ~/Desktop/codes/interpol_traj.f95

for i in `seq 1 10`
do

  echo "${i}" >mode
  ./a.out

done

for i in `seq 1 10`
do

  sed -n "5,1236p" max${i}.pdb >> interpol_${i}.pdb

done

ipython ~/Desktop/codes/overlay.py

exit


