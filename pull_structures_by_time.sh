#!/bin/bash -l

codesdir=/home/ebeyerle/Desktop/codes
mol='1UBQ' #Molecule
nres=$(sed -n "2p" protname.txt)
natoms=$(sed -n "4p" protname.txt)
mult=$(( ${natoms} + 7 ))
for mode in `seq 1 5`
do
  echo "$mode" > mode
  nfrs=$(wc -l anly_${mode}.dat | awk '{ print $1 }') 
  format="le4pd"

  echo ${format} > anly
  ncoords=$(wc -l coords_${mode}.dat | awk '{ print $1 }')
  echo "$ncoords"
  echo "$ncoords2"
  echo $ncoords > ncoords
  cp -v coords_${mode}.dat coords
  #Selcting an average conformation from an ensemble:
  f95 ${codesdir}/test_ensemble_finder.f95
  ./a.out
  echo "1" > select.inp
  gmx trjconv -f ${mol}.xtc -s ${mol}_pro1.tpr -fr frames_${mode}.ndx -o all.pdb < select.inp
  cp -v all.pdb tmp
  mkdir -v mode_${mode}_states
  while read line 
  do 
    counter=$(( ${counter} +1 )) 
    echo "$line"
    tline=$( echo "${line}*${mult}" | bc -l)
    echo "$tline"
    sed -n "1,${tline}p" tmp >${mol}_${counter}.pdb
    cp -v ${mol}_${counter}.pdb mode_${mode}_states/
    sed -i "1,${tline}d" tmp
  done <"nframes_${mode}"
  rm -rfv tmp
  echo "$counter"

  for state in `seq 1 ${counter}`
  do
    echo "1 1" | gmx rmsf -f ${mol}_${state}.pdb -s ${mol}_pro1.tpr -ox avg_mode_${mode}_${state}.pdb -fit no
    rm -rfv ${mol}_${state}.pdb
  done
  
  #Concatenate PDB files
  rm -rfv avg_mode_${i}.pdb
  cat avg_mode_${mode}_*.pdb >> avg_mode_${mode}.pdb
  rm -rfv ${mol}_${mode}.pdb
done

rm -rfv ./#*#

exit
