#!/bin/bash -l

codesdir=~/codes
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

  mkdir -v mode_${mode}_states
  counter=0
  echo "1" > select.inp
  cp -v frames_${mode}.ndx tmp.ndx
  while read line 
  do 
    counter=$(( ${counter} +1 )) 
    echo "$line"
    tline=$( echo "${line}*${mult}" | bc -l)
    echo "$tline"
    #sed -n "1,${tline}p" tmp >${mol}_${counter}.pdb
    #cp -v ${mol}_${counter}.pdb mode_${mode}_states/
    #sed -i "1,${tline}d" tmp
    sed -n "1,$(( 1+ ${line} ))p" tmp.ndx > frames_${mode}_${counter}.ndx
    cp -v frames_${mode}_${counter}.ndx mode_${mode}_states/
    gmx trjconv -f ../../${mol}.xtc -s ../../${mol}_2.tpr -fr frames_${mode}_${counter}.ndx -o ${mol}_${counter}.pdb < select.inp
    mv -v 1UBQ_${counter}.pdb mode_${mode}_states/
    sed -i "2,$(( 1+ ${line} ))d" tmp.ndx
  done <"nframes_${mode}"
  rm -rfv tmp
  rm -rfv tmp.ndx
  echo "$counter"
  cp -v all.pdb tmp

  for state in `seq 1 ${counter}`
  do
    echo "1 1" | gmx rmsf -f mode_${mode}_states/${mol}_${state}.pdb -s ../../${mol}_2.tpr -ox mode_${mode}_states/avg_mode_${mode}_${state}.pdb -fit no
    rm -rfv ${mol}_${state}.pdb
  done
  
  #Concatenate PDB files
  rm -rfv avg_mode_${i}.pdb
  cat avg_mode_${mode}_*.pdb >> avg_mode_${mode}.pdb
  rm -rfv ${mol}_${mode}.pdb
done

rm -rfv ./#*#

exit
