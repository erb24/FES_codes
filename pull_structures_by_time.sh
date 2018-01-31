#!/bin/bash -l

#path=/home/promano/Datasets/Oligonucleotides/Trinucleotides
mol='1UBQ' #Oligonucleotide
#nbases=3 #Number of nucleic acid bases
for mode in 1 2
do
  echo "$mode" > mode
  nfrs=$(wc -l anly_${mode}.dat | awk '{ print $1 }') #Really twice the number of frames for a trinucleotide
  #natoms=95 #Total number of atoms in the oligonucleotide
  format="le4pd"
  #rm -rfv protname.txt

  #for i in ${mol} ${nbases} ${nfrs} ${natoms}
  #do
  #  echo ${i} >> protname.txt
  #done

  echo ${format} > anly
  ncoords=$(wc -l coords_${mode}.dat | awk '{ print $1 }')
  #ncoords2=$(wc -l coords2.dat | awk '{ print $1 }')
  echo "$ncoords"
  echo "$ncoords2"
  echo $ncoords > ncoords
  #echo $ncoords2 > ncoords2
  cp -v coords.dat coords
  #cp -v coords2.dat coords2
  f95 ensemble_finder.f95
  ./a.out
  #mv -v frames_*.ndx ../../Analysis/
  #mv -v frames2_*.ndx ../../Analysis/
  echo "1" > select.inp

  for state in `seq 1 ${ncoords}`
  do
      nframes=$(wc -l frames_${state}.ndx | awk '{ print $1 }')
      echo "$nframes" > nframes_${state}
      gmx trjconv -f ${mol}.xtc -s ${mol}_pro1.tpr -fr frames_${state}.ndx -o ~/Desktop/${mol}_${state}.pdb < select.inp
  
  done

  f95 average_structures.f95 #Find average structure
  ./a.out

done

exit
