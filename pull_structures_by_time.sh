#!/bin/bash -l

codesdir=/home/ebeyerle/Desktop/codes
#path=/home/promano/Datasets/Oligonucleotides/Trinucleotides
mol='1UBQ' #Molecule
#nbases=3 #Number of nucleic acid bases
for mode in 1 2 3
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
  cp -v coords_${mode}.dat coords
#Selcting an average conformation from an ensemble:
  f95 ensemble_finder.f95
  ./a.out
  #mv -v frames_*.ndx ../../Analysis/
  #mv -v frames2_*.ndx ../../Analysis/
  #echo "1" > select.inp
  for state in `seq 1 ${ncoords}`
  do
      nframes=$(wc -l frames_${state}.ndx | awk '{ print $1 }')
      echo "$nframes" > nframes_${state}
      gmx trjconv -f ${mol}.xtc -s ${mol}_pro1.tpr -fr frames_${state}.ndx -o ${mol}_${state}.pdb < select.inp
  
  #f95 ~/Desktop/codes/average_structures.f95 #Find average structure
  #./a.out

  #Or use GROMACS:
  echo "1 1" | gmx covar -f ${mol}_${state}.pdb -s ${mol}_pro1.tpr -av avg_mode_${mode}_${state}.pdb -fit no
  rm -rfv ./#*#
  done

#Single structure from click:
  f95 ${codesdir}/structure_finder.f95
  ./a.out
  echo "1" > select.inp
  for state in `seq 1 ${ncoords}`
  do
      gmx trjconv -f ${mol}.xtc -s ${mol}_pro1.tpr -fr struct_${state}.ndx -o ${mol}_mode_${mode}_struct_${state}.pdb < select.inp
  done
  
#Concatenate PDB files
  rm -rfv avg_mode_${i}.pdb
  cat avg_mode_${mode}_*.pdb >> avg_mode_${i}.pdb
  rm -rfv ${mol}_struct${mode}.pdb
  cat ${mol}_mode_${mode}_struct_*.pdb >> ${mol}_struct${mode}.pdb
done

exit
