#!/bin/bash

##################################################
#   SET UP VASP MINIMIZATIONS                    #
#    (1) convert gro -> POSCAR                   #
#    (2) copy vasp files                         #
#    (3) pbs title                               #
##################################################

pwd=$(pwd)

# folder containing all the files needed by vasp
vaspfiles="00_vasp-files/01_NO_DIPOLE_CORR"
scriptfile="script-salviati.sh"

# array containing the number of directories for minimization
#list_dir=(Al-Al_2x Al-Al_Si-Si Al-Si_2x Al-Si_Al-Al Al-Si_Si-Si Si-Si_2x)
list_dir=(test)

# parameters for gro2poscar
filename_out="POSCAR"
filename_ref="slab_reference.CONTCAR" # file has to be in directory of gro2poscar.py
multiply_x=3
multiply_y=2
multiply_z=1
n_h2o=2

echo "I will prepare the VASP run in ${#list_dir[@]} directories..."
echo "  ... is this OK? (y/n)"

read continue

if [ ${continue} == "n" ]
then
  echo "*** [ABORT] ***"
  exit
fi

# check if the folder vasp-files exists
if [ ! -d ${vaspfiles} ]
then
  echo "[ERROR]. The folder ${vaspfiles} does not exist"
  echo "         Please create the folder in ${DFT_dir}. It should contain:"
  echo "           (+) slab_reference.CONTCAR"
  echo "           (+) INCAR"
  echo "           (+) KPOINTS"
  echo "           (+) POTCAR"
  echo "           (+) script-salviati.sh"
  exit
fi

# execute gro2poscar.py in each folder
for dir in ${list_dir[@]}
do
  cd $dir

  # make 01_NO_DIPOLE_CORR
  mkdir 01_NO_DIPOLE_CORR
  mv * 01_NO_DIPOLE_CORR/ &> /dev/null
  cd 01_NO_DIPOLE_CORR

  list_gro=($(ls *gro))

  if [ ${#list_gro[@]} -ne 1 ]
  then
    echo "I found more than one gro-file in ${dir}"
    echo "${list_gro[@]}"
    echo "Please remove all but one gro files"
    echo "*** [ABORT] ***"
    exit
  fi
  
  filename_in=${list_gro[0]} 
  
  if [ ! -f ${filename_in} ]
  then 
    echo "[ERROR]. gro-file ${filename_in} doesn't exist"
    exit
  fi

  # copy vasp-files
  cp ${pwd}/${vaspfiles}/* .

  # execute gro2poscar.py
  gro2poscar.py ${filename_in} ${filename_ref} ${multiply_x} ${multiply_y} ${multiply_z} ${n_h2o} ${filename_out}

  # set jobname in archer-vasp.sh
  sed -i -e "s/TITEL/${dir}/g" ${scriptfile}

  cd ${pwd}
done

