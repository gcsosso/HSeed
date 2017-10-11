#!/bin/bash

# [ THIS SCRIPT IS SUPPOSED TO RUN ON MY LOCAL MACHINE  ]
#    ==> no parallelization of minimization task
#
dir_gro_files="./00a_start-structures"                    # contains start structures
dir_gromacs_files="./00b_gromacs_files/01_LOCAL_MACHINE/" # contains forcefield, MINIMIZE files, topology files
dir_RSS_scripts="./00c_RSS_scripts/01_LOCAL_MACHINE/"     # contains run-script for minimization
dir_minimize="02_minimize"                                # minimization folder
filename_n_ice="N_ICE.dat"

##########################################################################################
#                            CHECK ARGUMENTS                                             # 
##########################################################################################  
if [ "$#" -ne 1 ]
then
  echo "[ERROR]. The number of arguments is wrong."
  echo "         Expected: $0 N_ICE"
  exit
fi

N_ICE=$1

##########################################################################################
#                            check if all folders are here                               # 
##########################################################################################  
if [ ! -d ${dir_gromacs_files} ]
then
  echo "[ERROR]. The folder ${dir_gromacs_files} does not exist"
  exit
fi

if [ ! -d ${dir_RSS_script} ]
then
  echo "[ERROR]. The folder ${dir_RSS_scripts} does not exist"
  exit
fi

if [ ! -d ${dir_gro_files} ]
then
  echo "[ERROR]. The folder ${dir_gro_files} does not exist"
  exit
fi

##########################################################################################
#                            move gro-files into subfolder ./minimize                    # 
##########################################################################################  
mkdir ${dir_minimize}
cp -r ${dir_gro_files}/*gro ${dir_minimize}

cd ${dir_minimize}

##########################################################################################
#                            move each gro into folder                                   # 
##########################################################################################  
for file in $(ls *.gro)
do
  # cut out gro file ending
  suffix=${file%.gro}

  mkdir $suffix

  mv $file $suffix/
done

##########################################################################################
#                            copy gromacs-files and run-files                            # 
##########################################################################################  
# copy gromacs files
cp -r ../${dir_gromacs_files}/* .

# copy run files
cp ../${dir_RSS_scripts}/02_run/* .

echo "${N_ICE}" > ${filename_n_ice}

##########################################################################################
#                            change back                                                 # 
##########################################################################################  
cd ..
