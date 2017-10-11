#!/bin/bash

# [ EXTRACT THE BEST STRUCTURES FROM E_ads_sorted.dat ]
# 
#   this script needs to be executed in folder containing ./minimize
#   the extracted structures are stored in results directory
#
#   n_padding_zeros ... how many digits for the counter 
#     ( 1 --> 0,     1,   2, ... )
#     ( 3 --> 000, 001, 002, ... )

dir_results="_RESULTS_"
dir_structures="extracted_structures"

# specify which structures you want to extract 
structures=(out_92 out_46)

if [ ! -d ${dir_results} ]
then
  echo "[ERROR]. The directory _RESULTS_ does not exist"
  echo "         Make sure to run ./01_get_energies.sh first"
  exit
fi

if [ $# -ne 1  ]
then
  echo "[ERROR]. Please specify how many structures you want to extract"
  echo "         $0 n_padding_zeros"
  exit
fi

pwd=$(pwd)

cd ${dir_results}
mkdir ${dir_structures}

# copy specified structures
counter=0
for grofile in ${structures[@]}
do
  # pad $counter with leading zeros
  printf -v c "%0$2d" ${counter}

  cp ${pwd}/minimize/${grofile}/minimize*gro ${dir_structures}
  
  # get filename of file just copied
  newest=$(ls -rtlh ${dir_structures} | tail -n 1 | awk '{print $9}' )

  mv "${dir_structures}/${newest}" "${dir_structures}/${c}_${newest}"
  
  # increase counter for next time
  let "counter+=1"
done
