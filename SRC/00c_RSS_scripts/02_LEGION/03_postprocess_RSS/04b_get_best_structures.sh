#!/bin/bash

# [ EXTRACT THE BEST STRUCTURES FROM E_ads_sorted.dat ]

dir_results="_RESULTS_"
dir_best_structures="best_structures"
dir_minimize="02_minimize"
file_energies_sorted="03_E_ads_sorted.dat"
prefix="min_SOL_"  # prefix added to grofile (to distinguish different relaxations)
                   # this is important for eg. only_HW relax followed by full_relax

if [ ! -d ${dir_results} ]
then
  echo "[ERROR]. The directory _RESULTS_ does not exist"
  echo "         Make sure to run ./03_get_energies.sh first"
  exit
fi

if [ $# -ne 2  ]
then
  echo "[ERROR]. Please specify how many structures you want to extract"
  echo "         $0 n_struct n_padding_zeros"
  exit
fi

pwd=$(pwd)

cd ${dir_results}

mkdir ${dir_best_structures}

best_structures=($(head -n $1 ${file_energies_sorted} | awk '{print $1}'))

# copy best structures
counter=0
for grofile in ${best_structures[@]}
do
  # pad $counter with leading zeros
  printf -v c "%0$2d" ${counter}

  cp ${pwd}/${dir_minimize}/${grofile}/${prefix}*gro ${dir_best_structures}

  # give it a little processing time, sometimes the next command (ls -rtlh) can be messed up
  sleep 0.1

  # get filename of file just copied
  newest=$(ls -rtlh ${dir_best_structures} | tail -n 1 | awk '{print $9}' )

  mv "${dir_best_structures}/${newest}" "${dir_best_structures}/${c}_${newest}"
  
  # increase counter for next time
  let "counter+=1"
done
