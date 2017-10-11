#!/bin/bash

#
# [ CALCULATE E_ADS OF MINIMIZES STRUCTURES ]
#   --> results stored in separate folder
# 
# two versions of E_ads are created: unsorted and sorted for easier post-processing
# the adsorption energy is calculated on the spot, using data provided in file_components
#
dir_minimize="02_minimize"
dir_results="_RESULTS_"
file_energy="E_tot.dat"
file_components="/home/gsosso/CODES/GIT/Hseed/EXAMPLE/COMPONENTS/E_Components.dat"
file_crystallite="00d_crystallite/DATA_CRYSTALLITE.dat"

if [ $# -ne 1 ]
then
  echo "[ERROR]. Number of arguments wrong!"
  echo "         $0 PREFIX_DIR"
  exit
fi

prefix_dir=$1

if [ ! -f ${file_components} ]
then
  echo "[ERROR]. ${file_components} does not exist."
  exit
fi

if [ ! -f ${file_crystallite} ]
then
  echo "[ERROR]. ${file_crystallite} does not exist."
  exit
fi

# read E_SOL (line 1), E_SF (line 2)
n_line=$(wc -l ${file_components} | awk '{print $1}')
if [ ${n_line} -ne 2 ]
then
  echo "[ERROR]. The number of lines in ${file_components} is wrong."
  echo "         I expect 2 lines (E_SOL, E_SF)"
  exit
fi
E_SOL=$(sed '1q;d' ${file_components})
E_SF=$(sed '2q;d' ${file_components})

# read n_sol_cryst (line 1), E_CRYST (line 2)
n_line=$(wc -l ${file_crystallite} | awk '{print $1}')
if [ ${n_line} -ne 2 ]
then
  echo "[ERROR]. The number of lines in ${file_crystallite} is wrong."
  echo "         I expect 2 lines (n_sol_cryst, E_CRYST)"
  exit
fi
n_sol_cryst=$(sed '1q;d' ${file_crystallite})
E_CRYST=$(sed '2q;d' ${file_crystallite})

pwd=$(pwd)

if [ ! -d ${dir_results} ]
then
  mkdir ${dir_results}
fi

cp ${dir_minimize}/${file_energy} ${dir_results}/01_${file_energy}

# check if file exists already
if [ -f ${dir_results}/02_E_ads_unsorted.dat ]
then
  rm ${dir_results}/02_E_ads_unsorted.dat 
fi

list_folders=($(ls -d ${dir_minimize}/${prefix_dir}_*))

echo "I found ${#list_folders[@]} folders to calculate E_ads in"
echo " ... do you want to continue? (y/n)"

read continue

if [ ${continue} != "y" ]
then
  echo "*** [ABORT] ***"
  exit
fi

for folder in ${list_folders[@]}
do
  grofile=($(ls ${folder}/min*.gro))

  # check if there is one and only one min*.gro in directory
  if [ ${#grofile[@]} -ne 1 ]
  then
    echo "[WARNING]. There are multiple min*.gro files in ${folder}."
    echo "           ${grofile[@]}"
    continue
  fi
 
  n_sol=$(grep "OW" ${grofile[0]} | wc -l)
  n_sol_CL=$(echo "print ${n_sol}-${n_sol_cryst}" | python)

  grofile_name=${grofile[0]##*/}
  grofile_name=${grofile_name:4:-4}
  E_TOTAL=$(grep ${grofile_name} ${dir_minimize}/${file_energy} | awk '{ print $2 }')

  if [ -z ${E_TOTAL} ]
  then
    echo "[WARNING]. Could not extract E_TOTAL for ${grofile_name}, string empty"
    continue
  fi

  # conversion kJ/mol -> eV: 1 kJ/mol == 0.0103642723 eV
  E_ADS=$(echo "print (${E_TOTAL}-${E_SF}-${n_sol}*${E_SOL})/${n_sol}*0.0103642723" | python )
  E_CL=$(echo "print (${E_TOTAL}-${E_SF}-${E_CRYST}-${n_sol_CL}*${E_SOL})/${n_sol_CL}*0.0103642723" | python )

  echo "${grofile_name} ${n_sol} ${n_sol_CL} ${E_TOTAL} ${E_ADS} ${E_CL}" >> ${dir_results}/02_E_ads_unsorted.dat

done

# sort energies
sort -g -k 5 ${dir_results}/02_E_ads_unsorted.dat > ${dir_results}/03_E_ads_sorted_TOTAL.dat
sort -g -k 6 ${dir_results}/02_E_ads_unsorted.dat > ${dir_results}/04_E_ads_sorted_CL.dat
