#!/bin/bash

# [ get orientation of seed on top of surface ]

# get current workind directory
pwd=$(pwd)

dir_minimize="02_minimize"
dir_results="_RESULTS_"
file_results="05_seed_orientation.dat"
src_measure_orientation="/home/pp/programs/python/measure_seed_orientation/measure_seed_orientation.py"

if [ $# -ne 3 ]
then
  echo "[ERROR]. Number of arguments ($#) wrong!"
  echo "         $0 polymorph face n_sol"
  echo "         polymorph ... {H, C}"
  echo "         face      ... H: {001, 100, 110}, C: {100, 111}"
  echo "         n_sol     ... target number of H2O"
  exit
fi

polymorph=$1
face=$2
n_sol=$3

# make sure folders are set-up correctly
list_folders=($(ls -d ${dir_minimize}/*_????))

echo "I found ${#list_folders[@]} folders to analyse"
echo " ... do you want to measure seed orientation of them? (y/n)"

read continue

if [ ${continue} == "n" ]
then
  echo "*** [ABORT] ***"
  exit
fi

if [ -f ${dir_results}/${file_results} ]
then
  rm ${dir_results}/${file_results}
fi

for folder in ${list_folders[@]}
do
  grofile=($(ls ${folder}/min*gro))
  if [ ${#grofile[@]} -ne 1 ]
  then
    echo "[WARNING]. More than one grofile in ${folder}: ${grofile[@]}"
    echo "${folder} NA NA" >> ${dir_results}/${file_results}
  fi

  result=$(python ${src_measure_orientation} ${grofile[0]} ${polymorph} ${face} ${n_sol})

  if [ $(echo ${result} | wc -w) != 2 ]
  then
    echo ""
    echo ${result}
    echo " --- ABORTING RUN ---"
    echo ""
    exit
  fi

  echo ${folder##*/} ${result} >> ${dir_results}/${file_results}
done
