#!/bin/bash

# [ get distribution of x,y,z coordinate of OW ]

# get current workind directory
pwd=$(pwd)

dir_minimize="02_minimize"
dir_results="_RESULTS_"

# make sure folders are set-up correctly
list_folders=($(ls -d ${dir_minimize}/*_????))

echo "I found ${#list_folders[@]} folders to analyse"
echo " ... do you want to get the x/y/z distribution of them? (y/n)"

read continue

echo "name x y z" > ${dir_results}/xyz_distribution.dat

if [ ${continue} == "n" ]
then
  echo "*** [ABORT] ***"
  exit
fi

for folder in ${list_folders[@]}
do
  # get data for the first N OW
  #data=$(grep -m 1 "OW" $folder/*.gro | awk '{print $4, $5, $6}')
  
  # get all OW data
  grep "OW" $folder/*.gro | awk -v name="$folder" '{print name, $4, $5, $6}' >> ${dir_results}/xyz_distribution.dat
done
