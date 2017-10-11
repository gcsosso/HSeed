#!/bin/bash

# [ THIS SCRIPT IS SUPPOSED TO RUN ON LEGION ]
#    ==> no parallelization of minimization task using strides of job-arrays
#
dir_gro_files="./00a_start-structures"                    # contains start structures
dir_gromacs_files="./00b_gromacs_files/02_LEGION"         # contains forcefield, MINIMIZE files, topology files
dir_RSS_scripts="./00c_RSS_scripts/02_LEGION"             # contains run-script for minimization
dir_minimize="02_minimize"                                # minimization folder
paramfile="run_params.txt"                                # contains folder names for array-stride
runfile="legion-gromacs_array.sh"                         # runfile for array-stride
filename_n_ice="N_ICE.dat"                                # file that stores the number of frozen molecules


##########################################################################################
#                            check if all arguments were specified                       # 
##########################################################################################  
if [ "$#" -ne 4 ]
then
  echo "[ERROR]. The number of arguments is wrong."
  echo "         Expected: $0 TITEL NJOBS STRIDE N_ICE"
  exit
fi

JOB_TITEL=$1
NJOBS=$2
STRIDE=$3
N_ICE=$4

pwd=$(pwd)

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

# generate paramfile that stored info for array-stride
if [ -f ${paramfile} ]
then
  rm ${paramfile}
fi

touch ${paramfile}

##########################################################################################
#                            move each gro into folder                                   # 
##########################################################################################  
for file in $(ls *.gro)
do
  # cut out gro file ending
  suffix=${file%.gro}

  mkdir $suffix

  mv $file $suffix/

  echo ${suffix} >> ${paramfile}
done

##########################################################################################
#                            copy gromacs-files and run-files                            # 
##########################################################################################  
# copy gromacs files
cp -r ../${dir_gromacs_files}/* .

# copy run files
cp ../${dir_RSS_scripts}/02_run/${runfile} .

# create N_ICE.dat file
echo "${N_ICE}" > ${filename_n_ice}

# replace space holders in runfile
sed -i "s/_TITEL_/${JOB_TITEL}/g" ${runfile}
sed -i "s/_NJOBS_/${NJOBS}/g" ${runfile}
sed -i "s/_STRIDE_/${STRIDE}/g" ${runfile}

##########################################################################################
#                            zip grofiles                                                # 
##########################################################################################  
cd ${pwd}/${dir_gro_files}

tar cfz start_structures.tar.gz *gro
rm *gro

##########################################################################################
#                            change back                                                 # 
##########################################################################################  
cd ${pwd}
