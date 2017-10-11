#!/bin/bash -l

# Batch script to run a GROMACS job on Legion with the upgraded software
# stack under SGE with Intel MPI. Updated Oct 2015.

# 1. Force bash as the executing shell.
#$ -S /bin/bash

# 2. Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=5:00:0

# 3. Request 1 gigabyte of RAM per process.
#$ -l mem=1G

# 4. Request 15 gigabyte of TMPDIR space per node (default is 10 GB)
#$ -l tmpfs=15G

# 5. Set the name of the job.
#$ -N _TITEL_

# 6. Select the MPI parallel environment, here: serial job
#$ -pe mpi 1

# 6. Set up the job array. 
#$ -t 1-_NJOBS_:_STRIDE_

# 7. Set the working directory of the job to the current directory
#     containing your input files.
#    This *has* to be somewhere in your Scratch space, or else your
#     job will go into the Eqw state.
#$ -cwd

# 8. load gromacs
module load gromacs/5.0.4/intel-2015-update2

# Parse parameters
paramfile="run_params.txt"

if [ ! -f ${paramfile} ]
then
  echo "[ERROR]. The parameter-file ${paramfile} does not exit!"
  echo " *** ABORT ***"
  exit 
fi
 
# set some variables
pwd=$(pwd)
file_topology="topol_cm0_331.top"
file_minimize="grompp_MINIMIZE.mdp"
file_index="index_cm0_331.ndx"
file_n_ice="N_ICE.dat"


for (( i=${SGE_TASK_ID}; i<${SGE_TASK_ID}+_STRIDE_; i++ ))
do
  dir=`sed -n ${i}p $paramfile | awk '{print $1}'`

  ########################################################################################
  #      check if directory exists                                                       #
  ########################################################################################
  if [ ! -d ${dir} ]
  then
    echo "[ERROR]. ${dir} does not exist. Skipping it..."
  else
    cd ${dir}
    
    ##########################################################################################
    #      prepare topology and index files                                                  # 
    ##########################################################################################  
    # get N_SOL and use it for topology file
    n_sol=$(grep "OW" ${dir}.gro | wc -l)
    sed "s/N_SOL/${n_sol}/g" ${pwd}/${file_topology} > topol.top

    # append crystallite info to index file
    n_cryst=$(head -n 1 ${pwd}/${file_n_ice})

    # get n_atoms
    n_atoms=$(sed '2q;d' ${dir}.gro)
    cp ${pwd}/${file_index} index.ndx

    idx_start=$((${n_atoms} - ${n_cryst}))
    idx_end=${n_atoms}

    idx_cryst=$(seq ${idx_start} ${idx_end})
    echo "[ ice ]" >> index.ndx
    echo ${idx_cryst} >> index.ndx
    
    ######################################################################################
    #      prepare tpr file                                                              #
    ######################################################################################
    echo " --- [PREPARING] ${dir} ---"
    grompp -f ${pwd}/${file_minimize} -c ${dir}.gro -n index.ndx -p topol.top -o min_${dir}.tpr &> /dev/null
    
    ######################################################################################
    #      run minimization                                                              #
    ######################################################################################
    echo " --- [MINIMIZING] ${dir} ---"
    gerun mdrun_mpi -deffnm min_${dir} &> /dev/null
    
    ######################################################################################
    #      collect energy                                                                #
    ######################################################################################
    # get energy of relaxed structure
    energy=""
    energy=$(grep "Potential Energy  =" min_${dir}.log | awk '{print $4}') 

    echo "${dir}     ${energy}" >> ${pwd}/E_tot.dat
    
    # delete unneccesary files
    rm *tpr *trr *edr mdout* *log ${dir}.gro index.ndx topol.top *.pdb

    # change back to pwd (safer than cd ..)
    cd $pwd 
  fi
done

