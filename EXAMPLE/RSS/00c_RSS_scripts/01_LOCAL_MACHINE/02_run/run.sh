#!/bin/bash

export GROMP='/home/gsosso/CODES/GROMACS/gromacs-5.0.4/BIN_MPI/bin/grompp_mpi'
export MDR='/home/gsosso/CODES/GROMACS/gromacs-5.0.4/BIN_MPI/bin/mdrun_mpi'

# get current workind directory
pwd=$(pwd)

##########################################################################################
#                            SET VARIABLES                                               # 
##########################################################################################  
# make sure folders are set-up correctly
list_folders=($(ls -d rss_*))

file_topology="topol.top"
file_minimize="grompp.mdp"
file_index="index.ndx"
file_n_ice="N_ICE.dat"

##########################################################################################
#                            CHECKS                                                      # 
##########################################################################################  
echo "I found ${#list_folders[@]} folders to minimize"
echo " ... do you want to run the minimization on them? (y/n)"

read continue

if [ ${continue} == "n" ]
then
  echo "*** [ABORT] ***"
  exit
fi

##########################################################################################
#                            LOOP THROUGH ALL FOLDERS                                    # 
##########################################################################################  
for folder in ${list_folders[@]}
do
  if [ -d ${folder} ]
  then
    # change back to pwd (just to make sure everything is OK)
    cd $pwd 
    
    # change into folder to minimize
    cd $folder

    # get gro filename
    grofile="${folder}.gro"  # get *gro filename
    suffix=${folder}         # get suffix
    
    ##########################################################################################
    #                            PREPARE FILES                                               # 
    ##########################################################################################  
    # get N_SOL and use it for topology file
    n_sol=$(grep "OW" ${folder}.gro | wc -l)
    sed "s/N_SOL/${n_sol}/g" ${pwd}/${file_topology} > topol_TMP.top

    # append crystallite info to index file
    n_cryst=$(head -n 1 ${pwd}/${file_n_ice})

    # get n_atoms
    n_atoms=$(sed '2q;d' ${grofile})
    cp ${pwd}/${file_index} index.ndx

    idx_start=$((${n_atoms} - ${n_cryst}))
    idx_end=${n_atoms}

    idx_cryst=$(seq ${idx_start} ${idx_end})
    echo "[ ice ]" >> index.ndx
    echo ${idx_cryst} > tmp ; sed -r 's/(([^[:blank:]]+[[:blank:]]+){10})/\1\n/g' tmp >> index.ndx ; rm -r -f tmp

    ##########################################################################################
    #                            RUN MINIMIZATION                                            # 
    ##########################################################################################  
    # create tpr
    $GROMP -f ${pwd}/${file_minimize} -c ${grofile} -p topol_TMP.top -n index.ndx -o min_$suffix.tpr -maxwarn 1 > /dev/null 2>&1
    
    # run minimization
    mpirun -np 4 $MDR -ntomp 1 -maxh 1 -deffnm min_$suffix > /dev/null 2>&1

    ##########################################################################################
    #                            COLLECT ENERGY                                              # 
    ##########################################################################################  
    # get energy of relaxed structure
    energy=""
    energy=$(grep "Potential Energy  =" min_$suffix.log | awk '{print $4}') 

    echo "${folder}     ${energy}" >> ${pwd}/E_tot.dat
    
    # delete unneccesary files
    rm *tpr *trr *edr mdout* *log ${grofile} topol_TMP.top index.ndx

    # change back to pwd (safer than cd ..)
    cd $pwd 
  fi
done
