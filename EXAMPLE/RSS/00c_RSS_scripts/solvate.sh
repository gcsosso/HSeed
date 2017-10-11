#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Illegal number of parameters, please specify grofilename and z_cutoff"
    exit
fi

if [ ! -f $1 ]
then
  echo "$1 does not exist"
  exit
fi

outfile="solvated.gro"
z_cutoff=$2
##########################################################################################
#                            check if height OK                                          # 
##########################################################################################  
box_z=$(tail -n 1 $1 | awk '{print $3}')

if [ ${box_z} != "6.32592" ]
then
  echo "[ERROR]. box_z needs to be 6.32592 (${box_z})"
  echo "  --- ABORT ---"
  exit
fi


##########################################################################################
#                            SOLVATE                                                     # 
##########################################################################################  
sed -i 's/SOL/ICE/g' $1
genbox -cp $1 -cs /usr/local/gromacs/share/gromacs/top/tip4p.gro -o ${outfile} -vdwd 0.4

##########################################################################################
#                            REMOVE SOL THAT ARE TOO LOW                                 # 
##########################################################################################  
n_line=0
while read line
do
  n_line=$((n_line+1))
  if [[ ${line} == *"SOL"* ]]
  then
    if [[ ${line} == *"OW"* ]]
    then
      z=$(echo ${line} | awk '{ print $NF }')
      flag_del=$(echo "$z < ${z_cutoff}" | bc)
      if [ ${flag_del} == "1" ]
      then
        n_line_finish=$((n_line+3))
        sed -i "${n_line},${n_line_finish}d" ${outfile}
        n_line=$((n_line-4))
      fi
    fi
  fi
done < ${outfile}


##########################################################################################
#                            CHANGE GRO-FILE DETAILS                                     # 
##########################################################################################  
n_lines_total=$(wc -l solvated.gro | awk '{print $1}')
n_atoms=$((n_lines_total-3))

sed -i '2d' ${outfile}
sed -i "2i${n_atoms}" ${outfile}

# boxdimension
sed -i "${n_lines_total}s/6.32592/9.00000/" ${outfile}

##########################################################################################
#                            GENERATE INDEX FILE                                         # 
##########################################################################################  
# this already includes a category [ ICE ]
make_ndx -f ${outfile} -o index.ndx <<< q

##########################################################################################
#                            ICE --> SOL                                                 # 
##########################################################################################  
sed -i 's/ICE/SOL/g' ${outfile}
