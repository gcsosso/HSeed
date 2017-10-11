#!/bin/bash

export GROMP='/home/gsosso/CODES/GROMACS/gromacs-5.0.4/BIN_MPI/bin/grompp_mpi'
export MDR='/home/gsosso/CODES/GROMACS/gromacs-5.0.4/BIN_MPI/bin/mdrun_mpi'


# create tpr
$GROMP -n index.ndx -o topol.tpr -maxwarn 1 > grompp.log 2>&1

# run minimization
mpirun -np 4 $MDR -ntomp 1 -maxh 1  > gromacs.log 2>&1

