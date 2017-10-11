#!/bin/bash

export nn=0

export GROMP='/home/gsosso/CODES/GROMACS/gromacs-5.0.4/BIN_MPI/bin/grompp_mpi'
export MDR='/home/gsosso/CODES/GROMACS/gromacs-5.0.4/BIN_MPI/bin/mdrun_mpi'

rm -r -f grompp.log job_craypath.sh mdout.mdp mdrun_mpi* report.dat  topol.tpr *.log *.ap2 *.rpt mdout.mdp gromacs.log grompp.log md.0.edr md.0.gro state.cpt topol.tpr traj.0.xtc dfs_surf.dat shoot.dat

$GROMP -f grompp.mdp -c conf.gro -p topol.top -o topol.tpr -n index.ndx -maxwarn 1 # > gromacs.log 2>&1
mpirun -np 16 $MDR -ntomp 1 -maxh 24 -s topol.tpr -o traj.${nn}.trr -x traj.${nn}.xtc -c md.${nn}.gro -e md.${nn}.edr -g md.${nn}.log # > gromacs.log 2>&1

exit 0

