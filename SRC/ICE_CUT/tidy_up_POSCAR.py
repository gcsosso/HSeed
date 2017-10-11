#!/usr/bin/python

import sys
import numpy as np
from rotations import rotation_matrix_v2v
import transform_coordinates as tc # convert frac2cart / cart2frac
from pbc_computations import pbc_vector, pbc_distance

from file_formats import format_GRO, format_POSCAR

"""
the goal of this tool is to crop away O, H and X that are were separated by periodic boundary conditions.
It also helps to put the atom naming / ordering into the correct structure... (by explicitly finding atoms that form a SOL)

it is supposed to be used as pre-processor for make_hemisphere.py [ create ice crystallites/lattices for RSS ]
"""
##########################################################################################
#                            MAKROS                                                      # 
##########################################################################################  
IDX_START = 0
IDX_END   = 1

# check if all necessary arguments are here
if (len(sys.argv) != 5):
  print "[ERROR]. Please specify all arguments according to:"
  print "         %s INFILE.POSCAR OUTFILE.POSCAR Z_MIN Z_MAX" % sys.argv[0]
  sys.exit()

infilename_poscar   = sys.argv[1]
outfilename_poscar  = sys.argv[2]
z_min               = float(sys.argv[3])
z_max               = float(sys.argv[4])

##########################################################################################
#                            READ POSCAR FILE                                            # 
##########################################################################################  
file_poscar = format_POSCAR()

file_poscar.read_file(infilename_poscar)

##########################################################################################
#                            REMOVE ATOMS TOO HIGH/LOW                                   # 
##########################################################################################  
n_remove = 0
list_remove = []
for idx, atom in enumerate(file_poscar.list_atoms):
  if atom[file_poscar.Z_COORD] < z_min:
    n_remove += 1
    list_remove.append(idx)
  elif atom[file_poscar.Z_COORD] > z_max:
    n_remove += 1
    list_remove.append(idx)

# order list
list_remove_sorted = sorted(list_remove)

# remove elements
for idx in reversed(list_remove_sorted):
  del file_poscar.list_atoms[idx]

##########################################################################################
#                            FIND SOL MOLECULES THAT BELONG TOGETHER                     # 
##########################################################################################  
"""
 ase is a sh*t tool and makes a complete mess in the way it orders the molecules!
 we therefore have to explicitly find atoms that form a water molecule
"""

list_SOL = []
for atom_O in file_poscar.list_atoms:
  if atom_O[file_poscar.ATOMNAME] == "O":
    list_H_idx = []
    list_X_idx = []
    for idx, atom in enumerate(file_poscar.list_atoms):
      if atom[file_poscar.ATOMNAME] == "H":
        d = pbc_distance( [atom_O[file_poscar.X_COORD], atom_O[file_poscar.Y_COORD], atom_O[file_poscar.Z_COORD]],           \
                          [atom[file_poscar.X_COORD], atom[file_poscar.Y_COORD], atom[file_poscar.Z_COORD]],                 \
                          file_poscar.lattice_vector_1, \
                          file_poscar.lattice_vector_2, \
                          file_poscar.lattice_vector_3  )

        if d < 1.2:
          list_H_idx.append(idx)
          
      elif atom[file_poscar.ATOMNAME] == "X":
        d = pbc_distance( [atom_O[file_poscar.X_COORD], atom_O[file_poscar.Y_COORD], atom_O[file_poscar.Z_COORD]],           \
                          [atom[file_poscar.X_COORD], atom[file_poscar.Y_COORD], atom[file_poscar.Z_COORD]],                 \
                          file_poscar.lattice_vector_1, \
                          file_poscar.lattice_vector_2, \
                          file_poscar.lattice_vector_3  )

        if d < 0.5:
          list_X_idx.append(idx)

    if len(list_H_idx) != 2:
      print "[ERROR]. Code 0"
      print list_H_idx
      sys.exit()
    if len(list_X_idx) != 1:
      print "[ERROR]. Code 1"
      sys.exit()

    list_SOL.append([ atom_O, file_poscar.list_atoms[list_H_idx[0]], file_poscar.list_atoms[list_H_idx[1]], file_poscar.list_atoms[list_X_idx[0]] ])

##########################################################################################
#                            RESORT POSCAR INTO ORDER: O, H, X                           # 
########################################################################################## 
list_atom_names_unique = []

# expect O H X O H X O H X ...
if (len(file_poscar.list_atom_numbers) % 3) != 0 or file_poscar.list_atom_names[0] != "O":
    print "[ERROR]. Code 2"
    sys.exit()

for an in file_poscar.list_atom_names:
  if an not in list_atom_names_unique:
    list_atom_names_unique.append(an)

if list_atom_names_unique != ["O", "H", "X"]:
  print "[ERROR]. Code 3"
  sys.exit()

# resort
list_atoms_sorted = []
for atomtype in range(4):
  for mol in list_SOL:
    list_atoms_sorted.append(mol[atomtype])

##########################################################################################
#                            REASSIGN POSCAR                                             # 
##########################################################################################  
file_poscar.titel = "Z-PCB modified POSCAR"

file_poscar.list_atoms = list_atoms_sorted
file_poscar.write_file(outfilename_poscar)
