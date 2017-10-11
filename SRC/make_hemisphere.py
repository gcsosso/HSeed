#!/usr/bin/python

"""
[ crop hemisphere out of POSCAR file, useful for creating heteronuc crystallites ]

#------------------------#
# prepare surface POSCAR #
#------------------------#
the POSCAR file is ideally generated using convert_gro-2-poscar and then using the ase to generate
a slab with appropriate orientation, see the following example script

   from ase import Atoms
   from ase.lattice.surface import surface
   from ase.io import *
   ice = read("file.POSCAR")
   sf = surface(ice, (1,1,0), 1)
   sf.center(vacuum=0, axis=2)
   write("SF.POSCAR", sf)

done in: /home/pp/work/Ice/Ice-Crystallites/02_heterogeneous/TIP4P_ICE/02_hemispheres

#------------------------#
# tidy up surface POSCAR #
#------------------------#
once, the POSCAR file has been created with ase, it might be necessary to tidy it up.
The order of atoms for example might be a mess (O H X O H X O H X ...)

do this with:
    /home/pp/programs/python/tidy_up_POSCAR/tidy_up_POSCAR.py

#-------------------#
# create hemisphere #
#-------------------#
then, specifying the number of SOL molecules that are wanted at the end, figure out the radius
of the hemispherical crystallite

also, produce some info about the crystallite:
  R, N_SOL, N_LATTICE, N_CRYSTALLITE

#---------------------------------------------------------#
# works to generate lattices for mW / TIP4P equally well! #
#---------------------------------------------------------#
"""

import sys
import numpy as np
from rotations import rotation_matrix_v2v
import transform_coordinates as tc # convert frac2cart / cart2frac
from pbc_computations import pbc_vector, pbc_distance

from file_formats import format_GRO, format_POSCAR

##########################################################################################
#                            MAKROS                                                      # 
##########################################################################################  
CENTER = 0
RADIUS = 1


# check if all necessary arguments are here
if (len(sys.argv) != 6):
  print "[ERROR]. Please specify all arguments according to:"
  print "         %s INFILE.POSCAR OUTFILE.gro OUTFILE.lattice nSOL_target FLAG_Z_LAYER_REMOVE" % sys.argv[0]
  sys.exit()

infilename_poscar   = sys.argv[1]
outfilename_gro     = sys.argv[2]
outfilename_lattice = sys.argv[3]
nSOL_target         = int(sys.argv[4])
flag_z_layer_remove = sys.argv[5] # True / False

"""
 all OW atoms below Z_LATTICE_CUTOFF will not be included in the gro-file,
 but will instead be used to create a lattice file that matches the crystallite.

 This flexibility in creating lattices is needed for cases were we might need 2 layer
 thick lattices.

 The value of the cutoff is obtained from directly inspecting the ice-surface used for the crop

 keep in mind that part of the cropped hemisphere will be used as lattice-file (up to Z_LATTICE_CUTOFF),
 whereas the rest will be used as crystallite (gro-file)
  ==> this way, the underlying lattice matches the structure of hemispherical crystallite
"""

##########################################################################################
#                            Ih hemisphere settings (TIP4P/Ice)                          # 
##########################################################################################  
#======================#
# 1L lattice, rest gro #
#======================#
#----------#
# Ih basal #
#----------#
#Z_LAYER_REMOVE   = 1.3 # all OW (and corresponding HW, MW) below this z will be removed before crop
#Z_LATTICE_CUTOFF = 5.8 # after shifting the structure down to remove uncomplete bottom layer:
#                       # this value is shifted down as well
#                       # ==> this means, its value can be read directly from inputfile.POSCAR

#------------------#
# Ih primary prism #
#------------------#
#Z_LAYER_REMOVE   = -1.0 # all OW (and corresponding HW, MW) below this z will be removed before crop
#Z_LATTICE_CUTOFF = 3.0  # after shifting the structure down to remove uncomplete bottom layer:
#                        # this value is shifted down as well
#                        # ==> this means, its value can be read directly from inputfile.POSCAR

#--------------------#
# Ih secondary prism #
#--------------------#
#Z_LAYER_REMOVE   = -1.0 # all OW (and corresponding HW, MW) below this z will be removed before crop
#Z_LATTICE_CUTOFF = 13.0 # after shifting the structure down to remove uncomplete bottom layer:
#                        # this value is shifted down as well
#                        # ==> this means, its value can be read directly from inputfile.POSCAR


#======================#
# 2L lattice, rest gro #
#======================#
#----------#
# Ih basal #
#----------#
#Z_LAYER_REMOVE   = 1.3 # all OW (and corresponding HW, MW) below this z will be removed before crop
#Z_LATTICE_CUTOFF = 9.5 # after shifting the structure down to remove uncomplete bottom layer:
#                       # this value is shifted down as well
#                       # ==> this means, its value can be read directly from inputfile.POSCAR

#------------------#
# Ih primary prism #
#------------------#
#Z_LAYER_REMOVE   = -1.0 # all OW (and corresponding HW, MW) below this z will be removed before crop
#Z_LATTICE_CUTOFF = 7.0  # after shifting the structure down to remove uncomplete bottom layer:
#                        # this value is shifted down as well
#                        # ==> this means, its value can be read directly from inputfile.POSCAR

#--------------------#
# Ih secondary prism #
#--------------------#
#Z_LAYER_REMOVE   = -1.0 # all OW (and corresponding HW, MW) below this z will be removed before crop
#Z_LATTICE_CUTOFF = 15.5 # after shifting the structure down to remove uncomplete bottom layer:
#                        # this value is shifted down as well
#                        # ==> this means, its value can be read directly from inputfile.POSCAR

##########################################################################################
#                            Ic hemisphere settings (TIP4P)                              # 
##########################################################################################  
#======================#
# 1L lattice, rest gro #
#======================#
#----------------#
# Ic 100/010/001 #
#----------------#
Z_LAYER_REMOVE   = -1.0 # all OW (and corresponding HW, MW) below this z will be removed before crop
#Z_LATTICE_CUTOFF = 12.5 # after shifting the structure down to remove uncomplete bottom layer:
#                        # this value is shifted down as well
#
# 2L lattice
Z_LATTICE_CUTOFF = 14.0 # ==> this means, its value can be read directly from inputfile.POSCAR
## 3L lattice
#Z_LATTICE_CUTOFF = 15.0 # after shifting the structure down to remove uncomplete bottom layer:
#                        # this value is shifted down as well
#                        # ==> this means, its value can be read directly from inputfile.POSCAR

#--------#
# Ic 111 #
#--------#
#Z_LAYER_REMOVE   = -1.0 # all OW (and corresponding HW, MW) below this z will be removed before crop
#Z_LATTICE_CUTOFF = 15.0 # after shifting the structure down to remove uncomplete bottom layer:
#                        # this value is shifted down as well
#                        # ==> this means, its value can be read directly from inputfile.POSCAR

#======================#
# 2L lattice, rest gro #
#======================#
#----------------#
# Ic 100/010/001 #
#----------------#
#Z_LAYER_REMOVE   = -1.0 # all OW (and corresponding HW, MW) below this z will be removed before crop
#Z_LATTICE_CUTOFF = 13.8 # after shifting the structure down to remove uncomplete bottom layer:
#                        # this value is shifted down as well
#                        # ==> this means, its value can be read directly from inputfile.POSCAR

#--------#
# Ic 111 #
#--------#
#Z_LAYER_REMOVE   = -1.0 # all OW (and corresponding HW, MW) below this z will be removed before crop
#Z_LATTICE_CUTOFF = 18.5 # after shifting the structure down to remove uncomplete bottom layer:
#                        # this value is shifted down as well
#                        # ==> this means, its value can be read directly from inputfile.POSCAR


##########################################################################################
#                            Ih hemisphere settings (mW)                                 # 
##########################################################################################  
#======================#
# FULL LATTICE, NO gro #
#======================#
#----------#
# Ih basal #
#----------#
#Z_LAYER_REMOVE   = 1.3  # all OW below this z will be removed before crop
#Z_LATTICE_CUTOFF = 50.0 # use all points within hemisphere as latticepoints 

#------------------#
# Ih primary prism #
#------------------#
#Z_LAYER_REMOVE   = -1.0 # all OW below this z will be removed before crop
#Z_LATTICE_CUTOFF = 50.0 # use all points within hemisphere as latticepoints 

#--------------------#
# Ih secondary prism #
#--------------------#
#Z_LAYER_REMOVE   = -1.0 # all OW below this z will be removed before crop
#Z_LATTICE_CUTOFF = 50.0 # use all points within hemisphere as latticepoints 


##########################################################################################
#                            Ic hemisphere settings (mW)                                 # 
##########################################################################################  
#======================#
# FULL LATTICE, NO gro #
#======================#
#----------------#
# Ic 100/010/001 #
#----------------#
#Z_LAYER_REMOVE   = -1.0 # all OW below this z will be removed before crop
#Z_LATTICE_CUTOFF = 50.0 # use all points within hemisphere as latticepoints 

#--------#
# Ic 111 #
#--------#
#Z_LAYER_REMOVE   = -1.0 # all OW below this z will be removed before crop
#Z_LATTICE_CUTOFF = 50.0 # use all points within hemisphere as latticepoints 


##########################################################################################
#                            READ POSCAR FILE                                            # 
##########################################################################################  
file_poscar = format_POSCAR()

file_poscar.read_file(infilename_poscar)

# consistency check, does #OW correspond to #MW and #HW/2

n_O = file_poscar.list_atom_numbers[0]
n_H = file_poscar.list_atom_numbers[1]
n_X = file_poscar.list_atom_numbers[2]

if n_O != n_X and n_O*2 != n_X and n_O*4 != len(file_poscar.list_atoms):
    print "[ERROR]. The water strucutre in", infilename_poscar, "does not make sense."
    print "         n_O, n_H and n_X do not have the correct ratio (assuming O,H,X order)"
    print n_O, n_H, n_X
    sys.exit()

# consistency check #2: is idx_O, idx_O+n_O, idx_O+2*n_O and idx_O+3*n_O a watermolecules?
#   --> we exploit this order down the line, so make sure its OK
#       what can happen is that PBC separate water molecules, in that case: remove them on the fly
for idx, atom in enumerate(file_poscar.list_atoms):
    if atom[file_poscar.ATOMNAME] == "O":

      # --- distance to HW1 ---
      d = pbc_distance( [ atom[file_poscar.X_COORD],                              \
                          atom[file_poscar.Y_COORD],                              \
                          atom[file_poscar.Z_COORD] ],                            \
                        [ file_poscar.list_atoms[idx+n_O][file_poscar.X_COORD],   \
                          file_poscar.list_atoms[idx+n_O][file_poscar.Y_COORD],   \
                          file_poscar.list_atoms[idx+n_O][file_poscar.Z_COORD] ], \
                        file_poscar.lattice_vector_1,                             \
                        file_poscar.lattice_vector_2,
                        file_poscar.lattice_vector_3) 

      if d > 1.5:
        print "[ERROR]. O ( idx:", idx, ") and H1 ( idx:", n_O+idx, ") are expected to belong to the same water molecule."
        print "         They are too far away from each other however, d=", d
        sys.exit()
      
      # --- distance to HW2 ---
      d = pbc_distance( [ atom[file_poscar.X_COORD],                                \
                          atom[file_poscar.Y_COORD],                                \
                          atom[file_poscar.Z_COORD] ],                              \
                        [ file_poscar.list_atoms[idx+2*n_O][file_poscar.X_COORD],   \
                          file_poscar.list_atoms[idx+2*n_O][file_poscar.Y_COORD],   \
                          file_poscar.list_atoms[idx+2*n_O][file_poscar.Z_COORD] ], \
                        file_poscar.lattice_vector_1,                               \
                        file_poscar.lattice_vector_2,
                        file_poscar.lattice_vector_3) 

      if d > 1.5:
        print "[ERROR]. O ( idx:", idx, ") and H2 ( idx:", 2*n_O+idx, ") are expected to belong to the same water molecule."
        print "         They are too far away from each other however, d=", d
        sys.exit()
      
      # --- distance to MW ---
      d = pbc_distance( [ atom[file_poscar.X_COORD],                                \
                          atom[file_poscar.Y_COORD],                                \
                          atom[file_poscar.Z_COORD] ],                              \
                        [ file_poscar.list_atoms[idx+3*n_O][file_poscar.X_COORD],   \
                          file_poscar.list_atoms[idx+3*n_O][file_poscar.Y_COORD],   \
                          file_poscar.list_atoms[idx+3*n_O][file_poscar.Z_COORD] ], \
                        file_poscar.lattice_vector_1,                               \
                        file_poscar.lattice_vector_2,
                        file_poscar.lattice_vector_3) 

      if d > 0.5:
        print "[ERROR]. O ( idx:", idx, ") and MW ( idx:", 3*n_O+idx, ") are expected to belong to the same water molecule."
        print "         They are too far away from each other however, d=", d
        sys.exit()

##########################################################################################
#                            MODIFY STRUCTURE                                            # 
##########################################################################################  
"""
 depending on the surface cut, the bottom layer might have not all SOL molecules in there
 because of PBC (in z)

 here we perform the following modifications to structure:
  (1) crop away unnecessary/incomplete bottom layer
  (2) shift the lowest lying OW to z=0
"""

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
#                            (1) crop away unwanted bottom layer                        #  
#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
list_remove = []

if flag_z_layer_remove in ["True", "true", "yes", "Yes", "1"]:
  print "\n[FLAG_Z_LAYER_REMOVE] ON. Removing all OW below", Z_LAYER_REMOVE
  print "\n"

  nO_remove = 0
  for idx, atom in enumerate(file_poscar.list_atoms):
    if atom[file_poscar.ATOMNAME] == "O":
      if atom[file_poscar.Z_COORD] < Z_LAYER_REMOVE:
        nO_remove += 1
                           
        list_remove.extend([idx, n_O+idx, 2*n_O+idx, 3*n_O+idx]) # we checked previously that this order is OK

  
  # order list
  list_remove_sorted = sorted(list_remove)

  # remove elements
  for idx in reversed(list_remove_sorted):
    del file_poscar.list_atoms[idx]

  n_O -= nO_remove
  n_H -= nO_remove*2
  n_X -= nO_remove

  if 4*n_O != len(file_poscar.list_atoms):
    print "[ERROR]. Could not remove Z_LAYER successfully."
    print "         The number of atoms makes no sense anymore"
    print "         expected:", 4*n_O, "found:", len(file_poscar.list_atoms)
    sys.exit()

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
#                            (2) shift lowest lying OW to z=0                           #  
#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
z_min = 999999999

for atom in file_poscar.list_atoms:
  if atom[file_poscar.ATOMNAME] == "O" and atom[file_poscar.Z_COORD] < z_min:
    z_min = atom[file_poscar.Z_COORD]

for atom in file_poscar.list_atoms:
  atom[file_poscar.Z_COORD] -= z_min

Z_LATTICE_CUTOFF -= z_min # this value comes from original POSCAR and needs to be changed accordingly

##########################################################################################
#                            DEFINE CROP REGION                                          # 
##########################################################################################  
crop_region = []
crop_region.append([0.0,0.0,0.0]) 

# --- figure out radius that brings you closest to desired N_SOL ---
"""
 (1) figure out a reasonable start value for the radius search
      in continuos case: N_SOL = 4/6*pi*r**3*rho_N
       with rho_N the number density of SOL molecules, we actively underestimate it slightly here
        rho_N ... 0.03

 (2) scan through values of stepsize, minimize deviation N_SOL - N_TARGET
"""
stepsize = 0.05


r_i = np.power( nSOL_target / (4.0/6.0 * np.pi * 0.03 ), 1.0/3.0 )

# find best R
FLAG_CONVERGED = False
nSOL_prev = -1
iteration = 0

print "******************************************************"
print "* FINDING OPTIMAL RADIUS TO ACCOMODATE", nSOL_target, "MOLECULES *"
print "******************************************************"

print "\ni N_SOL R"

while FLAG_CONVERGED == False:
  iteration += 1
  nSOL = 0

  # find out maximum size of supercell to create hemispherical cap of radius r_i
  n_x = int(np.ceil(r_i / file_poscar.lattice_vector_1[0]))
  n_y = int(np.ceil(r_i / file_poscar.lattice_vector_2[1]))

  # find number of SOL within r_i
  for idx, atom in enumerate(file_poscar.list_atoms):
    if atom[file_poscar.ATOMNAME] == "O":
      for i_x in range(-n_x, n_x):
        for i_y in range(-n_y, n_y):  
          if ( (atom[file_poscar.X_COORD] + i_x*file_poscar.lattice_vector_1[0] + i_y*file_poscar.lattice_vector_2[0] - crop_region[CENTER][0])**2 + \
               (atom[file_poscar.Y_COORD] + i_x*file_poscar.lattice_vector_1[1] + i_y*file_poscar.lattice_vector_2[1] - crop_region[CENTER][1])**2 + \
               (atom[file_poscar.Z_COORD] - crop_region[CENTER][2])**2 < \
               r_i**2):
            
            # count the total number of water molecules in lattice+crystallite:
            nSOL += 1

  print iteration, nSOL, r_i

  # --- exact match nSOL - target ---
  if nSOL == nSOL_target:
    FLAG_CONVERGED = True
    r = r_i
  
  # --- hemisphere too big ---
  elif nSOL > nSOL_target:
    # did we just find the boundary of r?
    if nSOL_prev < nSOL_target and nSOL_prev > 0:
      diff_prev = np.abs(nSOL_target-nSOL_prev)
      diff_now  = np.abs(nSOL_target-nSOL)

      if diff_prev <= diff_now:
        r = r_i - stepsize
      else:
        r = r_i

      FLAG_CONVERGED = True
    # keep looking
    else:
      nSOL_prev = nSOL
      r_i -= stepsize

  # --- hemisphere too small ---
  elif nSOL < nSOL_target:
    # did we just find the boundary of r?
    if nSOL_prev > nSOL_target and nSOL_prev > 0:
      diff_prev = np.abs(nSOL_target-nSOL_prev)
      diff_now  = np.abs(nSOL_target-nSOL)

      if diff_prev <= diff_now:
        r = r_i + stepsize
      else:
        r = r_i

      FLAG_CONVERGED = True
    # keep looking
    else:
      nSOL_prev = nSOL
      r_i += stepsize

print "[CONVERGENCE REACHED]\n"

crop_region.append(r)

##########################################################################################
#                            HEMISPHERE-CROP                                             # 
##########################################################################################  
"""
we only look for OW molecules inside hemisphere, the HW and MW will be added at the end

also screen through periodic images, because center is often at 0/0/0!
 ==> no periodicity in z because it is a slab!
"""
list_molecules_cryst   = []     # stores the atoms belonging to crystallite (gro)
list_molecules_lattice = []     # stores the atoms belonging to lattice

nSOL = 0 # total number of SOL molecules

# find out maximum size of supercell to create hemispherical cap
n_x = int(np.ceil(crop_region[RADIUS] / file_poscar.lattice_vector_1[0]))
n_y = int(np.ceil(crop_region[RADIUS] / file_poscar.lattice_vector_2[1]))

# check on the fly: is the water slab used to crop hemisphere thick enough?
z_O_max = -99999999
for idx, atom in enumerate(file_poscar.list_atoms):
  if atom[file_poscar.ATOMNAME] == "O":
    if atom[file_poscar.Z_COORD] > z_O_max:
      z_O_max = atom[file_poscar.Z_COORD]
if z_O_max - crop_region[CENTER][2] < crop_region[RADIUS]:
  print "[ERROR]. The slab is not thick (z) enough to be able to crop out hemisphere!"
  print "         max height (O):", z_O_max, "radius:", crop_region[RADIUS], "center[Z]:", crop_region[CENTER][2] 
  sys.exit()

for idx, atom in enumerate(file_poscar.list_atoms):
  if atom[file_poscar.ATOMNAME] == "O":
    for i_x in range(-n_x, n_x):
      for i_y in range(-n_y, n_y):  
        if ( (atom[file_poscar.X_COORD] + i_x*file_poscar.lattice_vector_1[0] + i_y*file_poscar.lattice_vector_2[0] - crop_region[CENTER][0])**2 + \
             (atom[file_poscar.Y_COORD] + i_x*file_poscar.lattice_vector_1[1] + i_y*file_poscar.lattice_vector_2[1] - crop_region[CENTER][1])**2 + \
             (atom[file_poscar.Z_COORD] - crop_region[CENTER][2])**2 < \
             crop_region[RADIUS]**2):
          
          # count the total number of water molecules in lattice+crystallite:
          nSOL += 1

          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
          #                            create lattice file                                        #  (center around origin)
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
          if atom[file_poscar.Z_COORD] < Z_LATTICE_CUTOFF:
            list_molecules_lattice.append( \
                                           [ "OW", atom[file_poscar.X_COORD] + i_x*file_poscar.lattice_vector_1[0] + i_y*file_poscar.lattice_vector_2[0] - crop_region[CENTER][0], \
                                                   atom[file_poscar.Y_COORD] + i_x*file_poscar.lattice_vector_1[1] + i_y*file_poscar.lattice_vector_2[1] - crop_region[CENTER][1], \
                                                   atom[file_poscar.Z_COORD] - crop_region[CENTER][2] ] )

          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
          #                            create crystallite (gro)                                   #  
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
          else:
            molecule = [ \
                         [ "OW", atom[file_poscar.X_COORD] + i_x*file_poscar.lattice_vector_1[0] + i_y*file_poscar.lattice_vector_2[0] - crop_region[CENTER][0],                          \
                                 atom[file_poscar.Y_COORD] + i_x*file_poscar.lattice_vector_1[1] + i_y*file_poscar.lattice_vector_2[1] - crop_region[CENTER][1],                          \
                                 atom[file_poscar.Z_COORD] - crop_region[CENTER][2]],                                                                                                     \
                         [ "HW1", \
                           file_poscar.list_atoms[n_O+idx][file_poscar.X_COORD] + i_x*file_poscar.lattice_vector_1[0] + i_y*file_poscar.lattice_vector_2[0] - crop_region[CENTER][0],     \
                           file_poscar.list_atoms[n_O+idx][file_poscar.Y_COORD] + i_x*file_poscar.lattice_vector_1[1] + i_y*file_poscar.lattice_vector_2[1] - crop_region[CENTER][1],     \
                           file_poscar.list_atoms[n_O+idx][file_poscar.Z_COORD] - crop_region[CENTER][2]                                                                             ],   \
                         [ "HW2", \
                           file_poscar.list_atoms[2*n_O+idx][file_poscar.X_COORD] + i_x*file_poscar.lattice_vector_1[0] + i_y*file_poscar.lattice_vector_2[0] - crop_region[CENTER][0],   \
                           file_poscar.list_atoms[2*n_O+idx][file_poscar.Y_COORD] + i_x*file_poscar.lattice_vector_1[1] + i_y*file_poscar.lattice_vector_2[1] - crop_region[CENTER][1],   \
                           file_poscar.list_atoms[2*n_O+idx][file_poscar.Z_COORD] - crop_region[CENTER][2]                                                                             ], \
                         [ "MW", \
                           file_poscar.list_atoms[3*n_O+idx][file_poscar.X_COORD] + i_x*file_poscar.lattice_vector_1[0] + i_y*file_poscar.lattice_vector_2[0] - crop_region[CENTER][0],   \
                           file_poscar.list_atoms[3*n_O+idx][file_poscar.Y_COORD] + i_x*file_poscar.lattice_vector_1[1] + i_y*file_poscar.lattice_vector_2[1] - crop_region[CENTER][1],   \
                           file_poscar.list_atoms[3*n_O+idx][file_poscar.Z_COORD] - crop_region[CENTER][2]                                                                             ]  ]

            # check for PBC within molecule, if necessary project to OW
            for a in molecule:
              if a[0] == "OW":
                ref_x = a[1]
                ref_y = a[2]
                ref_z = a[3]
              else:
                d = pbc_vector( [a[1], a[2], a[3]],           \
                                [ref_x, ref_y, ref_z],        \
                                file_poscar.lattice_vector_1, \
                                file_poscar.lattice_vector_2, \
                                file_poscar.lattice_vector_3  )

                a[1] = ref_x + d[0]
                a[2] = ref_y + d[1]
                a[3] = ref_z + d[2]

            list_molecules_cryst.append(molecule)

##########################################################################################
#                            GENERATE GROMACS FILE (CRYSTALLITE)                         # 
##########################################################################################  
# only generate gro if there are points of hemisphere that are not included in lattice
if list_molecules_cryst != []:
    file_gro = format_GRO()
    
    # it is a crystallite --> no more pbc necessary
    file_gro.lattice_vector_1 = [i/10.0 for i in file_poscar.lattice_vector_1]
    file_gro.lattice_vector_2 = [i/10.0 for i in file_poscar.lattice_vector_2]
    file_gro.lattice_vector_3 = [i/10.0 for i in file_poscar.lattice_vector_3]
    
    idx=0
    for mol_num, mol in enumerate(list_molecules_cryst):
      for atom in mol:
        idx += 1
        file_gro.list_atoms.append([mol_num, "SOL", atom[0], idx, atom[1]/10.0, atom[2]/10.0, atom[3]/10.0])
    
    file_gro.n_atoms=len(file_gro.list_atoms)
    
    # the titel in the gro-file is used to link the crystallite with the lattice
    # the link is made via index, linking the first atom in list of both

    dx = list_molecules_cryst[0][0][1] - list_molecules_lattice[0][1]
    dy = list_molecules_cryst[0][0][2] - list_molecules_lattice[0][2]
    dz = list_molecules_cryst[0][0][3] - list_molecules_lattice[0][3]
    
    file_gro.titel = "MODE_LOCK_INDEX 0 0 " + str(dx) + " " + str(dy) + " " + str(dz)
    
    file_gro.write_file(outfilename_gro)

##########################################################################################
#                            GENERATE LATTICE FILE                                       # 
##########################################################################################  
"""
 generate a generic lattice file here, with the following properties:

  (+) 360 rotation
  (+) _X_, _Y_, _Z_ translation (depends on the symmetry of the underlying substrate)
  (+) keep OW at exact lattice positions
  (+) wiggle space in z: +- 0.5 AA
  (+) linkage to crystallite
"""
file_lattice = open(outfilename_lattice, "w")

comment = "# lattice file generated by " + sys.argv[0] + " using " + infilename_poscar + "\n"
file_lattice.write(comment)

comment = "#  ==> generate crystallite and interlocked lattice\n\n"
file_lattice.write(comment)

comment = "# n_latticepoints\n"
file_lattice.write(comment)
file_lattice.write(str(len(list_molecules_lattice)) + "\n\n")

comment = "# specify the lattice [Angstrom]\n"
file_lattice.write(comment)

for idx,lp in enumerate(list_molecules_lattice):
  file_lattice.write(str(idx) + " " + str(lp[1]) + " " + str(lp[2]) + " " + str(lp[3]) + "\n")
    
file_lattice.write("\n")

comment = "# n_fixed_points\n0\n\n"
file_lattice.write(comment)

comment = "# fixed points\n\n"
file_lattice.write(comment)

comment = "# center of lattice\n# always at origin with this tool\n0.0 0.0 0.0\n\n"
file_lattice.write(comment)

comment = "# replicas\n# they are useless here because its only one nonperiodic seed\n"
file_lattice.write(comment)
comment = "COMPUTE 1 1  100.0 0.0 0.0    0.0 100.0 0.0     0.0 0.0 100.0\n\n"
file_lattice.write(comment)

comment = "# target position [Angstrom]\nTARGET_X TARGET_Y TARGET_Z\n\n"
file_lattice.write(comment)

comment = "# rotation angle of lattice around z [degree]\n360.0\n\n"
file_lattice.write(comment)

comment = "# displacement of WHOLE lattice points [Angstrom]\n#   x, y, z\n"
file_lattice.write(comment)
comment = "DISPLACE_X DISPLACE_Y DISPLACE_Z\n\n"
file_lattice.write(comment)

comment = "# displacement around each individual lattice point [Angstrom]\n#   xy,z\n0.0 0.0\n\n"
file_lattice.write(comment)

comment = "# scale lattice\n#   keep positions exactly as they are\nFalse\n"
file_lattice.write(comment)

file_lattice.close()

##########################################################################################
#                            SUMMARY                                                     # 
########################################################################################## 
print "R   N_SOL   N_LATTICE"
print crop_region[RADIUS], nSOL, len(list_molecules_lattice)
