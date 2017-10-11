#!/usr/bin/python

"""
This code produces a random structure of an adsorbate on a given
substrate

(+) INPUT: how many water do you wanna put on SF
(+) INPUT: substrate-structure (*gro file)
(+) adsorbate: water (SPC model)
(+) randomly rotate water molecule (angle, axis)
(+) three ways of positioning water molecules on slab
    (a) completely random: 
          for high coverages, chances are high that several water molecules overlap
          --> minimized structure will be crap
    (b) along pre-defined lattice: 
          user defines basic lattice, onto each lattice-point a water molecule will be placed
          (randomly rotated watermolecule, displaced a bit from latticepoint)
          --> this way water molecules will barely overlap, most of the structures will be more
          sensible than (1)
    (c) along custom lattice specified in file (filename_custom_lattice)
          user creates file for that. This is a good approach for water clusters on the surface
          ==> add rotation along z axis, since cluster should be able to rotate
"""

__author__    = "Pedevilla Philipp"
__copyright__ = "all rights reserved. 2014"
__status__    = "Production"

##########################################################################################
#                            libraries                                                   # 
##########################################################################################  
import math
import copy              # to copy lists
import sys               # to exit
import os.path           # to check if file exists

import numpy as np       # invert matrix
import random            # to generate random number

# --- my own libraries ---
import transform_coordinates as tc         # convert frac2cart / cart2frac
import construct_lattice as cl             # construct lattice for adsorbate positioning
import rotations as rot                    # construct roatation matrix
import file_formats as ff                  # read files [for laziness, not the entire code is using this]
from pbc_computations import pbc_distance  # pbc_distance calculation

##########################################################################################
#                            MAKROS                                                      # 
##########################################################################################  
coord_X = 0
coord_Y = 1
coord_Z = 2

# --- for structure atom ---
INDEX  = 3
X_CART = 4
Y_CART = 5
Z_CART = 6
X_FRAC = 7
Y_FRAC = 8
Z_FRAC = 9

# --- latticetypes ---
RECTANGULAR = 0
HEXAGONAL   = 1
IH_BASAL    = 2
IH_PRISM    = 3
CUSTOM      = 4 

# --- lattice positioning ---
RANDOM   = 0
CENTERED = 1

# --- others ---
NO  = 0
YES = 1

# --- lattice replica type ---
EXPLICIT = 1  # predefined positions for lattice replicas (e.g. if irregular)
COMPUTE  = 2  # calculate lattice replicas on the run

# --- occupy lattice points at random or equally ---
RANDOM_OCCUPATION = 0 # default
EQUAL_OCCUPATION  = 1

# --- how to fill lattices ---  (if more lattice points than water molecules)
FILL_LATTICE_1_FULLY = 1
FILL_LATTICE_2_FULLY = 2

# --- is lattice filled completely --- (then lattice points are not picked randomly to be filled)
FILLED_PARTIALLY  = 0
FILLED_COMPLETELY = 1

# --- water models ---
SPC       = 0
TIP4P_ICE = 1 # different TIP4P models have different r(OM)

##########################################################################################
#                  INPUT PARAMETER                                                       #
##########################################################################################

# (1) lattice_type
# (2) extra_displacement
# (3) flag_scale_lattice
# (4) multiply cell (in x,y and z)
# (5) custom lattice filename

#---------------------#
#     LATTICETYPE     #
#---------------------#
lattice_type = CUSTOM

#-------------------------#
#     LATTICE SCALING     #
#-------------------------#
# should the lattice be scaled to unit cell dimension
flag_scale_lattice = True

#-----------------------#
#     MULTIPLY_CELL     #
#-----------------------#
# for FF minimization
# [test]  multiply_x = 1   
# [test]  multiply_y = 1   
# [test]  multiply_z = 1   

# [001 microcline]   multiply_x = 3  
# [001 microcline]   multiply_y = 2  
# [001 microcline]   multiply_z = 1  

# [010 microcline]  multiply_x = 3    
# [010 microcline]  multiply_y = 4    
# [010 microcline]  multiply_z = 1    

# [100 microcline] multiply_x = 2    
# [100 microcline] multiply_y = 3    
# [100 microcline] multiply_z = 1    
# [100 microcline - melt-MD] multiply_x = 3    
# [100 microcline - melt-MD] multiply_y = 4    
# [100 microcline - melt-MD] multiply_z = 1    

# [-1-10 microcline]   multiply_x = 2    
# [-1-10 microcline]   multiply_y = 4    
# [-1-10 microcline]   multiply_z = 1    

multiply_x = 1   # kaolinite seed test  
multiply_y = 1   # kaolinite seed test  
multiply_z = 1   # kaolinite seed test  


#-----------------------------------------#
#     USE LATTICE / RANDOM GENERATION     #
#-----------------------------------------#
flag_lattice = True       # use lattice
#flag_lattice = False     # random placement

#---------------------------------#
#     CUSTOM_LATTICE FILENAME     #
#---------------------------------#
# *** test (empty slab) ***
#filename_custom_lattice = ["ice-Ih-basal_1L.lattice"]
#filename_custom_lattice = ["ice-Ih-basal_1L.lattice", "ice-Ih-basal_2L.lattice"]
#filename_custom_lattice = ["ice-Ih-prism-100_1L.lattice"]
#filename_custom_lattice = ["ice-Ih-prism-100_1L.lattice", "ice-Ih-prism-100_2L.lattice"]
#filename_custom_lattice = ["ice-Ih-prism-110_1L.lattice"]
#filename_custom_lattice = ["ice-Ih-prism-110_1L.lattice", "ice-Ih-prism-110_2L.lattice"]
#filename_custom_lattice = ["ice-Ih-prism-110_COMBINED_1L.lattice", "ice-Ih-prism-110_COMBINED_2L.lattice"]

# *** microcline (001) search ***
#filename_custom_lattice = ["ice-Ih-prism-100_1L.lattice"]

# *** microcline (010) search ***
#filename_custom_lattice = ["10_lattices-010/ice-Ih-basal_1L.lattice"]
#filename_custom_lattice = ["10_lattices-010/ice-Ih-prism-100_1L.lattice"]
#filename_custom_lattice = ["10_lattices-010/ice-Ih-prism-110_COMBINED_1L.lattice"]
#filename_custom_lattice = ["10_lattices_010/ice-Ih-prism-110_COMBINED_1L.lattice", "10_lattices_010/ice-Ih-prism-110_COMBINED_2L.lattice"]
#filename_custom_lattice = ["10_lattices_010/ice-Ih-basal_1L.lattice", "10_lattices_010/ice-Ih-basal_2L.lattice"]
#filename_custom_lattice = ["10_lattices_010/ice-Ih-prism-100_1L.lattice","10_lattices_010/ice-Ih-prism-100_2L.lattice"]


# *******************************
# *** microcline (100) search ***
# *******************************
#filename_custom_lattice = ["10a_100_lattices/ice-Ih-prism-100_1L.lattice"]
#filename_custom_lattice = ["10a_100_lattices/ice-basal_100_1L.lattice"]
#filename_custom_lattice = ["10a_100_lattices/rectangular_100.lattice"]
#filename_custom_lattice = ["10a_100_lattices/monomer_scan.lattice"]
#filename_custom_lattice = ["12_lattices_100/ice-basal_100_1L.lattice", "12_lattices_100/ice-basal_100_2L.lattice"]
#filename_custom_lattice = ["12_lattices_100/ice-Ih-prism-100_1L.lattice", "12_lattices_100/ice-Ih-prism-100_2L.lattice"]
#filename_custom_lattice = ["12_lattices_100/ice-Ih-prism-110_COMBINED_1L.lattice"]
#filename_custom_lattice = ["12_lattices_100/ice-Ih-prism-110_COMBINED_1L.lattice","12_lattices_100/ice-Ih-prism-110_COMBINED_2L.lattice"]
#filename_custom_lattice = ["12_lattices_100/ice-Ih-prism-100_4L_COMBINED.lattice"]

# *********************************
# *** microcline (-1-10) search ***
# *********************************
#filename_custom_lattice = ["11_lattices_-1-10/ice-Ih-basal_1L.lattice"]
#filename_custom_lattice = ["11_lattices_-1-10/ice-Ih-prism-100_1L.lattice"]
#filename_custom_lattice = ["11_lattices_-1-10/ice-Ih-prism-110_COMBINED_1L.lattice"]
#filename_custom_lattice = ["11_lattices_-1-10/ice-Ih-prism-110_COMBINED_1L.lattice", "11_lattices_-1-10/ice-Ih-prism-110_COMBINED_2L.lattice"]
#filename_custom_lattice = ["11_lattices_-1-10/ice-Ih-basal_1L.lattice","11_lattices_-1-10/ice-Ih-basal_2L.lattice"]
#filename_custom_lattice = ["11_lattices_-1-10/ice-Ih-prism-100_1L.lattice", "11_lattices_-1-10/ice-Ih-prism-100_2L.lattice"]
#filename_custom_lattice = ["11_lattices_-1-10/ice-Ih-prism-110_COMBINED_8L.lattice"]

##########################################################################################
#                            KAOLINITE                                                   # 
##########################################################################################  
# ************************************
# *** kaolinite basal face seeding ***
# ************************************
# --- box ---
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/1L/basal_small.lattice"]    # debug purpose
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/1L/basal_orthogonal_big.lattice"] # debug purpose
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/1L/basal_massive.lattice"] # debug purpose
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/1L/basal_orthogonal_big_shifted.lattice"] # debug purpose
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/5L/basal_orthogonal_big_shifted.lattice"] # seed lattice
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/spherical_cap/basal_massive_5L.lattice"] # big lattice from which we carve spherical cap
# --- spherical caps ---
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/spherical_cap/basal_spherecap_R-10.lattice"] # spherical cap, radius 10 AA, positions fixed
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/spherical_cap/basal_spherecap_R-10_rot-trans.lattice"] # spherical cap, radius 10 AA, positions flexible
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/spherical_cap/basal_spherecap_R-15.lattice"] # spherical cap, radius 15 AA
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/spherical_cap/basal_spherecap_R-15_rot-trans.lattice"] # spherical cap, radius 15 AA, positions flexible
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/spherical_cap/basal_spherecap_R-20.lattice"] # spherical cap, radius 20 AA
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/spherical_cap/basal_spherecap_R-20_rot-trans.lattice"] # spherical cap, radius 20 AA, positions flexible
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/spherical_cap/basal_spherecap_R-25.lattice"] # spherical cap, radius 25 AA
#filename_custom_lattice = ["13_seeds_kaolinite/01_basal_plane/spherical_cap/basal_spherecap_R-25_rot-trans.lattice"] # spherical cap, radius 25 AA, positions flexible

# ************************************
# *** kaolinite prism face seeding ***
# ************************************
# --- box ---
#filename_custom_lattice = ["13_seeds_kaolinite/02_prism_plane/1L/prism_small.lattice"] # debug purpose
#filename_custom_lattice = ["13_seeds_kaolinite/02_prism_plane/1L/prism_big.lattice"] # debug purpose
#filename_custom_lattice = ["13_seeds_kaolinite/02_prism_plane/1L/prism_massive.lattice"] # debug purpose, make superlattice
#filename_custom_lattice = ["13_seeds_kaolinite/02_prism_plane/1L/kao_prism_full_seed.lattice"] # debug purpose
#filename_custom_lattice = ["13_seeds_kaolinite/02_prism_plane/8L/kao_prism_full_seed_8L.lattice"] # seed lattice

# --- spherical cap ---
#filename_custom_lattice = ["13_seeds_kaolinite/02_prism_plane/spherical_cap/kao_prism_full_coverage.lattice"] # generate massive overlayer
#filename_custom_lattice = ["13_seeds_kaolinite/02_prism_plane/spherical_cap/prism_spherecap_R-10_rot-trans.lattice"] # spherical cap, radius 10 AA, positions flexible
#filename_custom_lattice = ["13_seeds_kaolinite/02_prism_plane/spherical_cap/prism_spherecap_R-15_rot-trans.lattice"] # spherical cap, radius 15 AA, positions flexible
#filename_custom_lattice = ["13_seeds_kaolinite/02_prism_plane/spherical_cap/prism_spherecap_R-20_rot-trans.lattice"] # spherical cap, radius 20 AA, positions flexible
#filename_custom_lattice = ["13_seeds_kaolinite/02_prism_plane/spherical_cap/prism_spherecap_R-25_rot-trans.lattice"] # spherical cap, radius 25 AA, positions flexible

##########################################################################################
#                            CHOLESTEROL                                                 # 
##########################################################################################  
# **************************
# *** cholesterol c0_213 *** (flat cholesterol phase)
# **************************
#filename_custom_lattice = ["14b_cholesterol-lattices/C0/ice-Ih-prism-100_1L.lattice"]
#filename_custom_lattice = ["14b_cholesterol-lattices/C0/ice-Ih-prism-100_4L.lattice"] # generate ice structure for collision mode cropping
#filename_custom_lattice = ["14b_cholesterol-lattices/C0/ice-Ih-prism-100_CROPPED-4L.lattice"]
#filename_custom_lattice = ["14b_cholesterol-lattices/C0/prim-pris_separated_lattices/ice-Ih-prism-100_CROPPED-CL.lattice", \
#                           "14b_cholesterol-lattices/C0/prim-pris_separated_lattices/ice-Ih-prism-100_CROPPED-OL.lattice"]
#filename_custom_lattice = ["14b_cholesterol-lattices/C0/prim-pris_separated_lattices/ice-Ih-prism-100_CROPPED-CL.lattice"] # CL
#filename_custom_lattice = ["14b_cholesterol-lattices/C0/prim-pris_separated_lattices/ice-Ih-prism-100_CROPPED-OL.lattice"] # create OL on empty slab

# ***************************
# *** cholesterol cm0_331 *** (flat cholesterol phase)
# ***************************
# --- Ih lattice, 150 SOL ---
#filename_custom_lattice = ["16_seeds_cholesterol/02_HEMISPHERES/cm0/Ih/01_basal/150_SOL/Ih_001.lattice"]
#filename_custom_lattice = ["16_seeds_cholesterol/02_HEMISPHERES/cm0/Ih/02_primary_prism/150_SOL/Ih_100.lattice"]
#filename_custom_lattice = ["16_seeds_cholesterol/02_HEMISPHERES/cm0/Ih/03_secondary_prism/150_SOL/Ih_110.lattice"]

# --- Ic lattice, 150 SOL ---
#filename_custom_lattice = ["16_seeds_cholesterol/02_HEMISPHERES/cm0/Ic/150_SOL/01_001/Ic_001_150-SOL.lattice"]
#filename_custom_lattice = ["16_seeds_cholesterol/02_HEMISPHERES/cm0/Ic/150_SOL/02_111/Ic_111_150-SOL.lattice"]

# **************************
# *** cholesterol c0_231 *** (buckled cholesterol phase)
# **************************
# --- Ih lattice, 150 SOL ---
#filename_custom_lattice = ["16_seeds_cholesterol/02_HEMISPHERES/c0/Ih/01_basal/150_SOL/Ih_001.lattice"]
#filename_custom_lattice = ["16_seeds_cholesterol/02_HEMISPHERES/c0/Ih/02_primary_prism/150_SOL/Ih_100.lattice"]
#filename_custom_lattice = ["16_seeds_cholesterol/02_HEMISPHERES/c0/Ih/03_secondary_prism/150_SOL/Ih_110.lattice"]

#filename_custom_lattice = ["16_seeds_cholesterol/02_HEMISPHERES/c0/Ic/001/150_SOL/Ic_001.lattice"]
#filename_custom_lattice = ["16_seeds_cholesterol/02_HEMISPHERES/c0/Ic/111/150_SOL/Ic_111.lattice"]

##########################################################################################
#                            MDHE                                                        # 
##########################################################################################  
# **************************
# *** MDHE testing       *** 
# **************************
#filename_custom_lattice = ["17_gabri_test_MDHE/Ih_001/test_150.lattice"]
filename_custom_lattice = ["./seed_Ic_001.lattice"]


##########################################################################################
#                            ONLY SEEDS                                                  # 
##########################################################################################  
# ************************
# *** only ice - seeds *** [ boxes ]
# ************************
#filename_custom_lattice = ["15_seeds_general/01_boxes/basal/basal_box_crystallite_080-OW.lattice"] # BASAL small box
#filename_custom_lattice = ["15_seeds_general/01_boxes/primary_prism/prim_pris_box_084-OW.lattice"] # PRIMARY PRISM small box
#filename_custom_lattice = ["15_seeds_general/01_boxes/secondary_prism/sec_pris_box_080-OW.lattice"] # SECONDARY PRISM small box

# --- shaped crystallites made with VESTA Ih ---
# [OK] filename_custom_lattice = ["15_seeds_general/00_homogeneous_crystallites/Ih/Ice-Ih_001_100_flat_648-O.lattice"]
# [OK] filename_custom_lattice = ["15_seeds_general/00_homogeneous_crystallites/Ih/Ice-Ih_001_100_elongated_648-O.lattice"]
# [OK] filename_custom_lattice = ["15_seeds_general/00_homogeneous_crystallites/Ih/Ice-Ih_001_100_similar_672-O.lattice"]

#filename_custom_lattice = ["15_seeds_general/00_homogeneous_crystallites/Ih/Ice-Ih_001_100_100_elongated_672-O.lattice"]
#filename_custom_lattice = ["15_seeds_general/00_homogeneous_crystallites/Ih/Ice-Ih_001_100_110_flat_612-O.lattice"]
#filename_custom_lattice = ["15_seeds_general/00_homogeneous_crystallites/Ih/Ice-Ih_001_100_100_similar_672-O.lattice"]


# --- shaped crystallites made with VESTA Ic ---
#filename_custom_lattice = ["15_seeds_general/00_homogeneous_crystallites/Ic/Ice_Ic_100_111_597-O.lattice"]
#filename_custom_lattice = ["15_seeds_general/00_homogeneous_crystallites/Ic/Ice_Ic_100_621-O.lattice"]
#filename_custom_lattice = ["15_seeds_general/00_homogeneous_crystallites/Ic/Ice_Ic_111_567-O.lattice"]
   
# *************
# *** DEBUG ***
# *************
#filename_custom_lattice = ["basal_COLLISION_TEST.lattice"]
#filename_custom_lattice = ["prism_COLLISION_TEST.lattice"]


if len(filename_custom_lattice) > 2:
  print "[ERROR]. At the moment only two lattices are supported"
  sys.exit()

#----------------------------------#
#     LATTICE POINT OCCUPATION     #
#----------------------------------#
flag_lattice_point_occupation = RANDOM_OCCUPATION
#flag_lattice_point_occupation = EQUAL_OCCUPATION     # useful for monomer site scan 

#-----------------------------------#
#     LINK LATTICE DISPLACEMENT     #
#-----------------------------------#
flag_link_lattice_displacement = True  # link lattice displacement (xy + z) of lattice 1/2 (bas+bas, pris+pris)
#flag_link_lattice_displacement = False # displace lattices 1/2 independendly (xy + z)

#--------------------------------#
#     LATTICE FILLING METHOD     #
#--------------------------------#
# if there are two lattices, and not enough water molecules for all lattices
# --> decide how to fill lattice points (lattice 1 completely filled / lattice 2 completely filled?)
#flag_fill_lattice_method = FILL_LATTICE_1_FULLY
flag_fill_lattice_method = FILL_LATTICE_2_FULLY

#---------------------------------#
#     LATTICEPOINT OCCUPATION     #
#---------------------------------#
# if all lattice points will be filled, then they are not picked randomly
# but are instead filled one after another
#   this here is just the initialization, the actual value will be determined on the run
flag_lattice_filling = FILLED_PARTIALLY

#---------------------#
#     WATER MODEL     #
#---------------------#
#water_model = SPC
water_model = TIP4P_ICE

#----------------------------------------------------------#
#     REMOVE LATTICEPOINTS THAT COLLIDE WITH SUBSTRATE     #
#----------------------------------------------------------#
flag_REMOVE_COLLISIONS = True
#flag_REMOVE_COLLISIONS = False

#-------------------------------------#
#     LATTICE COLLISION TRESHHOLD     #
#-------------------------------------#
LATTICE_COLLISION_TRESHHOLD = 2.0 # in Angstrom

##########################################################################################
#                            check if args are correct                                   # 
##########################################################################################  
filename_substrate  = ""
filename_extra      = ""
n_h2o               = 0
filename_out        = ""
n_structures        = 0

# read arguments
if (len(sys.argv) == 5):
  filename_substrate  = sys.argv[1]
  n_h2o               = int(sys.argv[2])
  filename_out        = sys.argv[3]
  n_structures        = int(sys.argv[4])
elif len(sys.argv) == 6:
  filename_substrate  = sys.argv[1]
  filename_extra      = sys.argv[2]   # extra GRO 
  n_h2o               = int(sys.argv[3])
  filename_out        = sys.argv[4]
  n_structures        = int(sys.argv[5])
elif len(sys.argv) == 7:
  filename_substrate  = sys.argv[1]
  filename_extra      = sys.argv[2]   # extra GRO 
  n_h2o               = int(sys.argv[3])
  filename_out        = sys.argv[4]
  n_structures_start  = int(sys.argv[5])
  n_structures_end    = int(sys.argv[6])
  n_structures        = n_structures_end - n_structures_start
else:
  print "[ERROR]. Please specify all arguments (5 or 6 or 7) according to:"
  print "         %s FILENAME_SUBSTRATE.gro [FILENAME_EXTRA.gro] n_h2o FILENAME_OUTPUT n_structures [n_structures_start n_structures_end] " % sys.argv[0]
  sys.exit()

if (lattice_type == CUSTOM):
  if os.path.isfile(filename_custom_lattice[0]) == False:
    print "[ERROR]. The file " + filename_custom_lattice[0] + " does not exist"
    sys.exit()
  if len(filename_custom_lattice) == 2 and os.path.isfile(filename_custom_lattice[1]) == False:
    print "[ERROR]. The file " + filename_custom_lattice[1] + " does not exist"
    sys.exit()


##########################################################################################
#                            substrate (*gro format)                                     # 
##########################################################################################  
# old way of manually reading file in
#  now there is a class to handle this, change in the future for readability :-)
file_input = open(filename_substrate, "r")

list_slab_atoms = []
line_number   = 0
n_atoms_found = 0
atom_number   = 0

flag_empty_slab = NO

lattice_v_1 = [0.0, 0.0, 0.0]
lattice_v_2 = [0.0, 0.0, 0.0]
lattice_v_3 = [0.0, 0.0, 0.0]

for line in file_input:
  line_number += 1
  if line_number == 1:   # titel
    continue

  elif line_number == 2: # n_atoms
    n_atoms = int(line.split()[0])
    
    # check if empty slab
    if n_atoms == 0:
      flag_empty_slab = YES
  
  else:                     # coords/lattice_vector
    list_line = line.split()
    # --- lattice vector ---
    if n_atoms_found == n_atoms:
      lattice_v_1 = [ float(list_line[0]) , 
                      float(list_line[3]) , 
                      float(list_line[4]) ]
      
      lattice_v_2 = [ float(list_line[5]) , 
                      float(list_line[1]) , 
                      float(list_line[6]) ]
      
      lattice_v_3 = [ float(list_line[7]) , 
                      float(list_line[8]) , 
                      float(list_line[2]) ]

    # --- atoms ---
    else:
      # --- gro has a well defined format, read it accordingly ---
      # http://manual.gromacs.org/current/online/gro.html
      #   "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
      # ==> we check on the run if it is read correctly
      
      tmp = line[0:5].split()
      if len(tmp) != 1:
        print "[ERROR]. Could not read gro-file. The format is not correct (molecule number)"
        print line
        sys.exit()
      res_index     = int(tmp[0])  # molecule number


      tmp = line[5:10].split()
      if len(tmp) != 1:
        print "[ERROR]. Could not read gro-file. The format is not correct (molecule name)"
        print line
        sys.exit()
      res_name    = tmp[0]    # molecule name
      
      tmp = line[10:15].split()
      if len(tmp) != 1:
        print "[ERROR]. Could not read gro-file. The format is not correct (atom type)"
        print line
        sys.exit()
      atom_name   =  tmp[0]   # atom type
      
      tmp = line[15:20].split()
      if len(tmp) != 1:
        print "[ERROR]. Could not read gro-file. The format is not correct (atom number)"
        print line
        sys.exit()
      atom_number = int(tmp[0])   # atom number
      
      tmp = line[20:28].split()
      if len(tmp) != 1:
        print "[ERROR]. Could not read gro-file. The format is not correct (x-coord)"
        print line
        sys.exit()
      x_coord     = float(tmp[0])  
      
      tmp = line[28:36].split()
      if len(tmp) != 1:
        print "[ERROR]. Could not read gro-file. The format is not correct (y-coord)"
        print line
        sys.exit()
      y_coord     = float(tmp[0])
      
      tmp = line[36:44].split()
      if len(tmp) != 1:
        print "[ERROR]. Could not read gro-file. The format is not correct (z-coord)"
        print line
        sys.exit()
      z_coord     = float(tmp[0])
      
      list_slab_atoms.append([ res_index    , 
                               res_name     ,
                               atom_name    ,
                               atom_number  , 
                               x_coord      , 
                               y_coord      ,
                               z_coord      ])
      n_atoms_found += 1

##########################################################################################
#                            extra file (*gro format) [e.g. OL]                          # 
##########################################################################################  
if filename_extra:
  file_extra = ff.format_GRO() 
  file_extra.read_file(filename_extra)  # read trajectory

  # for consistency: make sure that the dimensions (x,y,z) are the same for 
  # substrate and extra lattice --> for now only support 1x1 supercells, change if necessary
  if multiply_x != 1 or multiply_y != 1:
    print "[ERROR]. At the moment, extra_structures are only supported for 1x1 supercells."
    print "         Change implementation :-)"
    sys.exit()

  if abs(lattice_v_1[0] - file_extra.lattice_vector_1[0]) > 0.000001 or \
     abs(lattice_v_1[1] - file_extra.lattice_vector_1[1]) > 0.000001 or \
     abs(lattice_v_1[2] - file_extra.lattice_vector_1[2]) > 0.000001:
    #print "[WARNING]. The lattices between substrate and extra_structure do not match up perfectly."
    #print "           lattice_vector_1 differs..."
    #print "              ", lattice_v_1, file_extra.lattice_vector_1
    #print "           Setting file_extra.lattice_vector_1 to lattice_v_1"
    file_extra.lattice_vector_1 = lattice_v_1

  elif abs(lattice_v_2[0] - file_extra.lattice_vector_2[0]) > 0.000001 or \
       abs(lattice_v_2[1] - file_extra.lattice_vector_2[1]) > 0.000001 or \
       abs(lattice_v_2[2] - file_extra.lattice_vector_2[2]) > 0.000001:
    #print "[WARNING]. The lattices between substrate and extra_structure do not match up perfectly."
    #print "           lattice_vector_2 differs..."
    #print "              ", lattice_v_2, file_extra.lattice_vector_2
    #print "           Setting file_extra.lattice_vector_2 to lattice_v_2"
    file_extra.lattice_vector_2 = lattice_v_2

  elif abs(lattice_v_3[0] - file_extra.lattice_vector_3[0]) > 0.000001 or \
       abs(lattice_v_3[1] - file_extra.lattice_vector_3[1]) > 0.000001 or \
       abs(lattice_v_3[2] - file_extra.lattice_vector_3[2]) > 0.000001:
    #print "[WARNING]. The lattices between substrate and extra_structure do not match up perfectly."
    #print "           lattice_vector_3 differs..."
    #print "              ", lattice_v_3, file_extra.lattice_vector_3
    #print "           Setting file_extra.lattice_vector_3 to lattice_v_3"
    file_extra.lattice_vector_3 = lattice_v_3

  #------------------------------------------------------#
  # parse information about displacement for this struct # (from titel)
  #------------------------------------------------------#
  list_modes = ["MODE_EXPLICIT", "MODE_LOCK_INDEX"]
  # info contained in title line

  if not file_extra.titel.split():
      print "[ERROR]. Titel in extra_file is empty, but needs to contain info about LOCK_MODE..."
      sys.exit()

  extra_file_mode = file_extra.titel.split()[0]

  # check if mode is known
  if extra_file_mode not in list_modes:
    print "[ERROR]. Could not recognize mode (", extra_file_mode, "). It has to be one of the following:"
    print list_modes
    sys.exit()

  if extra_file_mode == "MODE_EXPLICIT":
    if len(file_extra.titel.split()) != 6:
      print "[ERROR]. The title in", filename_extra, "cannot be parsed!"
      print "         It has to contain exactly these entries: MODE_EXPLICIT, displace_x, displace_y, displace_z, min_z_lock, max_z_lock [in Angstrom]"
      print file_extra.titel
      sys.exit()

    extra_file_displ_x = float(file_extra.titel.split()[1])
    extra_file_displ_y = float(file_extra.titel.split()[2])
    extra_file_displ_z = float(file_extra.titel.split()[3])
    min_z_lock         = float(file_extra.titel.split()[4]) # minimum distance between lattice and extra_structure
    max_z_lock         = float(file_extra.titel.split()[5]) # maximum distance between lattice and extra_structure

  elif extra_file_mode == "MODE_LOCK_INDEX":
    if len(file_extra.titel.split()) != 6:
      print "[ERROR]. The title in", filename_extra, "cannot be parsed!"
      print "         It has to contain exactly these entries: MODE_LOCK_INDEX, idx_lock_lattice, idx_lock_extra, lock_dx, lock_dy, lock_dz [in Angstrom]"
      print " --- Current ---"
      print file_extra.titel
      sys.exit()

    # info: the index in lattice + extra file that are locked with each other +
    #       the displacement (lattice-->extra) between them (in Angstrom!)
    idx_lock_lattice = int(file_extra.titel.split()[1])
    idx_lock_extra   = int(file_extra.titel.split()[2])
    lock_dx          = float(file_extra.titel.split()[3]) # rx(extra) - rx(lattice)
    lock_dy          = float(file_extra.titel.split()[4]) # ry(extra) - ry(lattice)
    lock_dz          = float(file_extra.titel.split()[5]) # rz(extra) - rz(lattice)

  #----------------------------------------------#
  # find lowest OW in extra structure for z_lock #
  #----------------------------------------------#
  z_min_extra_file = 9999999
  for atom in file_extra.list_atoms:
    if atom[file_extra.ATOMNAME] == "OW" and atom[file_extra.Z_COORD] < z_min_extra_file:
      z_min_extra_file = atom[file_extra.Z_COORD]

##########################################################################################
#                            adsorbate                                                   # 
##########################################################################################  
#------------------------#
# define water structure # (SPC)
#------------------------#
SPC_OW  = np.array([0.0     , 0.0    , 0.0])
SPC_HW1 = np.array([0.1     , 0.0    , 0.0]) # d(O-H) = 0.1 nm
SPC_HW2 = np.array([-0.03333, 0.09428, 0.0]) # a(H-O-H) = 109.47 degree

#------------------------#
# define water structure # (TIP4P/ice)
#------------------------#
# the distance between OW and MW should be 0.1577 Angstrom, I checked manually that this is correct!
TIP4Pice_OW  = np.array([0.0            , 0.0           , 0.0])
TIP4Pice_HW1 = np.array([0.09572        , 0.0           , 0.0]) # d(O-H) = 0.9572 nm
TIP4Pice_HW2 = np.array([-0.02399872084 , 0.09266272064 , 0.0]) # a(H-O-H) = 104.52 degree
TIP4Pice_MW  = np.array([0.009652490016 , 0.01247085936 , 0.0]) # r = r1 + a*(r2-r1) + b*(r3-r1)
                                                             #   with a = b = 0.13458335

# random seed
random.seed()

##########################################################################################
#                            1 x LATTICE                                                 # 
##########################################################################################  
if lattice_type == CUSTOM and len(filename_custom_lattice) == 1:
  file_lattice = open(filename_custom_lattice[0], "r")
  
  flag_read_latticepoints   = False
  flag_read_fixedpoints     = False
  flag_read_lattice_center  = False
  flag_read_lattice_rot     = False 
  flag_target_positions     = False
  flag_replicas             = False
  flag_lattice_displacement = False   # displacement of lattice as a whole
  flag_extra_displacement   = False   # displacement around individual lattice points
  flag_read_scale_lattice   = False
  flag_read_all_params      = False
  
  line_number             = 0
  line_number_absolute    = 0
  custom_n_latticepoints  = -1
  custom_n_fixedpoints    = -1
  
  line_end_latticepoints = 0
  line_end_fixed_points  = 0
  
  list_latticepoints         = []
  list_fixedpoints           = []
  lattice_center             = []
  list_target_positions      = []
  lattice_multiplier_x       = 1.0
  lattice_multiplier_y       = 1.0
  angle_lattice_rot          = 0.0
  n_replicas                 = 1
  displacement_xy            = 0.0 # same displacement for x and y
  displacement_x             = 0.0 # different displacement for x and y
  displacement_y             = 0.0 # different displacement for x and y
  displacement_z             = 0.0
  extra_displacement_cart_xy = 0.0
  extra_displacement_cart_z  = 0.0
  
  for line in file_lattice:
    line_number_absolute += 1
  
    # check if line is blank
    if len(line) == 1: # only contains newline
      continue
  
    # check if line continues with # --> comment
    if line.split()[0][0] == "#":
      continue
  
    line_number += 1
  
    # --- N_LATTICEPOINTS ---
    if line_number == 1:
      custom_n_latticepoints = int(line.split()[0])
  
      line_end_latticepoints = line_number + custom_n_latticepoints
      flag_read_latticepoints = True
      continue
    
    # --- lattice point n ---
    elif flag_read_latticepoints == True:
      list_line = line.split()
      index = int(list_line[0])
      x = float(list_line[1])
      y = float(list_line[2])
      z = float(list_line[3])
  
      latticepoint = [ index, x, y, z ]
      list_latticepoints.append(latticepoint)
  
      # was it the last lattice point?
      if line_end_latticepoints == line_number:
        flag_read_latticepoints = False
        flag_read_n_fixed_points  = True
      continue
    
    # --- N_FIXEDPOINTS ---
    elif flag_read_n_fixed_points == True:
      custom_n_fixedpoints = int(line.split()[0])
  
      if custom_n_fixedpoints > 0:
        line_end_fixed_points = line_number + custom_n_fixedpoints
        flag_read_n_fixed_points = False
        flag_read_fixedpoints = True
      else:
        flag_read_n_fixed_points = False
        flag_read_lattice_center =  True
  
      continue

    # --- fixed points ---
    # they have an assigned atom name and are not available on the lattice
    elif flag_read_fixedpoints == True:
      list_line = line.split()
      index = int(list_line[0])
      atomname = list_line[1]
  
      fixed_point = [ index, atomname ]
  
      list_fixedpoints.append(fixed_point)
  
      if line_number == line_end_fixed_points:
        flag_read_fixedpoints = False
        flag_read_lattice_center = True
      continue
  
    # --- read lattice center ---
    elif flag_read_lattice_center == True:
      list_line = line.split()
      if len(list_line) != 3:
        print "[ERROR]. The lattice_center is specified with the wrong numbers of parameters"
        print "         Expected: 3, Found: ", len(list_line), "..."
        print "         Please check line number", line_number_absolute
        sys.exit()
  
      x = float(list_line[0])
      y = float(list_line[1])
      z = float(list_line[2])
  
      lattice_center = [ x, y, z ]
  
      flag_read_lattice_center = False
      flag_replicas            = True
      continue
  
    # --- replicas ---
    elif flag_replicas == True:

      if line.split()[0] == "EXPLICIT":
        replica_type = EXPLICIT           
        n_replicas = int(line.split()[1])
      elif line.split()[0] == "COMPUTE":
        replica_type = COMPUTE            
        lattice_multiplier_x = int(line.split()[1])
        lattice_multiplier_y = int(line.split()[2])
        n_replicas = lattice_multiplier_x * lattice_multiplier_y
        lattice_UC_1         = [ float(line.split()[3]),  float(line.split()[4]),  float(line.split()[5])  ]
        lattice_UC_2         = [ float(line.split()[6]),  float(line.split()[7]),  float(line.split()[8])  ]
        lattice_UC_3         = [ float(line.split()[9]), float(line.split()[10]), float(line.split()[11]) ]
      else:
        print "[ERROR]. The keyword ", line.split[1], "is not known."
        print "         Please specify either 'EXPLICIT' or 'COMPUTE'"
        sys.exit()
  
      flag_replicas         = False
      flag_target_positions = True
      continue
  
    # --- target positions ---
    # these are the target locations for the centers of the lattice center
    # in EXPLICIT replica mode: for each replica one target center has to be specified
    #    the target position needs to be specified in cartesian coordinates (as if it is in slab)
    # in COMPUTE mode: multiply_x multiply_y xx xy yx yy zz
    #    only one z component is enough, because it is a 1D film in xy
    elif flag_target_positions == True:
      list_line = line.split()
      if len(list_line) != 3:
        print "[ERROR]. The target position is specified with the wrong number of parameters"
        print "         Expected: 3 (EXPLICIT mode), Found: ", len(list_line), "..."
        print "         Please check line number", line_number_absolute
        sys.exit()

      x = float(list_line[0])
      y = float(list_line[1])
      z = float(list_line[2])
  
      target = [ x/10.0, y/10.0, z/10.0 ]
      
      list_target_positions.append(target)
  
      flag_target_positions  = False
      flag_read_lattice_rot  = True
      continue

    # --- rotate custom lattice around z ---
    elif flag_read_lattice_rot == True:
      angle_lattice_rot = float(line.split()[0])
      flag_read_lattice_rot = False
  
      flag_lattice_displacement = True
      continue
  
    # --- displacement of whole lattice (in Angstrom) ---
    elif flag_lattice_displacement == True:
      # depending on number of arguments in this line,
      # decide on whether or not xy displacement is the same for x and y or different

      # xy the same
      if len(line.split()) == 2:
        displacement_xy = float(line.split()[0])
        displacement_z = float(line.split()[1])
  
        displacement_xy = displacement_xy/10.0 # convert Angstrom -> nm
        displacement_z = displacement_z/10.0 # convert Angstrom -> nm

        displacement_x = displacement_xy
        displacement_y = displacement_xy
      
      # xy different
      elif len(line.split()) == 3:
        displacement_x = float(line.split()[0])
        displacement_y = float(line.split()[1])
        displacement_z = float(line.split()[2])
  
        displacement_x = displacement_x/10.0 # convert Angstrom -> nm
        displacement_y = displacement_y/10.0 # convert Angstrom -> nm
        displacement_z = displacement_z/10.0 # convert Angstrom -> nm

      else:
        print "[ERROR]. The number of arguments for FLAG_LATTICE_DISPLACEMENT is wrong"
        print "         There should be either 2 or 3 arguments..."
        print line
        sys.exit()

      flag_lattice_displacement = False
      flag_extra_displacement = True
      continue
    
    # --- displacement around individual lattice points
    elif flag_extra_displacement == True:
      extra_displacement_cart_xy = float(line.split()[0])
      extra_displacement_cart_z = float(line.split()[1])
  
      extra_displacement_cart_xy = extra_displacement_cart_xy /10.0 # convert Angstrom -> nm
      extra_displacement_cart_z  = extra_displacement_cart_z  /10.0 # convert Angstrom -> nm
      
      flag_extra_displacement = False
      flag_read_scale_lattice = True
      continue

    # --- scale lattice to slab UC ---
    if flag_read_scale_lattice == True:
      if line.split()[0].lower() == "true":
        flag_scale_lattice = True
      elif line.split()[0].lower() == "false":
        flag_scale_lattice = False
      else:
        print "[ERROR]. flag_scale_lattice can only be 'true' or 'false'"
        sys.exit()

      flag_read_scale_lattice = False
      flag_read_all_params = True
  file_lattice.close()
  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            check if everything went well                              #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  if flag_read_all_params == False:
    print "[ERROR]. Not all parameters were specified in lattice file"
    sys.exit()
  
  if len(list_latticepoints) != custom_n_latticepoints:
    print "[ERROR]. There is a mistake in the number of lattice-points"
    print "         n_latticepoints:", custom_n_latticepoints
    print list_latticepoints
    sys.exit()
  
  if len(list_fixedpoints) != custom_n_fixedpoints:
    print "[ERROR]. There is a mistake in the number of fixed-points"
    print "         n_fixedpoints:", custom_n_fixedpoints
    print list_fixedpoints
    sys.exit()
  
  list_latticepoints_occupy = []
  for latticepoint in list_latticepoints:
    flag_fixedpoint = False
    for fixedpoint in list_fixedpoints:
      # check if they have the same index
      if latticepoint[0] == fixedpoint[0]:
        flag_fixedpoint = True
    if flag_fixedpoint == False:
      list_latticepoints_occupy.append(latticepoint)
  
  if len(list_latticepoints_occupy) + len(list_fixedpoints) != custom_n_latticepoints:
    print "[ERROR]. Could not split up latticepoints (should not happen)"
    sys.exit()
  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            build lattice                                              #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  lattice_custom = cl.lattice_building_block()
  
  # *** if multiplying clusters ***
  if replica_type == EXPLICIT:
    # unit cell of lattice equal to unit cell of slab
    lattice_custom.set_a(lattice_v_1[0], lattice_v_1[1], lattice_v_1[2])
    lattice_custom.set_b(lattice_v_2[0], lattice_v_2[1], lattice_v_2[2])
    lattice_custom.set_c(lattice_v_3[0], lattice_v_3[1], lattice_v_3[2])
  # *** if periodic lattice *** (eg basal Ih phase)
  elif replica_type == COMPUTE:
    lattice_custom.set_a(lattice_UC_1[0], lattice_UC_1[1], lattice_UC_1[2])
    lattice_custom.set_b(lattice_UC_2[0], lattice_UC_2[1], lattice_UC_2[2])
    lattice_custom.set_c(lattice_UC_3[0], lattice_UC_3[1], lattice_UC_3[2])
  
  for latticepoint in list_latticepoints_occupy:
    # move lattice-center to (0 0 0)
    move_to_target_x = -lattice_center[0]
    move_to_target_y = -lattice_center[1]
    move_to_target_z = -lattice_center[2]
  
    frac_coord = tc.cart2frac (lattice_UC_1, lattice_UC_2, lattice_UC_3, [ latticepoint[1]+move_to_target_x, \
                                                                           latticepoint[2]+move_to_target_y, \
                                                                           latticepoint[3]+move_to_target_z ])
    lattice_custom.add_latticepoint(frac_coord[0], frac_coord[1], frac_coord[2])
  
  # build lattice in lattice coordinates
  cl.construct_lattice(lattice_custom, lattice_multiplier_x, lattice_multiplier_y)
  
  # the lattice vectors have to come in AA (new requirement for construct_lattice.py)
  # --> to make it compatible with eg LAMMPS RSS
  # right now they are specified in nm (from *gro file)
  lattice_v_1_AA = [ i*10.0 for i in lattice_v_1]
  lattice_v_2_AA = [ i*10.0 for i in lattice_v_2]
  lattice_v_3_AA = [ i*10.0 for i in lattice_v_3]
  
  lattice_custom.lattice2slab(lattice_v_1_AA, lattice_v_2_AA, lattice_v_3_AA, flag_scale_lattice)


#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
#                            define 2 x custom lattice                                  #  
#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
if lattice_type == CUSTOM and len(filename_custom_lattice) == 2:
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #                            read lattice 1                                             #  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  file_lattice_1 = open(filename_custom_lattice[0], "r")
  
  flag_read_latticepoints  = False
  flag_read_fixedpoints    = False
  flag_read_lattice_center = False
  flag_read_lattice_rot    = False
  flag_target_positions    = False
  flag_replicas            = False
  flag_lattice_displacement= False
  flag_extra_displacement  = False
  flag_read_scale_lattice  = False
  flag_read_all_params     = False
  
  line_number             = 0
  line_number_absolute    = 0
  custom_n_latticepoints  = [-1, -1]
  custom_n_fixedpoints    = [-1, -1]
  
  line_end_latticepoints = [0, 0]
  line_end_fixed_points  = [0, 0]
  
  list_latticepoints_1    = []
  list_latticepoints_2    = []
  list_fixedpoints_1      = []
  list_fixedpoints_2      = []
  lattice_center_1        = []
  lattice_center_2        = []
  list_target_positions_1 = []
  list_target_positions_2 = []
  angle_lattice_rot       = [0.0, 0.0]
  n_replicas              = [1, 1]
  replica_type            = [COMPUTE, COMPUTE]
  lattice_multiplier_x_1  = 0
  lattice_multiplier_y_1  = 0
  displacement_xy         = [0.0, 0.0]
  displacement_x          = [0.0, 0.0]
  displacement_y          = [0.0, 0.0]
  displacement_z          = [0.0, 0.0]
  extra_displacement_cart_xy = [0.0, 0.0]
  extra_displacement_cart_z  = [0.0, 0.0]
  flag_scale_lattice      = [ False, False ]
  
  
  for line in file_lattice_1:
    line_number_absolute += 1
  
    # check if line is blank
    if len(line) == 1: # only contains newline
      continue
  
    # check if line continues with # --> comment
    if line.split()[0][0] == "#":
      continue
  
    line_number += 1
  
    # --- N_LATTICEPOINTS ---
    if line_number == 1:
      custom_n_latticepoints[0] = int(line.split()[0])
  
      line_end_latticepoints[0] = line_number + custom_n_latticepoints[0]
      flag_read_latticepoints = True
      continue
    
    # --- lattice point n ---
    elif flag_read_latticepoints == True:
      list_line = line.split()
      index = int(list_line[0])
      x = float(list_line[1])
      y = float(list_line[2])
      z = float(list_line[3])
  
      latticepoint = [ index, x, y, z ]
      list_latticepoints_1.append(latticepoint)
  
      # was it the last lattice point?
      if line_end_latticepoints[0] == line_number:
        flag_read_latticepoints = False
        flag_read_n_fixed_points  = True
      continue
    
    # --- N_FIXEDPOINTS ---
    elif flag_read_n_fixed_points == True:
      custom_n_fixedpoints[0] = int(line.split()[0])
  
      if custom_n_fixedpoints[0] > 0:
        line_end_fixed_points = line_number + custom_n_fixedpoints[0]
        flag_read_n_fixed_points = False
        flag_read_fixedpoints = True
      else:
        flag_read_n_fixed_points = False
        flag_read_lattice_center =  True
  
      continue
  
    # --- fixed points ---
    # they have an assigned atom name and are not available on the lattice
    elif flag_read_fixedpoints == True:
      list_line = line.split()
      index = int(list_line[0])
      atomname = list_line[1]
  
      fixed_point = [ index, atomname ]
  
      list_fixedpoints_1.append(fixed_point)
  
      if line_number == line_end_fixed_points:
        flag_read_fixedpoints = False
        flag_read_lattice_center = True
      continue
  
    # --- read lattice center ---
    elif flag_read_lattice_center == True:
      list_line = line.split()
      if len(list_line) != 3:
        print "[ERROR]. The lattice_center is specified with the wrong numbers of parameters"
        print "         Expected: 3, Found: ", len(list_line), "..."
        print "         Please check line number", line_number_absolute
        sys.exit()
  
      x = float(list_line[0])
      y = float(list_line[1])
      z = float(list_line[2])
  
      lattice_center_1 = [ x, y, z ]
  
      flag_read_lattice_center = False
      flag_replicas            = True
      continue
  
    # --- replicas ---
    elif flag_replicas == True:
      if line.split()[0] == "EXPLICIT":
        replica_type[0] = EXPLICIT           
        n_replicas[0] = int(line.split()[1])
      elif line.split()[0] == "COMPUTE":
        replica_type[0] = COMPUTE            
        lattice_multiplier_x_1 = int(line.split()[1])
        lattice_multiplier_y_1 = int(line.split()[2])
        n_replicas[0] = lattice_multiplier_x_1 * lattice_multiplier_y_1
        lattice_UC_1_1       = [ float(line.split()[3]),  float(line.split()[4]),  float(line.split()[5])  ]
        lattice_UC_2_1       = [ float(line.split()[6]),  float(line.split()[7]),  float(line.split()[8])  ]
        lattice_UC_3_1       = [ float(line.split()[9]), float(line.split()[10]), float(line.split()[11]) ]
      else:
        print "[ERROR]. The keyword ", line.split[1], "is not known."
        print "         Please specify either 'EXPLICIT' or 'COMPUTE'"
        sys.exit()
  
      flag_replicas         = False
      flag_target_positions = True
      continue
  
    # --- target positionss ---
    # these are the target locations for the centers of the lattice center
    # for each replica one target center has to be specified
    # the target position needs to be specified in cartesian coordinates (as if it is in slab)
    elif flag_target_positions == True:
      list_line = line.split()
      if len(list_line) != 3:
        print "[ERROR]. The target position is spedified with the wrong number of parameters"
        print "         Expected: 3, Found: ", len(list_line), "..."
        print "         Please check line number", line_number_absolute
        sys.exit()
  
      x = float(list_line[0])
      y = float(list_line[1])
      z = float(list_line[2])
  
      target = [ x/10.0, y/10.0, z/10.0 ]
      
      list_target_positions_1.extend(target)
  
      flag_target_positions  = False
      flag_read_lattice_rot  = True
      continue
    
    # --- rotate custom lattice around z ---
    elif flag_read_lattice_rot == True:
      angle_lattice_rot[0] = float(line.split()[0])
      flag_read_lattice_rot = False
  
      flag_lattice_displacement = True
      continue
  
    # --- displacement of whole lattice (in Angstrom) ---
    elif flag_lattice_displacement == True:
      # xy the same
      if len(line.split()) == 2:
        displacement_xy[0] = float(line.split()[0])
        displacement_z[0] = float(line.split()[1])
  
        displacement_xy[0] = displacement_xy[0]/10.0 # convert Angstrom -> nm
        displacement_z[0] = displacement_z[0]/10.0 # convert Angstrom -> nm
        
        displacement_x[0] = displacement_xy[0]
        displacement_y[0] = displacement_xy[0]

        flag_lattice_displacement = False
        flag_extra_displacement = True
        continue
      
      # xy different
      elif len(line.split()) == 3:
        displacement_x[0] = float(line.split()[0])
        displacement_y[0] = float(line.split()[1])
        displacement_z[0] = float(line.split()[2])
  
        displacement_x[0] = displacement_x[0]/10.0 # convert Angstrom -> nm
        displacement_y[0] = displacement_y[0]/10.0 # convert Angstrom -> nm
        displacement_z[0] = displacement_z[0]/10.0 # convert Angstrom -> nm
        
        flag_lattice_displacement = False
        flag_extra_displacement = True
        continue
    
      else:
        print "[ERROR]. The number of arguments for FLAG_LATTICE_DISPLACEMENT is wrong"
        print "         There should be either 2 or 3 arguments..."
        print line
        sys.exit()

    # --- displacement around individual lattice points
    elif flag_extra_displacement == True:
      extra_displacement_cart_xy[0] = float(line.split()[0])
      extra_displacement_cart_z[0] = float(line.split()[1])
  
      extra_displacement_cart_xy[0] = extra_displacement_cart_xy[0] /10.0 # convert Angstrom -> nm
      extra_displacement_cart_z[0]  = extra_displacement_cart_z[0]  /10.0 # convert Angstrom -> nm
      
      flag_extra_displacement = False
      flag_read_scale_lattice = True
      continue

    # --- scale lattice to slab UC ---
    if flag_read_scale_lattice == True:
      if line.split()[0].lower() == "true":
        flag_scale_lattice[0] = True
      elif line.split()[0].lower() == "false":
        flag_scale_lattice[0] = False
      else:
        print "[ERROR]. flag_scale_lattice can only be 'true' or 'false'"
        sys.exit()

      flag_read_scale_lattice = False
      flag_read_all_params = True
  file_lattice_1.close()
  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            check if everything went well (lattice 1)                  #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  if flag_read_all_params == False:
    print "[ERROR]. Not all parameters were specified in lattice file"
    sys.exit()
  
  if len(list_latticepoints_1) != custom_n_latticepoints[0]:
    print "[ERROR]. There is a mistake in the number of lattice-points"
    print "         n_latticepoints:", custom_n_latticepoints[0]
    print list_latticepoints_1
    sys.exit()
  
  if len(list_fixedpoints_1) != custom_n_fixedpoints[0]:
    print "[ERROR]. There is a mistake in the number of fixed-points"
    print "         n_fixedpoints:", custom_n_fixedpoints[0]
    print list_fixedpoints_1
    sys.exit()
  
  list_latticepoints_occupy = []
  for latticepoint in list_latticepoints_1:
    flag_fixedpoint = False
    for fixedpoint in list_fixedpoints_1:
      # check if they have the same index
      if latticepoint[0] == fixedpoint[0]:
        flag_fixedpoint = True
    if flag_fixedpoint == False:
      list_latticepoints_occupy.append(latticepoint)
  
  if len(list_latticepoints_occupy) + len(list_fixedpoints_1) != custom_n_latticepoints[0]:
    print "[ERROR]. Could not split up latticepoints in lattice 1 (should not happen)"
    sys.exit()
  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            build lattice 1                                            #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  lattice_custom_1 = cl.lattice_building_block()
  
  # *** if multiplying clusters ***
  if replica_type[0] == EXPLICIT:
    # unit cell of lattice equal to unit cell of slab
    lattice_custom_1.set_a(lattice_v_1[0], lattice_v_1[1], lattice_v_1[2])
    lattice_custom_1.set_b(lattice_v_2[0], lattice_v_2[1], lattice_v_2[2])
    lattice_custom_1.set_c(lattice_v_3[0], lattice_v_3[1], lattice_v_3[2])
  # *** if periodic lattice *** (eg basal Ih phase)
  elif replica_type[0] == COMPUTE:
    lattice_custom_1.set_a(lattice_UC_1_1[0], lattice_UC_1_1[1], lattice_UC_1_1[2])
    lattice_custom_1.set_b(lattice_UC_2_1[0], lattice_UC_2_1[1], lattice_UC_2_1[2])
    lattice_custom_1.set_c(lattice_UC_3_1[0], lattice_UC_3_1[1], lattice_UC_3_1[2])
  
  for latticepoint in list_latticepoints_occupy:
    # move lattice-center to (0 0 0)
    move_to_target_x = -lattice_center_1[0]
    move_to_target_y = -lattice_center_1[1]
    move_to_target_z = -lattice_center_1[2]
  
    frac_coord = tc.cart2frac (lattice_UC_1_1, lattice_UC_2_1, lattice_UC_3_1, [ latticepoint[1]+move_to_target_x, \
                                                                                 latticepoint[2]+move_to_target_y, \
                                                                                 latticepoint[3]+move_to_target_z ])
    lattice_custom_1.add_latticepoint(frac_coord[0], frac_coord[1], frac_coord[2])
  
  # build lattice in lattice coordinates
  cl.construct_lattice(lattice_custom_1, lattice_multiplier_x_1, lattice_multiplier_y_1)
  
  # the lattice vectors have to come in AA (new requirement for construct_lattice.py)
  # --> to make it compatible with eg LAMMPS RSS
  # right now they are specified in nm (from *gro file)
  lattice_v_1_AA = [ i*10.0 for i in lattice_v_1]
  lattice_v_2_AA = [ i*10.0 for i in lattice_v_2]
  lattice_v_3_AA = [ i*10.0 for i in lattice_v_3]
  lattice_custom_1.lattice2slab(lattice_v_1_AA, lattice_v_2_AA, lattice_v_3_AA, flag_scale_lattice[0])
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #                            read lattice 2                                             #  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  file_lattice_2 = open(filename_custom_lattice[1], "r")

  # reset all the flags
  flag_read_latticepoints  = False
  flag_read_fixedpoints    = False
  flag_read_lattice_center = False
  flag_read_lattice_rot    = False
  flag_target_positions    = False
  flag_replicas            = False
  flag_lattice_displacement= False
  flag_extra_displacement  = False
  flag_read_scale_lattice  = False
  flag_read_all_params     = False
  
  line_number             = 0
  line_number_absolute    = 0
  
  
  for line in file_lattice_2:
    line_number_absolute += 1
  
    # check if line is blank
    if len(line) == 1: # only contains newline
      continue
  
    # check if line continues with # --> comment
    if line.split()[0][0] == "#":
      continue
  
    line_number += 1
  
    # --- N_LATTICEPOINTS ---
    if line_number == 1:
      custom_n_latticepoints[1] = int(line.split()[0])
  
      line_end_latticepoints[1] = line_number + custom_n_latticepoints[1]
      flag_read_latticepoints = True
      continue
    
    # --- lattice point n ---
    elif flag_read_latticepoints == True:
      list_line = line.split()
      index = int(list_line[0])
      x = float(list_line[1])
      y = float(list_line[2])
      z = float(list_line[3])
  
      latticepoint = [ index, x, y, z ]
      list_latticepoints_2.append(latticepoint)
  
      # was it the last lattice point?
      if line_end_latticepoints[1] == line_number:
        flag_read_latticepoints = False
        flag_read_n_fixed_points  = True
      continue
    
    # --- N_FIXEDPOINTS ---
    elif flag_read_n_fixed_points == True:
      custom_n_fixedpoints[1] = int(line.split()[0])
  
      if custom_n_fixedpoints[1] > 0:
        line_end_fixed_points = line_number + custom_n_fixedpoints[1]
        flag_read_n_fixed_points = False
        flag_read_fixedpoints = True
      else:
        flag_read_n_fixed_points = False
        flag_read_lattice_center =  True
  
      continue
  
    # --- fixed points ---
    # they have an assigned atom name and are not available on the lattice
    elif flag_read_fixedpoints == True:
      list_line = line.split()
      index = int(list_line[0])
      atomname = list_line[1]
  
      fixed_point = [ index, atomname ]
  
      list_fixedpoints_2.append(fixed_point)
  
      if line_number == line_end_fixed_points:
        flag_read_fixedpoints = False
        flag_read_lattice_center = True
      continue
  
    # --- read lattice center ---
    elif flag_read_lattice_center == True:
      list_line = line.split()
      if len(list_line) != 3:
        print "[ERROR]. The lattice_center is specified with the wrong numbers of parameters"
        print "         Expected: 3, Found: ", len(list_line), "..."
        print "         Please check line number", line_number_absolute
        sys.exit()
  
      x = float(list_line[0])
      y = float(list_line[1])
      z = float(list_line[2])
  
      lattice_center_2 = [ x, y, z ]
  
      flag_read_lattice_center = False
      flag_replicas            = True
      continue
  
    # --- replicas ---
    elif flag_replicas == True:
      if line.split()[0] == "EXPLICIT":
        replica_type[1] = EXPLICIT           
        n_replicas[1] = int(line.split()[1])
      elif line.split()[0] == "COMPUTE":
        replica_type[1] = COMPUTE            
        lattice_multiplier_x_2 = int(line.split()[1])
        lattice_multiplier_y_2 = int(line.split()[2])
        n_replicas[1] = lattice_multiplier_x_2 * lattice_multiplier_y_2
        lattice_UC_1_2       = [ float(line.split()[3]),  float(line.split()[4]),  float(line.split()[5])  ]
        lattice_UC_2_2       = [ float(line.split()[6]),  float(line.split()[7]),  float(line.split()[8])  ]
        lattice_UC_3_2       = [ float(line.split()[9]), float(line.split()[10]), float(line.split()[11]) ]
      else:
        print "[ERROR]. The keyword ", line.split[1], "is not known."
        print "         Please specify either 'EXPLICIT' or 'COMPUTE'"
        sys.exit()
  
      flag_replicas         = False
      flag_target_positions = True
      continue
  
    # --- target positionss ---
    # these are the target locations for the centers of the lattice center
    # for each replica one target center has to be specified
    # the target position needs to be specified in cartesian coordinates (as if it is in slab)
    elif flag_target_positions == True:
      list_line = line.split()
      if len(list_line) != 3:
        print "[ERROR]. The target position is spedified with the wrong number of parameters"
        print "         Expected: 3, Found: ", len(list_line), "..."
        print "         Please check line number", line_number_absolute
        sys.exit()
  
      x = float(list_line[0])
      y = float(list_line[1])
      z = float(list_line[2])
  
      target = [ x/10.0, y/10.0, z/10.0 ]
      

      list_target_positions_2.extend(target)
  
      flag_target_positions  = False
      flag_read_lattice_rot  = True
      continue
    
    # --- rotate custom lattice around z ---
    elif flag_read_lattice_rot == True:
      angle_lattice_rot[1] = float(line.split()[0])
      flag_read_lattice_rot = False
  
      flag_lattice_displacement = True
      continue
  
    # --- displacement of whole lattice (in Angstrom) ---
    elif flag_lattice_displacement == True:
      # xy the same
      if len(line.split()) == 2:
        displacement_xy[1] = float(line.split()[0])
        displacement_z[1] = float(line.split()[1])
  
        displacement_xy[1] = displacement_xy[1]/10.0 # convert Angstrom -> nm
        displacement_z[1] = displacement_z[1]/10.0 # convert Angstrom -> nm
        
        displacement_x[1] = displacement_xy[1]
        displacement_y[1] = displacement_xy[1]

        flag_lattice_displacement = False
        flag_extra_displacement = True
        continue
      
      # xy different
      elif len(line.split()) == 3:
        displacement_x[1] = float(line.split()[0])
        displacement_y[1] = float(line.split()[1])
        displacement_z[1] = float(line.split()[2])
  
        displacement_x[1] = displacement_x[1]/10.0 # convert Angstrom -> nm
        displacement_y[1] = displacement_y[1]/10.0 # convert Angstrom -> nm
        displacement_z[1] = displacement_z[1]/10.0 # convert Angstrom -> nm
        
        flag_lattice_displacement = False
        flag_extra_displacement = True
        continue
    
      else:
        print "[ERROR]. The number of arguments for FLAG_LATTICE_DISPLACEMENT is wrong"
        print "         There should be either 2 or 3 arguments..."
        print line
        sys.exit()
    
    # --- displacement around individual lattice points
    elif flag_extra_displacement == True:
      extra_displacement_cart_xy[1] = float(line.split()[0])
      extra_displacement_cart_z[1] = float(line.split()[1])
  
      extra_displacement_cart_xy[1] = extra_displacement_cart_xy[1] /10.0 # convert Angstrom -> nm
      extra_displacement_cart_z[1]  = extra_displacement_cart_z[1]  /10.0 # convert Angstrom -> nm
      
      flag_extra_displacement = False
      flag_read_scale_lattice = True
      continue

    # --- scale lattice to slab UC ---
    if flag_read_scale_lattice == True:
      if line.split()[0].lower() == "true":
        flag_scale_lattice[1] = True
      elif line.split()[0].lower() == "false":
        flag_scale_lattice[1] = False
      else:
        print "[ERROR]. flag_scale_lattice can only be 'true' or 'false'"
        sys.exit()

      flag_read_scale_lattice = False
      flag_read_all_params = True
      flag_read_all_params = True
  file_lattice_2.close()
  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            check if everything went well (lattice 2)                  #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  if flag_read_all_params == False:
    print "[ERROR]. Not all parameters were specified in lattice file"
    sys.exit()
  
  if len(list_latticepoints_2) != custom_n_latticepoints[1]:
    print "[ERROR]. There is a mistake in the number of lattice-points"
    print "         n_latticepoints:", custom_n_latticepoints[1]
    print list_latticepoints_2
    sys.exit()
  
  if len(list_fixedpoints_2) != custom_n_fixedpoints[1]:
    print "[ERROR]. There is a mistake in the number of fixed-points"
    print "         n_fixedpoints:", custom_n_fixedpoints[1]
    print list_fixedpoints_2
    sys.exit()
  
  list_latticepoints_occupy = []
  for latticepoint in list_latticepoints_2:
    flag_fixedpoint = False
    for fixedpoint in list_fixedpoints_2:
      # check if they have the same index
      if latticepoint[1] == fixedpoint[1]:
        flag_fixedpoint = True
    if flag_fixedpoint == False:
      list_latticepoints_occupy.append(latticepoint)
  
  if len(list_latticepoints_occupy) + len(list_fixedpoints_2) != custom_n_latticepoints[1]:
    print "[ERROR]. Could not split up latticepoints in lattice_2 (should not happen)"
    print len(list_latticepoints_occupy)
    print len(list_fixedpoints_2)
    print custom_n_latticepoints[1]
    sys.exit()
  

  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            rescale target position lattice 2                          #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  """ if lattice positions are scaled, then it should also be the case that the target position
       of lattice 2 is scaled appropriately...
       if eg. basal on top of basal, and the bottom+top basal lattice is scaled,
       but not the target position of the second basal --> there will be a mismatch between lattice 1 and lattice 2
       --> the target position will be scaled in the exact same way then latticepoints in cl.lattice2slab(...)
       don't forget to convert slab lattice from nm to AA
  """
  if flag_scale_lattice[1] == True:
    quotient_a = lattice_UC_1_2[0] / (lattice_v_1[0]*10.0) # scale x of a axis 
    list_target_positions_2[0] = list_target_positions_2[0]/quotient_a

    quotient_b = lattice_UC_2_2[1] / (lattice_v_2[1]*10.0) # scale y of b axis
    list_target_positions_2[1] = list_target_positions_2[1]/quotient_b

  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            build lattice 2                                            #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  lattice_custom_2 = cl.lattice_building_block()
  
  # *** if multiplying clusters ***
  if replica_type[1] == EXPLICIT:
    # unit cell of lattice equal to unit cell of slab
    lattice_custom_2.set_a(lattice_v_1[0], lattice_v_1[1], lattice_v_1[2])
    lattice_custom_2.set_b(lattice_v_2[0], lattice_v_2[1], lattice_v_2[2])
    lattice_custom_2.set_c(lattice_v_3[0], lattice_v_3[1], lattice_v_3[2])
  # *** if periodic lattice *** (eg basal Ih phase)
  elif replica_type[1] == COMPUTE:
    lattice_custom_2.set_a(lattice_UC_1_2[0], lattice_UC_1_2[1], lattice_UC_1_2[2])
    lattice_custom_2.set_b(lattice_UC_2_2[0], lattice_UC_2_2[1], lattice_UC_2_2[2])
    lattice_custom_2.set_c(lattice_UC_3_2[0], lattice_UC_3_2[1], lattice_UC_3_2[2])
  
  for latticepoint in list_latticepoints_occupy:
    # move lattice-center to (0 0 0)
    move_to_target_x = -lattice_center_1[0]
    move_to_target_y = -lattice_center_1[1]
    move_to_target_z = -lattice_center_1[2]
  

    frac_coord = tc.cart2frac (lattice_UC_1_2, lattice_UC_2_2, lattice_UC_3_2, [ latticepoint[1]+move_to_target_x, \
                                                                                 latticepoint[2]+move_to_target_y, \
                                                                                 latticepoint[3]+move_to_target_z ])
    lattice_custom_2.add_latticepoint(frac_coord[0], frac_coord[1], frac_coord[2])
  
  # build lattice in lattice coordinates
  cl.construct_lattice(lattice_custom_2, lattice_multiplier_x_2, lattice_multiplier_y_2)
  
  # the lattice vectors have to come in AA (new requirement for construct_lattice.py)
  # --> to make it compatible with eg LAMMPS RSS
  # right now they are specified in nm (from *gro file)
  lattice_v_1_AA = [ i*10.0 for i in lattice_v_1]
  lattice_v_2_AA = [ i*10.0 for i in lattice_v_2]
  lattice_v_3_AA = [ i*10.0 for i in lattice_v_3]
  lattice_custom_2.lattice2slab(lattice_v_1_AA, lattice_v_2_AA, lattice_v_3_AA, flag_scale_lattice[1])

##########################################################################################
#                            choose lattice to use                                       # 
##########################################################################################  
lattice   = []  # use if only one lattice is used
lattice_1 = []  # lattice 1/2
lattice_2 = []  # lattice 2/2

# --- custom lattice (from file) ---
if lattice_type == CUSTOM:
  if len(filename_custom_lattice) == 1:
    for latticepoint in lattice_custom.list_latticepoints:
      lattice.append(latticepoint.tolist())
  elif len(filename_custom_lattice) == 2:
    for latticepoint in lattice_custom_1.list_latticepoints:
      lattice_1.append(latticepoint.tolist())
    for latticepoint in lattice_custom_2.list_latticepoints:
      lattice_2.append(latticepoint.tolist())
  else:
    print "[ERROR]. Only two different lattices supported at the moment"
    sys.exit()

# --- undefined lattice
else:
  print "[ERROR]. lattice-type not defined", lattice_type
  sys.exit()

# choose way of positioning lattice
if lattice_type == CUSTOM:
  flag_lattice_location = CENTERED # center lattice to target
else:
  flag_lattice_location = RANDOM

# make sure there are not more water molecules than lattice-points
if len(filename_custom_lattice) == 1:
  if n_h2o > len(lattice):
    print "[ERROR]. There are more water molecules than lattice points"
    print "n_H2O:", n_h2o
    print "n_latticepoints:", len(lattice)
    sys.exit()

  # check if lattice is filled completely or not
  #  --> useful for the way lattice points are picked for occupation
  if n_h2o == len(lattice):
    flag_lattice_filling = FILLED_COMPLETELY
  elif n_h2o < len(lattice):
    flag_lattice_filling = FILLED_PARTIALLY

else:
  if n_h2o > (len(lattice_1) + len(lattice_2)):
    print "[ERROR]. There are more water molecules than lattice points"
    print "n_H2O:", n_h2o
    print "n_latticepoints:", len(lattice_1) + len(lattice_2)
    sys.exit()
  
  # check if lattice is filled completely or not
  #  --> useful for the way lattice points are picked for occupation
  if n_h2o == (len(lattice_1) + len(lattice_2)):
    flag_lattice_filling = FILLED_COMPLETELY
  elif n_h2o < (len(lattice_1) + len(lattice_2)):
    flag_lattice_filling = FILLED_PARTIALLY

# if equal occupation of lattice sites, check if correct input
if flag_lattice== True and flag_lattice_point_occupation == EQUAL_OCCUPATION:
  if n_structures % len(lattice) != 0:
    print "[ERROR]. For equal occupation of lattice-sites, n_structures/len(lattice) needs"
    print "         to be an integer"
    print "         remainder: ", n_structures % len(lattice)
    sys.exit()


##########################################################################################
#                            check if latticepoints are not too close                    # 
##########################################################################################  
"""
check if there are any latticepoints that are closer than some treshhold value
 --> if that is the case (eg some lattices that come from VESTA-crystallites): STOP
"""
lp_cart_tmp = []
for idx_1, lp_1 in enumerate(lattice):
  # in cartesian coords
  x1 = lp_1[0] * lattice_v_1[coord_X] + \
       lp_1[1] * lattice_v_2[coord_X] + \
       lp_1[2] * lattice_v_3[coord_X]

  y1 = lp_1[0] * lattice_v_1[coord_Y] + \
       lp_1[1] * lattice_v_2[coord_Y] + \
       lp_1[2] * lattice_v_3[coord_Y]

  z1 = lp_1[0] * lattice_v_1[coord_Z] + \
       lp_1[1] * lattice_v_2[coord_Z] + \
       lp_1[2] * lattice_v_3[coord_Z]

  lp_cart_tmp.append([x1, y1, z1])


for idx_1, lp_1 in enumerate(lp_cart_tmp):
  for idx_2, lp_2 in enumerate(lp_cart_tmp):
    if idx_1 <= idx_2:
      continue

    d = pbc_distance(lp_1, lp_2, lattice_v_1, lattice_v_2, lattice_v_3)

    if d < 0.2:
      print "[ERROR]. Two latticepoints are closer than the allowed treshhold!"
      print "         lp_1:", idx_1, lp_1
      print "         lp_2:", idx_2, lp_2
      print "         d   :", d
      print lattice_v_1, lattice_v_2, lattice_v_3
      sys.exit()

##########################################################################################
#                            generate n_structures water-structures                      # 
##########################################################################################  

try:
  n_structures_start
except:
  n_structures_start = 0

n = n_structures_start - 1 # start with n=0, messy because of backward compatibilty

while n < (n_structures_start + n_structures-1):
  n += 1
  if n%100 == 0 and n>0:
    print "\t(+) done:", n, "structures"

  # CORRECTION HERE?
  molecule_number = 2               # molecule 1: FELSP
  # molecule_number = 1               # debug: only H2O
  atom_index = len(list_slab_atoms) # atoms 1-n : FELSP
  
  list_water = []
  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  #                            generate water structure completely random                 #  
  #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
  if flag_lattice == False:
    for i in range(0, n_h2o):
      #--------------------#
      # decide orientation # (random: axis and angle)
      #--------------------#

      if water_model == SPC:
        # generate rotation axis (randomly in x,y and z)
        rot_axis = np.array([ random.uniform(-1.0,1.0) ,
                              random.uniform(-1.0,1.0) ,
                              random.uniform(-1.0,1.0) ])  
    
        # generate rotation angle
        rot_angle = random.uniform(0.0, 2*math.pi) # angle in radians
        
        oxygen_atom = [ molecule_number,
                        "SOL",
                        "OW",
                        atom_index+1,
                        np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[0],
                        np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[1],
                        np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[2]  ]
        
        hydrogen_atom_1 = [ molecule_number,
                            "SOL",
                            "HW1",
                            atom_index+2,
                            np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[0],
                            np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[1],
                            np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[2]  ]
        
        hydrogen_atom_2 = [ molecule_number,
                            "SOL",
                            "HW2",
                            atom_index+3,
                            np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[0],
                            np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[1],
                            np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[2]  ]
    
        atom_number     += 3
        molecule_number += 1
        
        #-----------------#
        # decide position # (in fractional coords)
        #-----------------#
        position_frac = [0.0, 0.0, 0.0]

        position_frac = [ random.uniform(0.0, 1.0) ,
                          random.uniform(0.0, 1.0) ]
     
        height_z = 2.0 # in nm (because of GROMACS)
        height_displacement = random.uniform(-0.10, 0.10) # in Angstrom
        
        # calculate displacement in cartesian coords
        displacement_cart_x = position_frac[coord_X] * lattice_v_1[coord_X] + \
                              position_frac[coord_Y] * lattice_v_2[coord_X]
        
        displacement_cart_y = position_frac[coord_X] * lattice_v_1[coord_Y] + \
                              position_frac[coord_Y] * lattice_v_2[coord_Y]
        
        displacement_cart_z = height_z + height_displacement

    
        #-----------------------------------#
        # move atom to its location on slab #
        #-----------------------------------#
        oxygen_atom[4] += displacement_cart_x
        oxygen_atom[5] += displacement_cart_y
        oxygen_atom[6] += displacement_cart_z
    
        hydrogen_atom_1[4] += displacement_cart_x
        hydrogen_atom_1[5] += displacement_cart_y
        hydrogen_atom_1[6] += displacement_cart_z
        
        hydrogen_atom_2[4] += displacement_cart_x
        hydrogen_atom_2[5] += displacement_cart_y
        hydrogen_atom_2[6] += displacement_cart_z
        
        #--------------------------#
        # generate water structure #
        #--------------------------#
        # list_water.append([oxygen_atom, hydrogen_atom_1, hydrogen_atom_2])
        list_water.append(oxygen_atom)
        list_water.append(hydrogen_atom_1)
        list_water.append(hydrogen_atom_2)
      else:
        print "[ERROR]. Only SPC water model supported currently"
        sys.exit()
  
  #######################################################################
  #                 generate water structure with lattice               #
  #######################################################################
  if flag_lattice == True:
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    #                            only one lattice to be filled                              #  
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if len(filename_custom_lattice) == 1:
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            decide orientation in xy plane                             #
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      # rotate around z (random value between 0 and angle_lattice_rot)
      # the lattice center was already moved to (0 0 0)

      # generate rotation angle (between 0.0 and angle_lattice_rot degree)
      rot_angle = random.uniform(0.0, angle_lattice_rot*math.pi/180.0) # angle in radians
      lattice_rotation_angle = rot_angle

      # deepcopy lattice, so it starts from same orientation in each run!
      # all lattice operations are done on lattice_tmp, lattice stores the original position
      lattice_tmp = copy.deepcopy(lattice)

      for latticepoint in lattice_tmp:
        # rotation axis --> z
        rot_axis = np.array([ 0.0 ,
                              0.0 ,
                              1.0 ])  

        # in cartesian coords
        x = latticepoint[0] * lattice_v_1[coord_X] + \
            latticepoint[1] * lattice_v_2[coord_X] + \
            latticepoint[2] * lattice_v_3[coord_X]
        
        y = latticepoint[0] * lattice_v_1[coord_Y] + \
            latticepoint[1] * lattice_v_2[coord_Y] + \
            latticepoint[2] * lattice_v_3[coord_Y]

        z = latticepoint[0] * lattice_v_1[coord_Z] + \
            latticepoint[1] * lattice_v_2[coord_Z] + \
            latticepoint[2] * lattice_v_3[coord_Z]
        
        before_rot = [ x, y, z ]
        after_rot = np.dot(rot.rotation_matrix(rot_axis, rot_angle), before_rot) 
        after_rot_frac = tc.cart2frac(lattice_v_1, lattice_v_2, lattice_v_3, after_rot)

        latticepoint[0] = after_rot_frac[0]
        latticepoint[1] = after_rot_frac[1]
        latticepoint[2] = after_rot_frac[2]
      
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            decide position of lattice                                 #  
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      displacement_lattice_cart_x = 0
      displacement_lattice_cart_y = 0
      displacement_lattice_cart_z = 0

      #--------------------#
      # RANDOM POSITIONING #
      #--------------------#
      if flag_lattice_location == RANDOM:
        position_frac = [ random.uniform(0.0, 1.0) ,
                          random.uniform(0.0, 1.0) ,
                          random.uniform(0.0, 1.0) ]
        
        # calculate displacement in cartesian coords
        displacement_lattice_cart_x = position_frac[coord_X] * lattice_v_1[coord_X] + \
                                      position_frac[coord_Y] * lattice_v_2[coord_X] + \
                                      position_frac[coord_Z] * lattice_v_3[coord_X]
        
        displacement_lattice_cart_y = position_frac[coord_X] * lattice_v_1[coord_Y] + \
                                      position_frac[coord_Y] * lattice_v_2[coord_Y] + \
                                      position_frac[coord_Z] * lattice_v_3[coord_Y]
        
        displacement_lattice_cart_z = position_frac[coord_X] * lattice_v_1[coord_Z] + \
                                      position_frac[coord_Y] * lattice_v_2[coord_Z] + \
                                      position_frac[coord_Z] * lattice_v_3[coord_Z]
      #----------------------#
      # USER DEFINED LATTICE #
      #----------------------#
      elif flag_lattice_location == CENTERED:
        # the user has to specify the center of the lattice as well as the target location of the lattice
        # ==> with this info the lattice will be placed onto the target location with some random displacement

        # --- CAREFUL, right now only one cluster replica is supported ---
        # the lattice center right now is at (0 0 0) [ from previous rotation ]
        displacement_lattice_cart_x = target[0] + random.uniform(-displacement_x, +displacement_x)
        displacement_lattice_cart_y = target[1] + random.uniform(-displacement_y, +displacement_y)
        displacement_lattice_cart_z = target[2] + random.uniform(-displacement_z, +displacement_z)
      
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    #                            TWO LATTICES to be filled                                  #  
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    if len(filename_custom_lattice) == 2:
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            decide orientation in xy plane                             #
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      # rotate around z (random value between 0 and angle_lattice_rot)
      # the lattice center was already moved to (0 0 0)

      # generate rotation angle (between 0.0 and angle_lattice_rot degree)
      rot_angle = [ random.uniform(0.0, angle_lattice_rot[0]*math.pi/180.0) , \
                    random.uniform(0.0, angle_lattice_rot[1]*math.pi/180.0) ] # angle in radians

      # rotate lattice 1
      lattice_1_tmp = copy.deepcopy(lattice_1)
      for latticepoint in lattice_1_tmp:
        # rotation axis --> z
        rot_axis = np.array([ 0.0 ,
                              0.0 ,
                              1.0 ])  

        # in cartesian coords
        x = latticepoint[0] * lattice_v_1[coord_X] + \
            latticepoint[1] * lattice_v_2[coord_X] + \
            latticepoint[2] * lattice_v_3[coord_X]
        
        y = latticepoint[0] * lattice_v_1[coord_Y] + \
            latticepoint[1] * lattice_v_2[coord_Y] + \
            latticepoint[2] * lattice_v_3[coord_Y]

        z = latticepoint[0] * lattice_v_1[coord_Z] + \
            latticepoint[1] * lattice_v_2[coord_Z] + \
            latticepoint[2] * lattice_v_3[coord_Z]
        
        before_rot = [ x, y, z ]

        after_rot = np.dot(rot.rotation_matrix(rot_axis, rot_angle[0]), before_rot) 
        after_rot_frac = tc.cart2frac(lattice_v_1, lattice_v_2, lattice_v_3, after_rot)

        latticepoint[0] = after_rot_frac[0]
        latticepoint[1] = after_rot_frac[1]
        latticepoint[2] = after_rot_frac[2]
      
      # rotate lattice 2
      lattice_2_tmp = copy.deepcopy(lattice_2)
      for latticepoint in lattice_2_tmp:
        # rotation axis --> z
        rot_axis = np.array([ 0.0 ,
                              0.0 ,
                              1.0 ])  

        # in cartesian coords
        x = latticepoint[0] * lattice_v_1[coord_X] + \
            latticepoint[1] * lattice_v_2[coord_X] + \
            latticepoint[2] * lattice_v_3[coord_X]
        
        y = latticepoint[0] * lattice_v_1[coord_Y] + \
            latticepoint[1] * lattice_v_2[coord_Y] + \
            latticepoint[2] * lattice_v_3[coord_Y]

        z = latticepoint[0] * lattice_v_1[coord_Z] + \
            latticepoint[1] * lattice_v_2[coord_Z] + \
            latticepoint[2] * lattice_v_3[coord_Z]
        
        before_rot = [ x, y, z ]
        after_rot = np.dot(rot.rotation_matrix(rot_axis, rot_angle[1]), before_rot) 
        after_rot_frac = tc.cart2frac(lattice_v_1, lattice_v_2, lattice_v_3, after_rot)

        latticepoint[0] = after_rot_frac[0]
        latticepoint[1] = after_rot_frac[1]
        latticepoint[2] = after_rot_frac[2]

      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            decide position of both lattices                           #  
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      displacement_lattice_cart_x = 0
      displacement_lattice_cart_y = 0
      displacement_lattice_cart_z = 0

      # ***** BOTH LATTICES TO RANDOM POSITION (independend of each other) *****
      if flag_lattice_location == RANDOM:
        position_frac = [ [ random.uniform(0.0, 1.0) , random.uniform(0.0, 1.0) , random.uniform(0.0, 1.0) ] ,
                          [ random.uniform(0.0, 1.0) , random.uniform(0.0, 1.0) , random.uniform(0.0, 1.0) ] ]
        
        # calculate displacement in cartesian coords
        displacement_lattice_cart_x = [ position_frac[0][coord_X] * lattice_v_1[coord_X] + \
                                        position_frac[0][coord_Y] * lattice_v_2[coord_X] + \
                                        position_frac[0][coord_Z] * lattice_v_3[coord_X] , \
                                        position_frac[1][coord_X] * lattice_v_1[coord_X] + \
                                        position_frac[1][coord_Y] * lattice_v_2[coord_X] + \
                                        position_frac[1][coord_Z] * lattice_v_3[coord_X] ]

        displacement_lattice_cart_y = [ position_frac[0][coord_X] * lattice_v_1[coord_Y] + \
                                        position_frac[0][coord_Y] * lattice_v_2[coord_Y] + \
                                        position_frac[0][coord_Z] * lattice_v_3[coord_Y] , \
                                        position_frac[1][coord_X] * lattice_v_1[coord_Y] + \
                                        position_frac[1][coord_Y] * lattice_v_2[coord_Y] + \
                                        position_frac[1][coord_Z] * lattice_v_3[coord_Y] ]

        displacement_lattice_cart_z = [ position_frac[0][coord_X] * lattice_v_1[coord_Z] + \
                                        position_frac[0][coord_Y] * lattice_v_2[coord_Z] + \
                                        position_frac[0][coord_Z] * lattice_v_3[coord_Z] , \
                                        position_frac[1][coord_X] * lattice_v_1[coord_Z] + \
                                        position_frac[1][coord_Y] * lattice_v_2[coord_Z] + \
                                        position_frac[1][coord_Z] * lattice_v_3[coord_Z] ]
        

      # ***** BOTH LATTICES TO TARGET POSITION *****
      elif flag_lattice_location == CENTERED:
        # the user has to specify the center of the lattice as well as the target location of the lattice
        # ==> with this info the lattice will be placed onto the target location with some random displacement

        # --- CAREFUL, right now only one cluster replica is supported ---
        # the lattice center right now is at (0 0 0) [ from previous rotation ]

        # there are two options, both lattices can be displaced for the exact same amount or both can be displaced independly
        # lattices displaced independendly
        if flag_link_lattice_displacement == False:
          displacement_lattice_cart_x = [ list_target_positions_1[0] + random.uniform(-displacement_x[0], +displacement_x[0]) , \
                                          list_target_positions_2[0] + random.uniform(-displacement_x[1], +displacement_x[1]) ]
          displacement_lattice_cart_y = [ list_target_positions_1[1] + random.uniform(-displacement_y[0], +displacement_y[0]) , \
                                          list_target_positions_2[1] + random.uniform(-displacement_y[1], +displacement_y[1]) ]
          displacement_lattice_cart_z = [ list_target_positions_1[2] + random.uniform(-displacement_z[0],  +displacement_z[0]) , \
                                          list_target_positions_2[2] + random.uniform(-displacement_z[1],  +displacement_z[1]) ]

        # lattices displaced in the same way
        # in this scenario, only displacement from lattice 1 counts
        elif flag_link_lattice_displacement == True:
          displacement = [ random.uniform(-displacement_x[0], +displacement_x[0]), \
                           random.uniform(-displacement_y[0], +displacement_y[0]), \
                           random.uniform(-displacement_z[0], +displacement_z[0]), ]

          displacement_lattice_cart_x = [ list_target_positions_1[0] + displacement[0] , \
                                          list_target_positions_2[0] + displacement[0] ]
          displacement_lattice_cart_y = [ list_target_positions_1[1] + displacement[1] , \
                                          list_target_positions_2[1] + displacement[1] ]
          displacement_lattice_cart_z = [ list_target_positions_1[2] + displacement[2] , \
                                          list_target_positions_2[2] + displacement[2] ]

                                        
      
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    #                            place the n_water molecules on lattice                     #  
    #-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#-#
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            ONLY ONE LATTICE                                           #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    if len(filename_custom_lattice) == 1: 
      # deep copy lattice into temporary lattice (used for this water structure only
      # in the deep copy, elements of the lattice which have been used, will be removed,
      # so that no latticepoint can be occupied twice
      # [DEBUG - NOW DONE BEFORE] lattice_tmp = copy.deepcopy(lattice)

      #--------------------------------#
      # MODE_LOCK_INDEX lattice coords #
      #--------------------------------#
      # if we are in MODE_LOCK_INDEX, then we need to figure out the location of the latticepoint
      # to which we want to lock the overlayer. 
      # This has to be done BEFORE any latticepoints are removed and latticepoints are displaced individually
      if filename_extra and extra_file_mode == "MODE_LOCK_INDEX":
          
          for lp_idx, latticepoint in enumerate(lattice_tmp):
            # get coordinates of latticepoint with specified index
            if lp_idx == idx_lock_lattice:
              LP_cart_x = latticepoint[0] * lattice_v_1[coord_X] + \
                          latticepoint[1] * lattice_v_2[coord_X] + \
                          latticepoint[2] * lattice_v_3[coord_X]
              
              LP_cart_y = latticepoint[0] * lattice_v_1[coord_Y] + \
                          latticepoint[1] * lattice_v_2[coord_Y] + \
                          latticepoint[2] * lattice_v_3[coord_Y]

              LP_cart_z = latticepoint[0] * lattice_v_1[coord_Z] + \
                          latticepoint[1] * lattice_v_2[coord_Z] + \
                          latticepoint[2] * lattice_v_3[coord_Z]
              
              # add in displacement of entire lattice
              LP_cart_x += displacement_lattice_cart_x
              LP_cart_y += displacement_lattice_cart_y
              LP_cart_z += displacement_lattice_cart_z
              
              coord_idx_lock_lattice = [LP_cart_x, LP_cart_y, LP_cart_z] # cartesian coords of lock-point in lattice
              break
  
      #--------------------------------------------------#
      # detect latticepoints that collide with substrate #
      #--------------------------------------------------#
      if flag_REMOVE_COLLISIONS == True:
        # loop over all latticepoints and remove if necessary
        list_lp_remove = []

        # --- find x, y and z boundaries of lattice ---
        x_min =  9999999
        x_max = -9999999
        y_min =  9999999
        y_max = -9999999
        z_min =  9999999
        z_max = -9999999
        for latticepoint in lattice_tmp:
          LP_cart_x = latticepoint[0] * lattice_v_1[coord_X] + \
                      latticepoint[1] * lattice_v_2[coord_X] + \
                      latticepoint[2] * lattice_v_3[coord_X]
          
          LP_cart_y = latticepoint[0] * lattice_v_1[coord_Y] + \
                      latticepoint[1] * lattice_v_2[coord_Y] + \
                      latticepoint[2] * lattice_v_3[coord_Y]

          LP_cart_z = latticepoint[0] * lattice_v_1[coord_Z] + \
                      latticepoint[1] * lattice_v_2[coord_Z] + \
                      latticepoint[2] * lattice_v_3[coord_Z]

          # add in displacement of entire lattice
          LP_cart_x += displacement_lattice_cart_x
          LP_cart_y += displacement_lattice_cart_y
          LP_cart_z += displacement_lattice_cart_z

          if LP_cart_x < x_min:
            x_min = LP_cart_x
          if LP_cart_y < y_min:
            y_min = LP_cart_y
          if LP_cart_z < z_min:
            z_min = LP_cart_z
          
          if LP_cart_x > x_max:
            x_max = LP_cart_x
          if LP_cart_y > y_max:
            y_max = LP_cart_y
          if LP_cart_z > z_max:
            z_max = LP_cart_z

        # --- find colliding latticepoints ---
        for atom_substrate in list_slab_atoms:
          # if substrate-atom outside boundaries: move on
          if atom_substrate[4] < x_min-LATTICE_COLLISION_TRESHHOLD/10.0:
            continue
          if atom_substrate[5] < y_min-LATTICE_COLLISION_TRESHHOLD/10.0:
            continue
          if atom_substrate[6] < z_min-LATTICE_COLLISION_TRESHHOLD/10.0:
            continue
          
          if atom_substrate[4] > x_max+LATTICE_COLLISION_TRESHHOLD/10.0:
            continue
          if atom_substrate[5] > y_max+LATTICE_COLLISION_TRESHHOLD/10.0:
            continue
          if atom_substrate[6] > z_max+LATTICE_COLLISION_TRESHHOLD/10.0:
            continue

          # if substrate atom within range, see if any latticepoint collides
          for lp_idx, latticepoint in enumerate(lattice_tmp):
            LP_cart_x = latticepoint[0] * lattice_v_1[coord_X] + \
                        latticepoint[1] * lattice_v_2[coord_X] + \
                        latticepoint[2] * lattice_v_3[coord_X]
            
            LP_cart_y = latticepoint[0] * lattice_v_1[coord_Y] + \
                        latticepoint[1] * lattice_v_2[coord_Y] + \
                        latticepoint[2] * lattice_v_3[coord_Y]

            LP_cart_z = latticepoint[0] * lattice_v_1[coord_Z] + \
                        latticepoint[1] * lattice_v_2[coord_Z] + \
                        latticepoint[2] * lattice_v_3[coord_Z]

            # add in displacement of entire lattice
            LP_cart_x += displacement_lattice_cart_x
            LP_cart_y += displacement_lattice_cart_y
            LP_cart_z += displacement_lattice_cart_z

            # calculate PBC distance between atom_substrate and latticepoint [UNITS: NM]
            d = pbc_distance ( p1 = [ LP_cart_x, LP_cart_y, LP_cart_z ]          , \
                               p2 = [ atom_substrate[4], atom_substrate[5], atom_substrate[6] ] , \
                               lattice_v_1 = lattice_v_1 , \
                               lattice_v_2 = lattice_v_2 , \
                               lattice_v_3 = lattice_v_3 )

            # treshhold input in AA, but gromacs: nm
            if d < LATTICE_COLLISION_TRESHHOLD/10.0:
              list_lp_remove.append(lp_idx)

      #--------------------------------------------------#
      # remove latticepoints that collide with substrate #
      #--------------------------------------------------#
      tmp_list = sorted(list(set(list_lp_remove)))
      list_lp_remove = tmp_list
      # loop through list backwards
      for idx_remove in reversed(list_lp_remove):
        del lattice_tmp[idx_remove]

      # because we potentially removed colliding lattice-points, there now might be more 
      # n_h2o than actual latticepoints --> in that case fill up every latticepoint
      n_h2o_tmp = n_h2o
      if n_h2o > len(lattice_tmp):
        n_h2o_tmp = len(lattice_tmp)

      # we have the option to not specify total number of h2o, but how many
      # latticepoints we want to keep empty
      #   --> user input for this: n_h2o < 0
      if n_h2o < 0:
        n_h2o_tmp = len(lattice_tmp) + n_h2o

      for i in range(0, n_h2o_tmp):
        #--------------------------------------------------#
        # decide which latticepoint this water should take #
        #--------------------------------------------------#
        # elements which were used of lattice are removed in lattice_tmp
        latticepoint_index = -1

        # *** generate random integer ***
        if flag_lattice_point_occupation == RANDOM_OCCUPATION:
          # if lattice is gonna be filled completely --> not random choise
          if flag_lattice_filling == FILLED_COMPLETELY:
            latticepoint_index = 0 # always fill the first element (which will be popped afterwards)
          elif flag_lattice_filling == FILLED_PARTIALLY:
            latticepoint_index = random.randint(0, len(lattice_tmp)-1)
          else:
            print "[ERROR]. flag_lattice_filling =", flag_lattice_filling, "is not known."
            sys.exit()
        
        # *** equal distribution in lattice sites ***
        elif flag_lattice_point_occupation == EQUAL_OCCUPATION:
          # how many structures for each lattice points?
          n_struct_per_latticepoint = int(n_structures / len(lattice))

          # decide which lattice-point to use
          latticepoint_index = int(n/n_struct_per_latticepoint) 
        
        # remove the chosen latticepoint so it does not get occupied twice  
        latticepoint = lattice_tmp.pop(latticepoint_index)

        # in cartesian coords
        displacement_latticepoint_cart_x = latticepoint[0] * lattice_v_1[coord_X] + \
                                           latticepoint[1] * lattice_v_2[coord_X] + \
                                           latticepoint[2] * lattice_v_3[coord_X]
        
        displacement_latticepoint_cart_y = latticepoint[0] * lattice_v_1[coord_Y] + \
                                           latticepoint[1] * lattice_v_2[coord_Y] + \
                                           latticepoint[2] * lattice_v_3[coord_Y]

        displacement_latticepoint_cart_z = latticepoint[0] * lattice_v_1[coord_Z] + \
                                           latticepoint[1] * lattice_v_2[coord_Z] + \
                                           latticepoint[2] * lattice_v_3[coord_Z]


        #--------------------#
        # decide orientation # of water molecule (random: axis and angle)
        #--------------------#
        # generate rotation axis (randomly in x,y and z)
        rot_axis = np.array([ random.uniform(-1.0,1.0) ,
                              random.uniform(-1.0,1.0) ,
                              random.uniform(-1.0,1.0) ])  
      
        # generate rotation angle
        rot_angle = random.uniform(0.0, 2*math.pi) # angle in radians

        # initialize
        oxygen_atom     = []
        hydrogen_atom_1 = []
        hydrogen_atom_2 = []
        virtual_site    = []
        
        if water_model == SPC:
          oxygen_atom = [ molecule_number,
                          "SOL",
                          "OW",
                          atom_index+1,
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[0],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[1],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[2]  ]
          
          hydrogen_atom_1 = [ molecule_number,
                              "SOL",
                              "HW1",
                              atom_index+2,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[2]  ]
          
          hydrogen_atom_2 = [ molecule_number,
                              "SOL",
                              "HW2",
                              atom_index+3,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[2]  ]
          
          atom_number     += 3
 
        elif water_model == TIP4P_ICE:
          oxygen_atom = [ molecule_number,
                          "SOL",
                          "OW",
                          atom_index+1,
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_OW)[0],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_OW)[1],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_OW)[2]  ]
          
          hydrogen_atom_1 = [ molecule_number,
                              "SOL",
                              "HW1",
                              atom_index+2,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW1)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW1)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW1)[2]  ]
          
          hydrogen_atom_2 = [ molecule_number,
                              "SOL",
                              "HW2",
                              atom_index+3,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW2)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW2)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW2)[2]  ]
          virtual_site   = [ molecule_number, 
                              "SOL",
                              "MW",
                              atom_index+4,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_MW)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_MW)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_MW)[2]  ]
      
          atom_number     += 4
        else:
          print "[ERROR]. Currently only SPC and TIP4P/ice supported to fill lattice points"
          sys.exit()
        
        molecule_number += 1

        #---------------------------------------#
        # decide displacement from latticepoint #
        #---------------------------------------#
        extra_displacement_lattice_cart_x = random.uniform(-extra_displacement_cart_xy, +extra_displacement_cart_xy)
        extra_displacement_lattice_cart_y = random.uniform(-extra_displacement_cart_xy, +extra_displacement_cart_xy)
        extra_displacement_lattice_cart_z = random.uniform(-extra_displacement_cart_z,  +extra_displacement_cart_z)
        
        #-----------------------------------#
        # move atom to its location on slab #
        #-----------------------------------#
        #                 whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
        oxygen_atom[4] += displacement_lattice_cart_x + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
        oxygen_atom[5] += displacement_lattice_cart_y + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
        oxygen_atom[6] += displacement_lattice_cart_z + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             
        
        #                      whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
        hydrogen_atom_1[4] +=  displacement_lattice_cart_x + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
        hydrogen_atom_1[5] +=  displacement_lattice_cart_y + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
        hydrogen_atom_1[6] +=  displacement_lattice_cart_z + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             
   
        #                     whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
        hydrogen_atom_2[4] += displacement_lattice_cart_x + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
        hydrogen_atom_2[5] += displacement_lattice_cart_y + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
        hydrogen_atom_2[6] += displacement_lattice_cart_z + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             
        
        if water_model == TIP4P_ICE:
          #                   whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
          virtual_site[4] += displacement_lattice_cart_x + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
          virtual_site[5] += displacement_lattice_cart_y + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
          virtual_site[6] += displacement_lattice_cart_z + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             
        
        #--------------------------#
        # generate water structure #
        #--------------------------#
        # list_water.append([oxygen_atom, hydrogen_atom_1, hydrogen_atom_2])

        list_water.append(oxygen_atom)
        list_water.append(hydrogen_atom_1)
        list_water.append(hydrogen_atom_2)

        if water_model == TIP4P_ICE:
          list_water.append(virtual_site)

    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            TWO LATTICES                                               #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 

    if len(filename_custom_lattice) == 2: 
      # deep copy lattice into temporary lattice (used for this water structure only
      # in the deep copy, elements of the lattice which have been used, will be removed,
      # so that no latticepoint can be occupied twice
      # [DEBUG - NOW DONE PREVIOUSLY] lattice_1_tmp = copy.deepcopy(lattice_1)
      # [DEBUG - NOW DONE PREVIOUSLY] lattice_2_tmp = copy.deepcopy(lattice_2)

      #----------------#
      # security check #
      #----------------#
      if flag_REMOVE_COLLISIONS == True:
        print "[ERROR]. Collision removal is only supported for 1 lattice at the moment."
        print "         Come back and implement this if necessary"
        sys.exit()

      #---------------------------------------#
      # figure out a way how to fill lattices #
      #---------------------------------------#
      n_h2o_1 = 0
      n_h2o_2 = 0

      # fill lattice 1 completely, the rest goes onto lattice 2
      if flag_fill_lattice_method == FILL_LATTICE_1_FULLY:
        n_h2o_1 = len(lattice_1)
        n_h2o_2 = n_h2o - n_h2o_1
      elif flag_fill_lattice_method == FILL_LATTICE_2_FULLY:
        n_h2o_2 = len(lattice_2)
        n_h2o_1 = n_h2o - n_h2o_2
      
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            fill lattice_1                                             #  
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      for i in range(0, n_h2o_1):
        #--------------------------------------------------#
        # decide which latticepoint this water should take #
        #--------------------------------------------------#
        # elements which were used of lattice are removed in lattice_tmp
        # generate random integer
        
        # [OLD, DEL?] latticepoint_index = random.randint(0, len(lattice_1_tmp)-1)
        
        # if lattice is gonna be filled completely --> not random choise
        if flag_lattice_filling == FILLED_COMPLETELY:
          latticepoint_index = 0 # always fill the first element (which will be popped afterwards)
        elif flag_lattice_filling == FILLED_PARTIALLY:
          latticepoint_index = random.randint(0, len(lattice_1_tmp)-1)
        else:
          print "[ERROR]. flag_lattice_filling =", flag_lattice_filling, "is not known."
          sys.exit()
        

        latticepoint = lattice_1_tmp.pop(latticepoint_index)

        # in cartesian coords
        displacement_latticepoint_cart_x = latticepoint[0] * lattice_v_1[coord_X] + \
                                           latticepoint[1] * lattice_v_2[coord_X] + \
                                           latticepoint[2] * lattice_v_3[coord_X]
        
        displacement_latticepoint_cart_y = latticepoint[0] * lattice_v_1[coord_Y] + \
                                           latticepoint[1] * lattice_v_2[coord_Y] + \
                                           latticepoint[2] * lattice_v_3[coord_Y]

        displacement_latticepoint_cart_z = latticepoint[0] * lattice_v_1[coord_Z] + \
                                           latticepoint[1] * lattice_v_2[coord_Z] + \
                                           latticepoint[2] * lattice_v_3[coord_Z]

        #--------------------#
        # decide orientation # (random: axis and angle)
        #--------------------#
        # generate rotation axis (randomly in x,y and z)
        rot_axis = np.array([ random.uniform(-1.0,1.0) ,
                              random.uniform(-1.0,1.0) ,
                              random.uniform(-1.0,1.0) ])  
      
        # generate rotation angle
        rot_angle = random.uniform(0.0, 2*math.pi) # angle in radians

        
        # initialize
        oxygen_atom     = []
        hydrogen_atom_1 = []
        hydrogen_atom_2 = []
        virtual_site    = []
        
        if water_model == SPC:
          oxygen_atom = [ molecule_number,
                          "SOL",
                          "OW",
                          atom_index+1,
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[0],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[1],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[2]  ]
          
          hydrogen_atom_1 = [ molecule_number,
                              "SOL",
                              "HW1",
                              atom_index+2,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[2]  ]
          
          hydrogen_atom_2 = [ molecule_number,
                              "SOL",
                              "HW2",
                              atom_index+3,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[2]  ]
      
          atom_number     += 3
        elif water_model == TIP4P_ICE:
          oxygen_atom = [ molecule_number,
                          "SOL",
                          "OW",
                          atom_index+1,
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_OW)[0],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_OW)[1],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_OW)[2]  ]
          
          hydrogen_atom_1 = [ molecule_number,
                              "SOL",
                              "HW1",
                              atom_index+2,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW1)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW1)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW1)[2]  ]
          
          hydrogen_atom_2 = [ molecule_number,
                              "SOL",
                              "HW2",
                              atom_index+3,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW2)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW2)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW2)[2]  ]
          
          virtual_site   = [ molecule_number, 
                              "SOL",
                              "MW",
                              atom_index+4,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_MW)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_MW)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_MW)[2]  ]
      
          atom_number     += 4
        else:
          print "[ERROR]. Currently only SPC and TIP4P/ice supported to fill lattice points"
          sys.exit()
        
        
        molecule_number += 1

        #---------------------------------------#
        # decide displacement from latticepoint #
        #---------------------------------------#
        extra_displacement_lattice_cart_x = random.uniform(-extra_displacement_cart_xy[0], + extra_displacement_cart_xy[0])
        extra_displacement_lattice_cart_y = random.uniform(-extra_displacement_cart_xy[0], + extra_displacement_cart_xy[0])
        extra_displacement_lattice_cart_z = random.uniform(-extra_displacement_cart_z[0], +extra_displacement_cart_z[0])
        
        #-----------------------------------#
        # move atom to its location on slab #
        #-----------------------------------#

        #                   whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
        oxygen_atom[4] += displacement_lattice_cart_x[0] + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
        oxygen_atom[5] += displacement_lattice_cart_y[0] + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
        oxygen_atom[6] += displacement_lattice_cart_z[0] + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             
        
        #                      whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
        hydrogen_atom_1[4] +=  displacement_lattice_cart_x[0] + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
        hydrogen_atom_1[5] +=  displacement_lattice_cart_y[0] + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
        hydrogen_atom_1[6] +=  displacement_lattice_cart_z[0] + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             
   
        #                     whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
        hydrogen_atom_2[4] += displacement_lattice_cart_x[0] + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
        hydrogen_atom_2[5] += displacement_lattice_cart_y[0] + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
        hydrogen_atom_2[6] += displacement_lattice_cart_z[0] + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             

        if water_model == TIP4P_ICE:
          #                   whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
          virtual_site[4] += displacement_lattice_cart_x[0] + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
          virtual_site[5] += displacement_lattice_cart_y[0] + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
          virtual_site[6] += displacement_lattice_cart_z[0] + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             
        
        #--------------------------#
        # generate water structure #
        #--------------------------#
        # list_water.append([oxygen_atom, hydrogen_atom_1, hydrogen_atom_2])
        list_water.append(oxygen_atom)
        list_water.append(hydrogen_atom_1)
        list_water.append(hydrogen_atom_2)
        
        if water_model == TIP4P_ICE:
          list_water.append(virtual_site)

      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            fill lattice_2                                             #  
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      for i in range(0, n_h2o_2):
        #--------------------------------------------------#
        # decide which latticepoint this water should take #
        #--------------------------------------------------#
        # elements which were used of lattice are removed in lattice_tmp
        # generate random integer
        
        # [OLD, DEL?] latticepoint_index = random.randint(0, len(lattice_2_tmp)-1)

        # if lattice is gonna be filled completely --> not random choise
        if flag_lattice_filling == FILLED_COMPLETELY:
          latticepoint_index = 0 # always fill the first element (which will be popped afterwards)
        elif flag_lattice_filling == FILLED_PARTIALLY:
          latticepoint_index = random.randint(0, len(lattice_2_tmp)-1)
        else:
          print "[ERROR]. flag_lattice_filling =", flag_lattice_filling, "is not known."
          sys.exit()

        latticepoint = lattice_2_tmp.pop(latticepoint_index)

        # in cartesian coords
        displacement_latticepoint_cart_x = latticepoint[0] * lattice_v_1[coord_X] + \
                                           latticepoint[1] * lattice_v_2[coord_X] + \
                                           latticepoint[2] * lattice_v_3[coord_X]
        
        displacement_latticepoint_cart_y = latticepoint[0] * lattice_v_1[coord_Y] + \
                                           latticepoint[1] * lattice_v_2[coord_Y] + \
                                           latticepoint[2] * lattice_v_3[coord_Y]

        displacement_latticepoint_cart_z = latticepoint[0] * lattice_v_1[coord_Z] + \
                                           latticepoint[1] * lattice_v_2[coord_Z] + \
                                           latticepoint[2] * lattice_v_3[coord_Z]

        #--------------------#
        # decide orientation # (random: axis and angle)
        #--------------------#
        # generate rotation axis (randomly in x,y and z)
        rot_axis = np.array([ random.uniform(-1.0,1.0) ,
                              random.uniform(-1.0,1.0) ,
                              random.uniform(-1.0,1.0) ])  
      
        # generate rotation angle
        rot_angle = random.uniform(0.0, 2*math.pi) # angle in radians

        # initialize
        oxygen_atom     = []
        hydrogen_atom_1 = []
        hydrogen_atom_2 = []
        virtual_site    = []
        
        if water_model == SPC:
          oxygen_atom = [ molecule_number,
                          "SOL",
                          "OW",
                          atom_index+1,
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[0],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[1],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_OW)[2]  ]
          
          hydrogen_atom_1 = [ molecule_number,
                              "SOL",
                              "HW1",
                              atom_index+2,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW1)[2]  ]
          
          hydrogen_atom_2 = [ molecule_number,
                              "SOL",
                              "HW2",
                              atom_index+3,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), SPC_HW2)[2]  ]
      
          atom_number     += 3
        elif water_model == TIP4P_ICE:
          oxygen_atom = [ molecule_number,
                          "SOL",
                          "OW",
                          atom_index+1,
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_OW)[0],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_OW)[1],
                          np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_OW)[2]  ]
          
          hydrogen_atom_1 = [ molecule_number,
                              "SOL",
                              "HW1",
                              atom_index+2,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW1)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW1)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW1)[2]  ]
          
          hydrogen_atom_2 = [ molecule_number,
                              "SOL",
                              "HW2",
                              atom_index+3,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW2)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW2)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_HW2)[2]  ]
          
          virtual_site   = [ molecule_number, 
                              "SOL",
                              "MW",
                              atom_index+4,
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_MW)[0],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_MW)[1],
                              np.dot(rot.rotation_matrix(rot_axis, rot_angle), TIP4Pice_MW)[2]  ]
      
          atom_number     += 4
        else:
          print "[ERROR]. Currently only SPC and TIP4P/ice supported to fill lattice points"
          sys.exit()
        
        
        molecule_number += 1

        #---------------------------------------#
        # decide displacement from latticepoint #
        #---------------------------------------#
        extra_displacement_lattice_cart_x = random.uniform(-extra_displacement_cart_xy[1], + extra_displacement_cart_xy[1])
        extra_displacement_lattice_cart_y = random.uniform(-extra_displacement_cart_xy[1], + extra_displacement_cart_xy[1])
        extra_displacement_lattice_cart_z = random.uniform(-extra_displacement_cart_z[1], +extra_displacement_cart_z[1])
        
        #-----------------------------------#
        # move atom to its location on slab #
        #-----------------------------------#
        #                   whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
        oxygen_atom[4] += displacement_lattice_cart_x[1] + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
        oxygen_atom[5] += displacement_lattice_cart_y[1] + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
        oxygen_atom[6] += displacement_lattice_cart_z[1] + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             
        
        #                      whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
        hydrogen_atom_1[4] +=  displacement_lattice_cart_x[1] + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
        hydrogen_atom_1[5] +=  displacement_lattice_cart_y[1] + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
        hydrogen_atom_1[6] +=  displacement_lattice_cart_z[1] + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             
   
        #                     whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
        hydrogen_atom_2[4] += displacement_lattice_cart_x[1] + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
        hydrogen_atom_2[5] += displacement_lattice_cart_y[1] + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
        hydrogen_atom_2[6] += displacement_lattice_cart_z[1] + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             
        
        if water_model == TIP4P_ICE:
          #                   whole lattice displacement     individual lattice point displ     extra displacement around latticepoint
          virtual_site[4] += displacement_lattice_cart_x[0] + displacement_latticepoint_cart_x + extra_displacement_lattice_cart_x 
          virtual_site[5] += displacement_lattice_cart_y[0] + displacement_latticepoint_cart_y + extra_displacement_lattice_cart_y
          virtual_site[6] += displacement_lattice_cart_z[0] + displacement_latticepoint_cart_z + extra_displacement_lattice_cart_z             
        
        #--------------------------#
        # generate water structure #
        #--------------------------#
        # list_water.append([oxygen_atom, hydrogen_atom_1, hydrogen_atom_2])
        list_water.append(oxygen_atom)
        list_water.append(hydrogen_atom_1)
        list_water.append(hydrogen_atom_2)

        if water_model == TIP4P_ICE:
          list_water.append(virtual_site)
  

  ##########################################################################################
  #                            POSITION EXTRA_STRUCTURE IF NECESSARY                       # 
  ##########################################################################################  
  if filename_extra:
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            EXPLICIT MODE                                              #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    if extra_file_mode == "MODE_EXPLICIT":
      # decide displacement x/y/z
      displacement_extra_struct_x = random.uniform(-extra_file_displ_x, +extra_file_displ_x)
      displacement_extra_struct_y = random.uniform(-extra_file_displ_y, +extra_file_displ_y)
      displacement_extra_struct_z = random.uniform(-extra_file_displ_z, +extra_file_displ_z)

      
      #--------------------------------------------------#
      # check if the structure falls within z_lock range #
      #--------------------------------------------------#
      # --> if lattice & extra_file are too far away or too close (in z), then disregard this struct and generate new one
      # find z_max of lattice
      lattice_z_max = -999999
      for atom in list_water:
        if atom[2] == "OW" and atom[Z_CART] > lattice_z_max:
          lattice_z_max = atom[Z_CART]

      if (z_min_extra_file+displacement_extra_struct_z/10.0)-lattice_z_max < min_z_lock/10.0:
        # [DEBUG] print "TOO CLOSE", 10*((z_min_extra_file+displacement_extra_struct_z/10.0)-lattice_z_max)
        n -= 1  # generate new structure, because this one is bad
        continue
      if (z_min_extra_file+displacement_extra_struct_z/10.0)-lattice_z_max > max_z_lock/10.0:
        # [DEBUG] print "TOO FAR", 10*((z_min_extra_file+displacement_extra_struct_z/10.0)-lattice_z_max)
        n -= 1 # generate new structure because this one is bad
        continue
      
      # displace all extra atoms accordingly
      file_extra_moved = copy.deepcopy(file_extra)
      for atom in file_extra_moved.list_atoms:
        atom[file_extra_moved.X_COORD] += displacement_extra_struct_x/10.0
        atom[file_extra_moved.Y_COORD] += displacement_extra_struct_y/10.0
        atom[file_extra_moved.Z_COORD] += displacement_extra_struct_z/10.0

    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            INDEX_LOCK                                                 #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    if extra_file_mode == "MODE_LOCK_INDEX":
      
      file_extra_moved = copy.deepcopy(file_extra)
      #---------------------------------------------------------------#
      # rotate extra-file (around z) with the same angle than lattice #
      #---------------------------------------------------------------#
      for idx_atom_extra, atom in enumerate(file_extra_moved.list_atoms):
        # rotation axis --> z
        rot_axis = np.array([ 0.0 ,
                              0.0 ,
                              1.0 ])  

        # coords are in nm at the moment, but this doesn't matter for rotation...
        # also, the center of rotation does not matter here, because we will lock indizes together afterwards...
        before_rot = [ atom[file_extra_moved.X_COORD], atom[file_extra_moved.Y_COORD], atom[file_extra_moved.Z_COORD] ]
        after_rot = np.dot(rot.rotation_matrix(rot_axis, lattice_rotation_angle), before_rot) 
      
        atom[file_extra_moved.X_COORD] = after_rot[0]
        atom[file_extra_moved.Y_COORD] = after_rot[1]
        atom[file_extra_moved.Z_COORD] = after_rot[2]

        if idx_atom_extra == idx_lock_extra:
          coord_idx_lock_extra = [ atom[file_extra_moved.X_COORD], atom[file_extra_moved.Y_COORD], atom[file_extra_moved.Z_COORD] ] # cartesian coords of lock-point in extra file

      #-----------------------#
      # lock indizes together #
      #-----------------------#
      # need to also rotate lock_dx, lock_dy, lock_dz
      lock_rot = np.dot(rot.rotation_matrix(rot_axis, lattice_rotation_angle), [lock_dx, lock_dy, lock_dz]) 

      # by displacing OL appropriately
      for atom in file_extra_moved.list_atoms:
        atom[file_extra_moved.X_COORD] += coord_idx_lock_lattice[0] - coord_idx_lock_extra[0] + lock_rot[0]/10.0
        atom[file_extra_moved.Y_COORD] += coord_idx_lock_lattice[1] - coord_idx_lock_extra[1] + lock_rot[1]/10.0
        atom[file_extra_moved.Z_COORD] += coord_idx_lock_lattice[2] - coord_idx_lock_extra[2] + lock_rot[2]/10.0

  ##########################################################################################
  #                            PUT WATER TOGETHER WITH SLAB                                # 
  ##########################################################################################  
  # we need to explicitly deepcopy them
  # else a change in list_atoms also affects list_slab_atoms / list_water
  list_atoms = copy.deepcopy(list_slab_atoms) + copy.deepcopy(list_water)
 
  ##########################################################################################
  #                            GENERATE SUPERCELL                                          # 
  ##########################################################################################  
  # print "\t(+) create %dx%dx%d supercell ..." % (multiply_x, multiply_y, multiply_z)
 
  #---------------------------------------------#
  # convert cartesian coords to fractional ones #
  #---------------------------------------------#
  for atom in list_atoms:
    coord = [ atom[X_CART], atom[Y_CART], atom[Z_CART] ]
    frac_coord = tc.cart2frac (lattice_v_1, lattice_v_2, lattice_v_3, coord)
  
    # save fractional coordinates
    atom.append(frac_coord[0])
    atom.append(frac_coord[1])
    atom.append(frac_coord[2])

  #--------------------------#
  # increase lattice vectors #
  #--------------------------#
  lattice_v_1_sc = [ coordinate*multiply_x for coordinate in lattice_v_1 ]
  lattice_v_2_sc = [ coordinate*multiply_y for coordinate in lattice_v_2 ]
  lattice_v_3_sc = [ coordinate*multiply_z for coordinate in lattice_v_3 ]
   
  # check if supercell is big enough for cutoff (1nm) --> ie 2nm long
  if lattice_v_1_sc[0] < 2.0 or lattice_v_2_sc[1] < 2.0:
    print "[ERROR]. The cell is too small for a 1nm cutoff"
    print sys.exit()

  #--------------------------------#
  # generate additional boxes in X #
  #--------------------------------#
  list_additional_boxes = []
  index_increase = 0
  if (multiply_x > 1):
    for i in range(1,multiply_x):
      index_increase += len(list_atoms)
      list_atoms_new_box = copy.deepcopy(list_atoms)
      for atom in list_atoms_new_box:
        atom[X_FRAC]      += i
        atom[INDEX]  += index_increase
      
      list_additional_boxes.append(list_atoms_new_box)
  
  #-----------------------------#
  # merge additional boxes in X #
  #-----------------------------#
  for box in list_additional_boxes:
    list_atoms += box
  
  #-------------------------------#
  # adjust fractional coordinates # (multiplying a cell by 2, divides the fractional coords by 2)
  #-------------------------------#
  for atom in list_atoms:
    atom[X_FRAC] = atom[X_FRAC] / multiply_x
  
  #--------------------------------#
  # generate additional boxes in Y #
  #--------------------------------#
  list_additional_boxes = []
  index_increase = 0
  if (multiply_y > 1):
    for i in range(1,multiply_y):
      index_increase += len(list_atoms)
      list_atoms_new_box = copy.deepcopy(list_atoms)
      for atom in list_atoms_new_box:
        atom[Y_FRAC]     += i
        atom[INDEX] += index_increase
      
      list_additional_boxes.append(list_atoms_new_box)
  
  #-----------------------------#
  # merge additional boxes in Y #
  #-----------------------------#
  for box in list_additional_boxes:
    list_atoms += box
  
  #-------------------------------#
  # adjust fractional coordinates # (multiplying a cell by 2, divides the fractional coords by 2)
  #-------------------------------#
  for atom in list_atoms:
    atom[Y_FRAC] = atom[Y_FRAC] / multiply_y
  
  #--------------------------------#
  # generate additional boxes in Z #
  #--------------------------------#
  list_additional_boxes = []
  index_increase = 0
  if (multiply_z > 1):
    for i in range(1,multiply_z):
      index_increase += len(list_atoms)
      list_atoms_new_box = copy.deepcopy(list_atoms)
      for atom in list_atoms_new_box:
        atom[Z_FRAC]     += i
        atom[INDEX] += index_increase
      
      list_additional_boxes.append(list_atoms_new_box)
  
  #-----------------------------#
  # merge additional boxes in Z #
  #-----------------------------#
  for box in list_additional_boxes:
    list_atoms += box
  
  #-------------------------------#
  # adjust fractional coordinates # (multiplying a cell by 2, divides the fractional coords by 2)
  #-------------------------------#
  for atom in list_atoms:
    atom[Z_FRAC] = atom[Z_FRAC] / multiply_z

  #-------------------#
  # convert frac2cart #
  #-------------------#
  for atom in list_atoms:
    coord = tc.frac2cart(lattice_v_1_sc, lattice_v_2_sc, lattice_v_3_sc, [atom[X_FRAC], atom[Y_FRAC], atom[Z_FRAC]])

    atom[X_CART] = coord[0]
    atom[Y_CART] = coord[1]
    atom[Z_CART] = coord[2]


  ####################################################################
  #                          write gro output file                   #
  ####################################################################
  # rjust --> fill string up with "0" from right until length 4
  length_max_num = 5

  filename_output = filename_out + "_" + str(n).rjust(length_max_num, "0") + ".gro"
  file_output = open(filename_output,"w")
  file_output.write("created by " + sys.argv[0] + "\n")
  if filename_extra:
    file_output.write(str(len(list_atoms)+len(file_extra_moved.list_atoms)) + "\n")
  else:
    file_output.write(str(len(list_atoms)) + "\n")

  # re-order index (first slab, then SOL)
  atom_index = 1

  molecule_number = 0
  # --- write slab ---
  for atom in list_atoms:
    line = "" 
    if atom[2] != "OW" and atom[2] != "HW1" and atom[2] != "HW2" and atom[2] != "MW":
      # keep molecule number up to date
      if atom[0] > molecule_number:
        molecule_number = atom[0]

      line = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (atom[0], 
                                                 atom[1], 
                                                 atom[2], 
                                                 atom_index, 
                                                 atom[4],
                                                 atom[5], 
                                                 atom[6]  )
      atom_index += 1
    
    file_output.write(line)
  
  # --- write SOL ---
  for atom in list_atoms:
    line = "" 
    if atom[2] == "OW" or atom[2] == "HW1" or atom[2] == "HW2" or atom[2] == "MW":
      # increase molecule number counter for each "OW"
      if atom[2] == "OW":
        molecule_number += 1

      line = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (molecule_number, 
                                                 atom[1], 
                                                 atom[2], 
                                                 atom_index, 
                                                 atom[4],
                                                 atom[5], 
                                                 atom[6]  )
      atom_index += 1
    
    file_output.write(line)
  
  # --- write exta_file if necessary ---
  if filename_extra:
    for atom in file_extra_moved.list_atoms:
      line = "" 
      if atom[2] == "OW" or atom[2] == "HW1" or atom[2] == "HW2" or atom[2] == "MW":
        # increase molecule number counter for each "OW"
        if atom[2] == "OW":
          molecule_number += 1

        line = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % (molecule_number, 
                                                   atom[1], 
                                                   atom[2], 
                                                   atom_index, 
                                                   atom[4],
                                                   atom[5], 
                                                   atom[6]  )
        atom_index += 1
      
      file_output.write(line)

  # --- write box dimensions ---
  line = "%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n" % (lattice_v_1_sc[0],   \
                                                                               lattice_v_2_sc[1],   \
                                                                               lattice_v_3_sc[2],   \
                                                                               lattice_v_1_sc[1],   \
                                                                               lattice_v_1_sc[2],   \
                                                                               lattice_v_2_sc[0],   \
                                                                               lattice_v_2_sc[2],   \
                                                                               lattice_v_3_sc[0],   \
                                                                               lattice_v_3_sc[1]    )
  file_output.write(line)
  file_output.close()
