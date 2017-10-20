# HSeed - Heterogeneous Seeding Approach for Molecular Simulation #
# ############################################################### #
# 20.10.2017

BUILD is the first step.

A working (?) example can be found in .EXAMPLE/BUILD/EXAMPLE_2 (Metaldehyde)

It all starts with a .gro file, containing the surface only (no water). e.g. mdhe_100.gro

 Get a .POSCAR containing the ice face we want to seed (e.g. Ic_SF_001_READY.POSCAR in the FACES dir)

 VMD typically gets the box wrong! Cut a slice so that all water molecules are nice and whole via VMD, 
 save it in a new POSCAR, e.g. Ic_001_T30_CUT.POSCAR, sort them (fort. POSCAR to xyz)
 Double the title line (the first one!) just after cell dimensions in the resulting POSCAR

 Clean the poscar via tidy_up_POSCAR.py, selecting the z cutoffs
 e.g. tidy up: ./tidy_up_POSCAR.py Ic_001_T30.POSCAR Ic_001_T30_CLEAN.POSCAR 16 56
 Check number of O, H, and X. Check whether every O as Two H...

 Select the right face in make_hemisphere.py (comment/scomment relevant bits)
 Select the desired number of layers to be considered as "lattice" in make_hemisphere.py 
 by modyfying Z_LATTICE_CUTOFF and Z_LAYER_REMOVE. Evrything between Z_LAYER_REMOVE and Z_LATTICE_CUTOFF would 
 be considered as lattice (stuff that changes) as opposed to seed (fixed). Play around with those 
 two thrshhold in order to have reasonable guesses.
 ../SRC/make_hemisphere.py Ic_SF_001_READY.POSCAR seed_Ic_001.gro seed_Ic_001.lattice 500 True

 Try out seeds with empty cells : EMPTY_100_x_100.gro in SRC
 Specify the lattice file in generate_random_structure.py. e.g. seed_Ic_001.lattice
 Fix the cell if orthorombic
 Fix these parameters in the .lattice file:
 TARGET_X TARGET_Y TARGET_Z e.g. 
 DISPLACE_X DISPLACE_Y DISPLACE_Z
 ./SRC/generate_random_structure.py mdhe_100.gro seed_Ic_001.gro 112 rss 10
 threshold chosen to remove lattice points (LATTICE_COLLISION_TRESHHOLD)
 make sure you got the face you wanted
