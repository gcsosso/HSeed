import os.path
import math
import sys
import numpy as np

""" 
the functions included here are:
  (1) get height above cube
  (2) get lateral profile above cube
  (3) get height over slab
  (4) get lateral profile over slab
  (5) integrate height profile (slab)
  (6) integrate height profile (cube) --> TODO
  (7) 2-D probability distribution
  (8) free energy profile
"""

##########################################################################################
# (1)                        GET HEIGHT ABOVE CUBE                                       #
##########################################################################################
def height_from_cube( trajfilename, \
                      outfilename,  \
                      nframes,      \
                      natoms,       \
                      nwater,       \
                      x0,           \
                      y0,           \
                      z0,           \
                      l=20.0):
    """ 
    [ADAPTED FROM STEVE]

    Calculates the height of a water molecule from the cube provided
    that the water molecule resides above one of the six cube surfaces.

    trajfilename = .gro trajectory file
    outfilename  = output file - list of heights from the cube
    nframes      = self explanatory
    natoms       = number of atoms in trajectory
    nwater       = number of water molecules in trajectory
    x0,y0,z0     = centre of the cube
    l            = edge-length of cube

    """
    
    if os.path.isfile(outfilename+"_px.dat") or \
       os.path.isfile(outfilename+"_mx.dat") or \
       os.path.isfile(outfilename+"_py.dat") or \
       os.path.isfile(outfilename+"_my.dat") or \
       os.path.isfile(outfilename+"_pz.dat") or \
       os.path.isfile(outfilename+"_mz.dat") :
        # don't allow overwriting of files
        print "Not allowed to overwrite files!"
        print "attempted file to overwrite: %s" % outfilename
        sys.exit(1)
    else:
        trajfile = open(trajfilename,'r')

        # create six outfiles, one for each cube-face
        # names: +x (px), -x (mx), +y (py), -y (my), ..
        outfile_px  = open(outfilename+"_px.dat",'w')
        outfile_mx  = open(outfilename+"_mx.dat",'w')
        outfile_py  = open(outfilename+"_py.dat",'w')
        outfile_my  = open(outfilename+"_my.dat",'w')
        outfile_pz  = open(outfilename+"_pz.dat",'w')
        outfile_mz  = open(outfilename+"_mz.dat",'w')

        line_number = 1    
        for frame in range(nframes):
            # read the two header lines in the .gro file
            line = trajfile.readline()
            line = trajfile.readline()

            line_number += 2

            # construct empty lists containing the atom positions 
            xlist = []
            ylist = []
            zlist = []

            # now loop through each water
            for i in range(natoms):

                line = trajfile.readline()
                line_number += 1
                info = line.split()

                # sometimes (n>9999) 2nd and 3rd column not separated in gro
                if len(info) == 6:
                  if info[1][:2] == "mW":
                    x = float(info[3])
                    y = float(info[4])
                    z = float(info[5])
                  else:
                    continue
                elif len(info) == 5:
                  if info[1][:2] == "mW":
                    x = float(info[2])
                    y = float(info[3])
                    z = float(info[4])
                  else:
                    continue
                else:
                  print "[ERROR]. Can't read gro file,", trajfilename
                  print "line: ", line_number
                  print line
                  sys.exit()

                xlist.append(x)
                ylist.append(y)
                zlist.append(z)

            # check if all water molecules were read in
            if len(xlist) != len(ylist) != len(zlist) != nwater:
              print "[ERROR]. Could not read the specified number of water molecules"
              print len(xlist), len(ylist), len(zlist)
              sys.exit()

            # read in the cell-vectors
            line = trajfile.readline()
            info = line.split()
            cellx = float(info[0])
            celly = float(info[1])
            cellz = float(info[2])

            # write one outputfile for each face
            for idx_water in range(nwater):
                # make {x0,y0,z0} our origin
                # the factor 10 converts to Angstrom from nm
                x = 10.0*xlist[idx_water] - x0
                y = 10.0*ylist[idx_water] - y0
                z = 10.0*zlist[idx_water] - z0

                # --- FACE +X ---
                # is water above +x surface of cube?
                if (x > l/2)              and \
                   (y > -l/2 and y < l/2) and \
                   (z > -l/2 and z < l/2)     :
                    height = abs(x-l/2)
                    outfile_px.write(str(height) + "\n")
            
                # --- FACE -X ---
                # is water below -x surface of cube?
                if (x < -l/2)             and \
                   (y > -l/2 and y < l/2) and \
                   (z > -l/2 and z < l/2)     :
                    height = abs(x+l/2)
                    outfile_mx.write(str(height) + "\n")
            
                # --- FACE +Y ---
                # is water above +y surface of cube?
                if (y > l/2)              and \
                   (x > -l/2 and x < l/2) and \
                   (z > -l/2 and z < l/2)     :
                    height = abs(y-l/2)
                    outfile_py.write(str(height) + "\n")
            
                # --- FACE -Y ---
                # is water below -y surface of cube?
                if (y < -l/2)             and \
                   (x > -l/2 and x < l/2) and \
                   (z > -l/2 and z < l/2)     :
                    height = abs(y+l/2)
                    outfile_my.write(str(height) + "\n")
            
                # --- FACE +Z ---
                # is water above +x surface of cube?
                if (z > l/2)              and \
                   (y > -l/2 and y < l/2) and \
                   (x > -l/2 and x < l/2)     :
                    height = abs(z-l/2)
                    outfile_pz.write(str(height) + "\n")
          
                # --- FACE -Z ---
                # is water above +x surface of cube?
                if (z > -l/2)             and \
                   (y > -l/2 and y < l/2) and \
                   (x > -l/2 and x < l/2)     :
                    height = abs(z+l/2)
                    outfile_mz.write(str(height) + "\n")

        trajfile.close()
        outfile_px.close()
        outfile_mx.close()
        outfile_py.close()
        outfile_my.close()
        outfile_pz.close()
        outfile_mz.close()

##########################################################################################
# (2)                        GET LATERAL POSITION OVER CUBE                              #
##########################################################################################
def lateral_profile_cube( trajfilename,  \
                          outfilename,   \
                          nframes,       \
                          natoms,        \
                          nwater,        \
                          x0,            \
                          y0,            \
                          z0,            \
                          l=20.0,        \
                          height_lo=0.0, \
                          height_hi=4.0) :

    """ 
    [ADAPTED FROM STEVE]

    Gives the lateral positions of the water molecules provided that
    they reside a particular height from the cube

    six outputfiles are created, one for each surface (px, mx, py, my, pz, mz)

    trajfilename     = .gro trajectory file
    outfilename      = output file - list of lateral positons in the plane of the GNF
    nframes          = self explanatory
    natoms           = total number of atoms in traj
    nwater           = self explanatory
    x0,y0,z0         = centre of the cube
    l                = defines the edge length of the cube
    zcut_lo, zcut_hi = defines the liquid layer of interest

    """
    
    if os.path.isfile(outfilename+"_px.dat") or \
       os.path.isfile(outfilename+"_mx.dat") or \
       os.path.isfile(outfilename+"_py.dat") or \
       os.path.isfile(outfilename+"_my.dat") or \
       os.path.isfile(outfilename+"_pz.dat") or \
       os.path.isfile(outfilename+"_mz.dat") :
        # don't allow overwriting of files
        print "Not allowed to overwrite files!"
        print "attempted file to overwrite: %s" % outfilename
        sys.exit(1)
    else:
        trajfile = open(trajfilename,'r')

        # create six outfiles, one for each cube-face
        # names: +x (px), -x (mx), +y (py), -y (my), ..
        outfile_px  = open(outfilename+"_px.dat",'w')
        outfile_mx  = open(outfilename+"_mx.dat",'w')
        outfile_py  = open(outfilename+"_py.dat",'w')
        outfile_my  = open(outfilename+"_my.dat",'w')
        outfile_pz  = open(outfilename+"_pz.dat",'w')
        outfile_mz  = open(outfilename+"_mz.dat",'w')

        for frame in range(nframes):
            # read the two header lines in the .gro file
            line = trajfile.readline()
            line = trajfile.readline()

            # construct empty lists containing the atom positions 
            xlist = []
            ylist = []
            zlist = []

            # now loop through each water
            for i in range(natoms):

                line = trajfile.readline()
                info = line.split()

                # sometimes (n>9999) 2nd and 3rd column not separated in gro
                if len(info) == 6:
                  if info[1][:2] == "mW":
                    x = float(info[3])
                    y = float(info[4])
                    z = float(info[5])
                  else:
                    continue
                elif len(info) == 5:
                  if info[1][:2] == "mW":
                    x = float(info[2])
                    y = float(info[3])
                    z = float(info[4])
                  else:
                    continue
                else:
                  print "[ERROR]. Can't read gro file,", trajfilename
                  print "line: ", line_number
                  print line
                  sys.exit()

                xlist.append(x)
                ylist.append(y)
                zlist.append(z)
            
            # check if all water molecules were read in
            if len(xlist) != len(ylist) != len(zlist) != nwater:
              print "[ERROR]. Could not read the specified number of water molecules"
              print len(xlist), len(ylist), len(zlist)
              sys.exit()

            # read in the cell-vectors
            line = trajfile.readline()
            info = line.split()
            cellx = float(info[0])
            celly = float(info[1])
            cellz = float(info[2])

            # write one outputfile for each face
            for idx_water in range(nwater):
                # make {x0,y0,z0} our origin
                # the factor 10 converts to Angstrom from nm
                x = 10.0*xlist[idx_water] - x0
                y = 10.0*ylist[idx_water] - y0
                z = 10.0*zlist[idx_water] - z0

                # --- FACE +X ---
                # is water in specified height interval surface of cube +x?
                if (x >= l/2+height_lo and x <= l/2+height_hi) and \
                   (y >= -l/2 and y <= l/2) and \
                   (z >= -l/2 and z <= l/2)     :
                    outfile_px.write(str(y) + " " + str(z) + "\n")
            
                # --- FACE -X ---
                if (x <= -l/2-height_lo and x >= -l/2-height_hi) and \
                   (y >= -l/2 and y <= l/2) and \
                   (z >= -l/2 and z <= l/2)     :
                    outfile_mx.write(str(y) + " " + str(z) + "\n")
            
                # --- FACE +Y ---
                if (y >= l/2+height_lo and y <= l/2+height_hi) and \
                   (x >= -l/2 and x <= l/2) and \
                   (z >= -l/2 and z <= l/2)     :
                    outfile_py.write(str(x) + " " + str(z) + "\n")
            
                # --- FACE -Y ---
                if (y <= -l/2-height_lo and y >= -l/2-height_hi) and \
                   (x >= -l/2 and x <= l/2) and \
                   (z >= -l/2 and z <= l/2)     :
                    outfile_my.write(str(x) + " " + str(z) + "\n")
            
                # --- FACE +Z ---
                if (z >= l/2+height_lo and z <= l/2+height_hi) and \
                   (y >= -l/2 and y <= l/2) and \
                   (x >= -l/2 and x <= l/2)     :
                    outfile_pz.write(str(y) + " " + str(x) + "\n")
          
                # --- FACE -Z ---
                if (z <= -l/2-height_lo and z >= -l/2-height_hi) and \
                   (y >= -l/2 and y <= l/2) and \
                   (x >= -l/2 and x <= l/2)     :
                    outfile_mz.write(str(y) + " " + str(x) + "\n")

        trajfile.close()
        outfile_px.close()
        outfile_mx.close()
        outfile_py.close()
        outfile_my.close()
        outfile_pz.close()
        outfile_mz.close()

##########################################################################################
# (3)                        HEIGHT OVER SLAB (INCL. ALL-ATOM)                           # 
##########################################################################################  
def height_from_slab( trajfilename, \
                      outfilename,  \
                      nframes,      \
                      natoms,       \
                      nwater,       \
                      slab_height = 13.8):
    """ 
    [ADAPTED FROM STEVE]

    Calculates the height of a water molecule above a surface slab
    all-atom trajectory is evaluated

    trajfilename = .gro trajectory file, water: type: SOL, atomname: OW
    outfilename  = output file - list of heights above SF
    nframes      = self explanatory
    natoms       = number of atoms in trajectory (don't forget dummy sites)
    nwater       = number of water molecules (=number of OW) in trajectory
    slab_height  = height of slab in Angstrom -->> heights of water molecules relative to that
    """
    
    if os.path.isfile(outfilename+".dat") :
        # don't allow overwriting of files
        print "Not allowed to overwrite files!"
        print "attempted file to overwrite: %s.dat" % outfilename
        sys.exit(1)
    else:
        trajfile = open(trajfilename,'r')

        # create outfile
        outfile  = open(outfilename+".dat",'w')

        line_number = 1    
        for frame in range(nframes):
            # read the two header lines in the .gro file
            line = trajfile.readline()
            line = trajfile.readline()

            line_number += 2

            # construct empty lists containing the atom positions 
            xlist = []
            ylist = []
            zlist = []

            # now loop through each water
            for i in range(natoms):

                line = trajfile.readline()
                line_number += 1
                info = line.split()

                # sometimes (n>9999) 2nd and 3rd column not separated in gro
                if len(info) == 6:
                  if info[1][:2] == "OW":
                    x = float(info[3])
                    y = float(info[4])
                    z = float(info[5])
                  else:
                    continue
                elif len(info) == 5:
                  if info[1][:2] == "OW":
                    x = float(info[2])
                    y = float(info[3])
                    z = float(info[4])
                  else:
                    continue
                else:
                  print "[ERROR]. Can't read gro file,", trajfilename
                  print "line: ", line_number
                  print line
                  sys.exit()

                xlist.append(x)
                ylist.append(y)
                zlist.append(z)

            # check if all water molecules were read in
            if len(xlist) != len(ylist) != len(zlist) != nwater:
              print "[ERROR]. Could not read the specified number of water molecules"
              print len(xlist), len(ylist), len(zlist)
              sys.exit()

            # read in the cell-vectors
            line = trajfile.readline()
            info = line.split()
            cellx = float(info[0])
            celly = float(info[1])
            cellz = float(info[2])
  
            # get height of water molecule (OW)
            for idx_water in range(nwater):
              # convert height from nm to Angstrom
              z = 10.0*zlist[idx_water]

              # calculate height relative to slab-height
              height = abs(z-slab_height)

              outfile.write(str(height) + "\n")
            
        trajfile.close()
        outfile.close()

##########################################################################################
# (4)                        GET LATERAL POSITION OVER SLAB (ALL-ATOM)                   #
##########################################################################################
def lateral_profile_slab( trajfilename,            \
                          outfilename,             \
                          nframes,                 \
                          natoms,                  \
                          nwater,                  \
                          slab_height=13.8,        \
                          height_lo=0.0,           \
                          height_hi=4.0,           \
                          analyt=0,              ) :

    """ 

    Gives the lateral positions of the water molecules over slab
    for all-atom simulation
    for analysis --> oxygen ("OW") position used
    

    trajfilename     = .gro trajectory file
    outfilename      = output file
    nframes          = self explanatory
    natoms           = total number of atoms in traj
    nwater           = self explanatory
    slab_height      = defines the height of the slab
    zcut_lo, zcut_hi = defines the liquid layer of interest
    analyt           = look at water ("OW")               ... analyt = 0
                               hydroxy-groups ("oh, ohs") ... analyt = 1
                               only "oh"                  ... analyt = 2
                               only "ohs"                 ... analyt = 3

    """
    
    if os.path.isfile(outfilename+".dat") :
        # don't allow overwriting of files
        print "Not allowed to overwrite files!"
        print "attempted file to overwrite: %s.dat" % outfilename
        sys.exit(1)
    else:
        trajfile = open(trajfilename,'r')

        # create outfile
        outfile  = open(outfilename+".dat",'w')

        for frame in range(nframes):
            # read the two header lines in the .gro file
            line = trajfile.readline()
            line = trajfile.readline()

            # construct empty lists containing the atom positions 
            xlist = []
            ylist = []
            zlist = []

            # now loop through each water
            for i in range(natoms):

                line = trajfile.readline()
                info = line.split()

                # sometimes (n>9999) 2nd and 3rd column not separated in gro
                if len(info) == 6:
                  # --- analyse water ---
                  if (analyt == 0):
                    if info[1][:2] == "OW":
                      x = float(info[3])
                      y = float(info[4])
                      z = float(info[5])
                    else:
                        continue
                  # --- analyse hydroxy groups (oh + ohs) ---
                  elif (analyt == 1):
                    # ohs contains oh in its first two letters
                    if info[1][:2] == "oh":
                      x = float(info[2])
                      y = float(info[3])
                      z = float(info[4])
                    else:
                        continue
                  # --- analyse hydroxy groups (oh) ---
                  elif (analyt == 2):
                    if info[1][:2] == "oh" and info[1][:3] != "ohs":
                      x = float(info[3])
                      y = float(info[4])
                      z = float(info[5])
                    else:
                        continue
                  # --- analyse hydroxy groups (ohs) ---
                  elif (analyt == 3):
                    if info[1][:3] == "ohs":
                      x = float(info[3])
                      y = float(info[4])
                      z = float(info[5])
                    else:
                        continue
                elif len(info) == 5:
                  # --- analyse water ---
                  if (analyt == 0):
                    if info[1][:2] == "OW":
                      x = float(info[2])
                      y = float(info[3])
                      z = float(info[4])
                    else:
                      continue
                  # --- analyse hydroxy groups (oh + ohs) ---
                  elif (analyt == 1):
                    # ohs contains oh in its first two letters
                    if info[1][:2] == "oh":
                      x = float(info[2])
                      y = float(info[3])
                      z = float(info[4])
                    else:
                        continue
                  # --- analyse hydroxy groups (oh) ---
                  elif (analyt == 2):
                    if info[1][:2] == "oh" and info[1][:3] != "ohs":
                      x = float(info[2])
                      y = float(info[3])
                      z = float(info[4])
                    else:
                        continue
                  # --- analyse hydroxy groups (ohs) ---
                  elif (analyt == 3):
                    if info[1][:3] == "ohs":
                      x = float(info[2])
                      y = float(info[3])
                      z = float(info[4])
                    else:
                        continue
                else:
                  print "[ERROR]. Can't read gro file,", trajfilename
                  print "line: ", line_number
                  print line
                  sys.exit()

                xlist.append(x)
                ylist.append(y)
                zlist.append(z)
            
            # check if all molecules were read in
            if (len(xlist) != nwater) or \
               (len(ylist) != nwater) or \
               (len(zlist) != nwater):
              print "[ERROR]. Could not read the specified number of water molecules (", nwater, ")"
              print len(xlist), len(ylist), len(zlist)
              sys.exit()

            # read in the cell-vectors
            line = trajfile.readline()
            info = line.split()

            cell_a = []
            cell_b = []
            cell_c = []
            
            # orthorombic cell
            if len(info) == 3:
              cell_a = [ float(info[0]) , 0.0 , 0.0 ]
              cell_b = [ 0.0 , float(info[1]) , 0.0 ]
              cell_c = [ float(info[2]) ]
            # triclinic cell
            elif len(info) == 9:
              cell_a = [ float(info[0]), float(info[3]), float(info[4]) ]
              cell_b = [ float(info[5]), float(info[1]), float(info[6]) ]
              cell_c = [ float(info[7]), float(info[8]), float(info[2]) ]
            else:
              print "[ERROR]. The lattice vector could not be read"
              print line
              sys.exit()


            # write outputfile
            for idx_water in range(nwater):
                # --- if necessary: change xy from rhomboid to rectangle ---

                # this procedure does not have any effect if cell is rectangular to start with
                # for every value outside rectangle -> transform coord into rec
                if x < 0:
                  x += cell_a[0]
                if x > cell_a[0]:
                  x -= cell_a[0]
                if y < 0:
                  y += cell_b[1]
                if y > cell_b[1]:
                  y -= cell_b[1]

                # the center of the rectangle: a/2 , b/2
                x = 10.0*(xlist[idx_water] - cell_a[0]/2.0)
                y = 10.0*(ylist[idx_water] - cell_b[1]/2.0)
                z = 10.0*zlist[idx_water] -  slab_height
                
                if (z>=height_lo) and (z<height_hi):
                  outfile.write(str(x) + " " + str(y) + "\n")
          
        trajfile.close()
        outfile.close()

##########################################################################################
# (5)                        INTEGRATE HEIGHT PROFILE (CUBE)                             # 
##########################################################################################  

##########################################################################################
# (6)                        INTEGRATE HEIGHT PROFILE (SLAB)                             # 
##########################################################################################  
def integrate_height_profile( trajfilename, \
                              outfilename,  \
                              nframes,      \
                              natoms,       \
                              nwater,       \
                              slab_height,  \
                              h_low,        \
                              h_high       ):
    """ 

    Integrate the water molecules in a specified layer 

    trajfilename = .gro trajectory file
    outfilename  = outputfile, containing n_water in specified peak for every frame
    nframes      = self explanatory (necessary for average value and standard deviation
    natoms       = number of atoms in trajectory (don't forget dummy sites)
    nwater       = number of water molecules (=number of OW) in trajectory
    slab_height  = slab height to which the h is calculated relatively (in Angstrom)
    h_low        = layer start (in Angstrom)
    h_high       = layer end   (in Angstrom)
    """
    
    if os.path.isfile(outfilename+".dat") :
        # don't allow overwriting of files
        print "Not allowed to overwrite files!"
        print "attempted file to overwrite: %s.dat" % outfilename
        sys.exit(1)
    else:
        trajfile = open(trajfilename,'r')

        # create outfile
        outfile  = open(outfilename+".dat",'w')

        line_number = 1    
        for frame in range(nframes):
            # read the two header lines in the .gro file
            line = trajfile.readline()
            line = trajfile.readline()

            line_number += 2

            # construct empty lists containing the atom positions 
            xlist = []
            ylist = []
            zlist = []
            water_in_layer = 0

            # now loop through each water
            for i in range(natoms):

                line = trajfile.readline()
                line_number += 1
                info = line.split()

                # sometimes (n>9999) 2nd and 3rd column not separated in gro
                if len(info) == 6:
                  if info[1][:2] == "OW":
                    x = float(info[3])
                    y = float(info[4])
                    z = float(info[5])
                  else:
                    continue
                elif len(info) == 5:
                  if info[1][:2] == "OW":
                    x = float(info[2])
                    y = float(info[3])
                    z = float(info[4])
                  else:
                    continue
                else:
                  print "[ERROR]. Can't read gro file,", trajfilename
                  print "line: ", line_number
                  print line
                  sys.exit()

                xlist.append(x)
                ylist.append(y)
                zlist.append(z)

            # check if all water molecules were read in
            if len(xlist) != len(ylist) != len(zlist) != nwater:
              print "[ERROR]. Could not read the specified number of water molecules"
              print len(xlist), len(ylist), len(zlist)
              sys.exit()

            # read in the cell-vectors
            line = trajfile.readline()
            info = line.split()
            cellx = float(info[0])
            celly = float(info[1])
            cellz = float(info[2])
  
            # get height of water molecule (OW)
            for idx_water in range(nwater):
              # convert height from nm to Angstrom
              z = 10.0*zlist[idx_water]

              # calculate height relative to slab-height
              height = abs(z-slab_height)

              # check if water is between h_low and h_high
              if height >= h_low and height < h_high:
                water_in_layer += 1

            # write water in layer for this frame into outfile
            outfile.write(str(water_in_layer) + "\n")
            
        trajfile.close()
        outfile.close()

##########################################################################################
# (7)                        CALCULATE 2-D PROBABILITY DISTRIBUTION                      # 
##########################################################################################  
def Lateral_Prob( infilename,  \
                  outfilename, \
                  xmin,        \
                  xmax,        \
                  ymin,        \
                  ymax,        \
                  dx,          \
                  dy):
    """ 
    [CODE FROM STEVE]

    infilename  : name of input file. Two column whitespace
                  separated. Accepts "#" as comment.
    outfilename : name of output file, given as x, y, z data.
    xmin, ymin  : specifies lower limit of distribution.
    xmax, ymax  : specifies upper limit of distribution.
    dx, dy      : bin widths.

    """

    if os.path.isfile(outfilename):
        # don't allow overwriting of files
        print "Not allowed to overwrite files!"
        print "attempted file to overwrite: %s" % outfilename
        sys.exit(1)
    else:
        infile  = open(infilename ,'r')
        outfile = open(outfilename,'w')
        
        nbin_x = math.trunc((xmax-xmin)/dx) + 1
        nbin_y = math.trunc((ymax-ymin)/dy) + 1

        prob = np.empty( (nbin_x, nbin_y) )
        prob.fill(10**-10)
        
        ndat = 0.0

        outfile.write("# xmin = " + str(xmin) + ", ymin = " + str(ymin) + 
                      ", xmax = " + str(xmax) + ", ymax = " + str(ymax) +
                      ", dx = "   + str(dx)   + ", dy = "   + str(dy) + " Angstrom\n")

        # Loop over all lines in the input file
        for line in infile:
            info = line.split()

            # check for comment line
            if info[0][0] == "#":
                continue
            
            # increase the total number of data points by 1
            ndat += 1.0

            x = float(info[0])
            y = float(info[1])

            # calculate the appropriate bin
            xbin = math.trunc( (x-xmin)/dx )
            ybin = math.trunc( (y-ymin)/dy )


            # check xbin/ybin
            # if the xy plane is not perfectly dividable into the grid specified,
            # errors might occur (index out of bounds)
            if xbin == nbin_x: 
              xbin-= 1
            if ybin == nbin_y: 
              ybin-= 1

            if xbin > nbin_x-1:
              print "[ERROR]. xbin out of bounds."
              print xbin, nbin_x
              sys.exit()
            if ybin > nbin_y-1:
              print "[ERROR]. ybin out of bounds."
              print ybin, nbin_y
              sys.exit()

            prob[xbin][ybin] += 1.0
        
        # normalize the distribution
        prob = prob/ndat

        logP = np.empty( (nbin_x, nbin_y) )
        logP = (-1)*np.log(prob)

        zero_logP = np.amin(logP)

        # print the probability distribution to the outfile
        for i in range(nbin_x):

            bc_x = xmin + i*dx

            for j in range(nbin_y):

                bc_y = ymin + j*dy

                logP[i][j] -= zero_logP 

                outfile.write(str(bc_x) + " " + str(bc_y) + " " + str(prob[i][j]) + " " + str(logP[i][j]) + "\n")

            outfile.write("\n")

        outfile.close()
        infile.close()

##########################################################################################
# (8)                        FREE ENERGY DISTRIBUTION                                    # 
##########################################################################################  
def Free_Energy( infilename,  \
                 outfilename, \
                 T):
    """
    calculate the free energy from probability distribution

    F = kT *ln(P)
     with P ... probability distribution
     ==> with this we get only some part of the free energy of the system, namely the configurational xy free energy
     ==> free energy cost associated with displacing from minimum in eg lateral probability profile

    units for F here: kcal/(mol)

    k_b = 0.0019872041 kcal/(mol*K)

    infilename  ... infile (eg 2D lateral probability dist)
    outfilename ... outfile
    T           ... temperature of simulation  
    """

    if os.path.isfile(outfilename):
        # don't allow overwriting of files
        print "Not allowed to overwrite files!"
        print "attempted file to overwrite: %s" % outfilename
        sys.exit(1)
    else:
        infile  = open(infilename ,'r')
        outfile = open(outfilename,'w')
        
        outfile.write("# free energy profile from " + infilename + " ; temperature = " + str(T) + "\n")

        k_b = 0.0019872041
                      
        # Loop over all lines in the input file
        for line in infile:
            info = line.split()

            # check for empty line --> needed for gnuplot
            if len(info) == 0:
              outfile.write("\n")
              continue

            # check for comment line
            if info[0][0] == "#":
                continue

            free_energy = k_b * T * float(info[3])
            outfile.write( info[0] + " " + info[1] +  " " + str(free_energy) + "\n")

    outfile.close()
    infile.close()
