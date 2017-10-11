import os.path
import math
import sys
import numpy as np
import re   # strip whitespaces (used for reading GRO format)

"""

different file formats are defined here as classes

supported at the moment:
  (+) init   (LAMMPS)
  (+) xyz    (commonly used)
  (+) POSCAR (VASP)
  (+) gro    (GROMACS)
  (+) nx4a   (proton - disordered ice following ice rules)
  (+) pdb    (very crude version!!!)

"""

##########################################################################################
#                            INIT (LAMMPS)                                               # 
##########################################################################################  
class format_init:
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            __init__                                                   #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def __init__(self):
    self.titel       = "no titel specified" # titel
    self.n_atoms     = 0                    # number of atoms
    self.n_atomtypes = 0                    # number of atomtypes
    self.masses      = []                   # [atomtype, mass]
    self.box_x       = []                   # x_start, x_end
    self.box_y       = []                   # y_start, y_end
    self.box_z       = []                   # z_start, z_end
    self.box_tilt    = []                   # yx, xz, yz
    self.list_atoms  = []                   # [atomnumber, atomtype, x, y, z, proj_x, proj_y, proj_z]
  
  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            check_data()                                               #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def check_data(self):
    """
    check if the data in class makes sense

    returns 0 if all the data is OK
    returns 1 if there is a problem with the data
    """

    if self.n_atoms == 0:
      print "[ERROR]. n_atoms == 0"
      sys.exit()
    if self.n_atomtypes == 0:
      print "[ERROR]. n_atomtypes == 0"
      sys.exit()
    if len(self.masses) != self.n_atomtypes:
      print "[ERROR]. len(masses) != n_atomtypes"
      sys.exit()
    if len(self.list_atoms) != self.n_atoms:
      print "[ERROR]. len(list_atoms) != n_atoms"
      print len(self.list_atoms)
      print self.n_atoms
      sys.exit()
    if len(self.box_x) != 2:
      print "[ERROR]. x box dimensions missing"
      sys.exit()
    if len(self.box_y) != 2:
      print "[ERROR]. y box dimensions missing"
      sys.exit()
    if len(self.box_z) != 2:
      print "[ERROR]. z box dimensions missing"
      sys.exit()
    if len(self.box_tilt) != 0:
      if len(self.box_tilt) != 3:
        print "[ERROR]. box tilt wrong"
        sys.exit()


  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            write_file                                                 #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def write_file(self, outfilename):
    """
    write outputfile
    """
    outfile  = open(outfilename, 'w')

    # check if all parameters are there
    flag_data_ok = self.check_data()

    # write generic info
    outfile.write(self.titel+"\n\n")
    outfile.write(str(self.n_atoms) + " atoms\n\n")
    outfile.write(str(self.n_atomtypes) + " atom types\n\n")

    # write box
    outfile.write(str(self.box_x[0]) + " "  + str(self.box_x[1]) + " xlo xhi\n")
    outfile.write(str(self.box_y[0]) + " "  + str(self.box_y[1]) + " ylo yhi\n")
    outfile.write(str(self.box_z[0]) + " "  + str(self.box_z[1]) + " zlo zhi\n")

    if len(self.box_tilt) != 0:
      # in LAMMPS, the tilt factor can not skew the box more than half the box length
      # eg: xlo=2, xhi=12 --> x-length=10 --> tilt factor needs to be in between -5 and 5
      # in this example, configurations with tilt -15, -5, 5, 15, ... ARE EQUIVALENT
      #   if this happens -> automatically corrected here
      x_len = self.box_x[1] - self.box_x[0]
      self.box_tilt[0]=self.box_tilt[0] - int(self.box_tilt[0]/x_len)*x_len

      if self.box_tilt[0] > x_len/2.0:
        self.box_tilt[0] -= x_len
      elif self.box_tilt[0] < -x_len/2.0:
        self.box_tilt[0] += x_len

      outfile.write(str(self.box_tilt[0]) + " " + str(self.box_tilt[1]) + " " +  str(self.box_tilt[2]) + " xy xz yz\n")

    outfile.write("\n")

    # write masses
    outfile.write("Masses\n\n")
    for i in self.masses:
      outfile.write(str(i[0]) + " " + str(i[1]) + "\n")

    outfile.write("\n")
    
    # write atoms
    outfile.write("Atoms\n\n")
    for i in self.list_atoms:
      outfile.write( str(i[0]) + " " + \
                     str(i[1]) + " " + \
                     str(i[2]) + " " + \
                     str(i[3]) + " " + \
                     str(i[4]) + " " + \
                     str(i[5]) + " " + \
                     str(i[6]) + " " + \
                     str(i[7]) + "\n")
    # close file  
    outfile.close()

##########################################################################################
#                            PDB (COMMON FORMAT)                                         # 
##########################################################################################  
class format_pdb:
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            __init__                                                   #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def __init__(self):
    self.n_atoms     = 0                    # number of atoms
    self.list_atoms  = []                   # [ [atomname, x, y, z], ... ] for each frame
    self.lattice_vector_1  = []             # line 1
    self.lattice_vector_2  = []             # line 1
    self.lattice_vector_3  = []             # line 1

    # makros
    self.ATOM_NAME = 0
    self.X_COORD   = 1
    self.Y_COORD   = 2
    self.Z_COORD   = 3

  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            read_file                                                  #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def read_file(self, infilename):
    """
    read infile
    """

    n_line = 0
    with open(infilename, "r") as f:
      for line in f:
        data = line.split()
          
        # read pbc
        if n_line == 0:
          # ignore element [0]
          len_a = float(data[1])
          len_b = float(data[2])
          len_c = float(data[3])
          alpha = float(data[4]) # should be 90 in slabs 
          beta  = float(data[5]) # should be 90 in slabs
          gamma = float(data[6]) # angle between a and b

          if (alpha != 90.0) or (beta != 90.0):
            print "[ERROR]. file-formats: read_xyz()"
            print "         In slab models, alpha and beta should both be 90 degree!", alpha, beta
            sys.exit()

          self.lattice_vector_1 = [ len_a, 0.0, 0.0 ]
          self.lattice_vector_2 = [ round(len_b*np.cos(gamma*np.pi/180.0),10), round(len_b*np.sin(gamma*np.pi/180.0),10), 0.0 ]
          self.lattice_vector_3 = [0.0, 0.0, len_c]

        # reached EOF
        elif data[0] == "END":
          break # done reading the file

        # read atom info
        else:
          self.list_atoms.append([ data[2], float(data[6]), float(data[7]), float(data[8])  ])

        n_line += 1

      f.close()
      self.n_atoms = len(self.list_atoms)

##########################################################################################
#                            XYZ (COMMON FORMAT)                                         # 
##########################################################################################  
class format_xyz:
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            __init__                                                   #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def __init__(self):
    self.titel       = "no titel specified" # titel
    self.n_atoms     = 0                    # number of atoms in 1 frame
    self.list_atoms  = []                   # [ [atomname, x, y, z], ... ] for each frame
    self.n_frames    = 0                    # number of frames (<0 means manual detection)
    self.lattice_vector_1  = []             # set manually by user! (or specify via pbc set in title)
    self.lattice_vector_2  = []             # set manually by user! (or specify via pbc set in title)
    self.lattice_vector_3  = []             # set manually by user! (or specify via pbc set in title)

    # makros
    self.ATOM_NAME = 0
    self.X_COORD   = 1
    self.Y_COORD   = 2
    self.Z_COORD   = 3

  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            read_file                                                  #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def read_file(self, infilename, n_frames=0, n_atoms=0, flag_read_extra_info=False, flag_read_UC_dim=False):
    """
    read infile

    n_frames ... if set to 0, then n_frames will be automatically determined on the run

    the xyz file might contain extra info (written by e.g. selector class)
    this is useful for further processing with analyser class

    to see the format of this --> consult write_file method

    flag_read_UC_dim ... read unit cell dimension from xyz file
                         if set to true, the unit cell has to be specified as pbc set { A B C alpha beta gamma }
                         gamma --> angle between A and B
                         alpha, beta --> for now should be 0, because it is angle to c (slab)
                         possible formats:
                            pbc set { A B C alpha beta gamma }    --> split() will result in 10 elements
                            A B C alpha beta gamma                --> split() will result in 6 elements
    """

    # if number of atoms not specified --> get it from first line in *xyz
    n_atoms_file=0
    with open(infilename, "r") as f:
      line_1 = f.readline()
      n_atoms_file = int(line_1)
    
    if (n_atoms > 0) and (n_atoms != n_atoms_file):
      print "[ERROR]. file formats - read xyz"
      print "         The number of atoms in the file does not match argument"
      print "           n_atoms", n_atoms
      print "           n_atoms_file", n_atoms_file
      sys.exit()

    n_atoms = n_atoms_file

    # calculate n_frames from file if n_frames<0
    # this is useful if analysis is performed on set of trajectories with different n_frames

    # get number of lines
    n_lines = sum(1 for line in open(infilename, 'r'))
    
    if n_frames > 0:
      if n_lines != ( n_frames*(n_atoms+2) ):
        print "[ERROR]. The specified n_frames and n_atoms does not match up with the number of lines."
        print "         n_lines found   : ", n_lines
        print "         n_lines expected: ", n_frames*(n_atoms+2)
        sys.exit()
    else:
      if n_lines % (n_atoms+2) != 0:
        print "[ERROR]. Cannot determine n_frames, n_atoms is not compatible with any interger number of frames"
        sys.exit()

      n_frames = n_lines/(n_atoms+2)
        
    infile  = open(infilename, 'r')
 
    self.n_frames = n_frames # set n_frames

    for frame in range(n_frames):
      # monitor progress
      #if int(n_frames*0.25)==frame:
      #  print "         25 % done"
      #if int(n_frames*0.50)==frame:
      #  print "         50 % done"
      #if int(n_frames*0.75)==frame:
      #  print "         75 % done   "

      # line 1 ... number of atoms
      line = infile.readline()
      data = line.split()
      self.n_atoms = int(data[0])
      
      # line 2 ... titel of frame
      line = infile.readline()
      data = line.split()
      self.titel = line


      # read unit cell dimension from title
      if flag_read_UC_dim == True:
        # check if UC dimension contained in title
        list_dim = self.titel.split()
        if len(list_dim)==6:
          len_a = float(list_dim[0])
          len_b = float(list_dim[1])
          len_c = float(list_dim[2])
          alpha = float(list_dim[3]) # should be 90 in slabs 
          beta  = float(list_dim[4]) # should be 90 in slabs
          gamma = float(list_dim[5]) # angle between a and b

          self.lattice_vector_1 = [ len_a, 0.0, 0.0 ]
          self.lattice_vector_2 = [ round(len_b*np.cos(gamma*np.pi/180.0),10), round(len_b*np.sin(gamma*np.pi/180.0),10), 0.0 ]
          self.lattice_vector_3 = [0.0, 0.0, len_c]

          if (alpha != 90.0) or (beta != 90.0):
            print "[ERROR]. file-formats: read_xyz()"
            print "         In slab models, alpha and beta should both be 90 degree!", alpha, beta
            sys.exit()

        # --- check if UC info here in frame by frame basis
        else:
          # find "pbc set {" in string
          str_pos = line.find("pbc set {")
          if str_pos < 0:
            print "[ERROR]. file-formats: read_xyz()"
            print "         Could not find 'pbc set {' in title. No unit cell dimension set!"
            sys.exit()
          
          # get values in between { } brackets
          str_pbc = line[str_pos+9:-2]
          data = str_pbc.split()

          len_a = float(data[0])
          len_b = float(data[1])
          len_c = float(data[2])
          alpha = float(data[3]) # should be 90 in slabs 
          beta  = float(data[4]) # should be 90 in slabs
          gamma = float(data[5]) # angle between a and b

          self.lattice_vector_1 = [ len_a, 0.0, 0.0 ]
          self.lattice_vector_2 = [ round(len_b*np.cos(gamma*np.pi/180.0),10), round(len_b*np.sin(gamma*np.pi/180.0),10), 0.0 ]
          self.lattice_vector_3 = [0.0, 0.0, len_c]

          if (alpha != 90.0) or (beta != 90.0):
            print "[ERROR]. file-formats: read_xyz()"
            print "         In slab models, alpha and beta should both be 90 degree!", alpha, beta
            sys.exit()
      
      # atom coordinates in this frame
      frame_atoms = []
      for i in range(n_atoms):
        line = infile.readline()
        data = line.split()
        new_atom = []

        if flag_read_extra_info == False:
          #            name     x               y               z
          new_atom = [ data[0], float(data[1]), float(data[2]), float(data[3])  ]
        else:
          # if dummy atom
          if data[0] == "XX":
            new_atom = [ data[0], float(data[1]), float(data[2]), float(data[3])  ]
          else:
            # read actual xyz data
            #            name     x               y               z
            new_atom = [ data[0], float(data[1]), float(data[2]), float(data[3])  ]

            # read extra info
            extra_info = data[4:]
            
            if len(extra_info[1:]) % 3 != 0:
              print "[ERROR]. The format of extra info is not correct at the moment."
              print "         Please read file_format description"
              sys.exit()

            list_donors = []
            list_acceptors = []

            for idx,partner in enumerate(extra_info[1::3]):
              if partner=="A":
                list_acceptors.append(extra_info[1+idx*3:1+idx*3+3])
              if partner=="D":
                list_donors.append(extra_info[1+idx*3:1+idx*3+3])

            list_extra_info = [ int(extra_info[0]) ] + [list_acceptors] + [list_donors]
          new_atom.extend(list_extra_info)

        frame_atoms.append(new_atom)
        
      self.list_atoms.append(frame_atoms)

    infile.close()

  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            write_file                                                  #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def write_file(self, outfilename, flag_write_extra_info=False):
    """
    write xyz file

    flag_write_extra_info: sometimes it is useful to append extra information to the atoms
                           for example: after some selectors --> keep track of H-bonding network

       the format of this extra info is the following:
       it is first of all just another element in row of atom
       atom[4] ... index of atom in original traj on which bonding was computed
                   this allows to compare traj easily with original traj
       atom[5] ... A or D
                   A ... atom accepts hydrogen bond
                   D ... atom donates hydrogen bond
                   this is followed by exactly two values before new info is written
       atom[6] ... index to O from which hydrogen bond is accepted/donated
                   this index is the same than the index of that particular atom in original traj
       atom[7] ... index of H completing the hydrogen bond
                   in case of A --> H from e.g. hydroxyl group
                   in case of D --> H from H2O that donates the Hbond
       atom[8] ... A/D again and so on (followed by [9], and [10]
    """
    outfile  = open(outfilename, 'w')

    # check if all parameters are there
    #flag_data_ok = self.check_data()

    line=""

    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            one single frame                                           #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    if isinstance(self.list_atoms[0][0], str):
      # if the first element of list_atoms[0] is a string --> atomname --> 1 frame
      outfile.write(str(len(self.list_atoms)) + "\n")
      if self.titel.endswith("\n"):
        outfile.write(self.titel)
      else:
        outfile.write(self.titel+"\n")

      for i in self.list_atoms:
        if flag_write_extra_info==False:
          if len(i) != 4:
            print "[ERROR]. write_file (xyz) failed / single-frame"
            print "         Cannot recognize xyz format, len(atom) != 4"
            sys.exit()
          outfile.write( str(i[0]) + "     " + \
                         format(round(i[1],6), "+.6f") + "   " + \
                         format(round(i[2],6), "+.6f") + "   " + \
                         format(round(i[3],6), "+.6f") + "\n")
        else:
          print "[ERROR]. write_file (xyz) failed"
          print "         Writing extra info not supported for single frame files"
          sys.exit()

    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            multiple frames                                            #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    elif isinstance(self.list_atoms[0][0], list):
      # if the first element of list_atoms[0] is a list --> frames --> trajectory
      for frame in self.list_atoms:
        # write generic info

        outfile.write(str(self.n_atoms) + "\n")
        outfile.write("\n")
        for i in frame:
          if flag_write_extra_info==False:
            if len(i) != 4:
              print "[ERROR]. write_file (xyz) failed / multiple-frames"
              print "         Cannot recognize xyz format, len(atom) != 4"
              print i
              sys.exit()
            outfile.write( str(i[0]) + " " + \
                           str(i[1]) + " " + \
                           str(i[2]) + " " + \
                           str(i[3]) + "\n")
          else:
            if len(i) < 4:
              print "[ERROR]. write_file (xyz) failed"
              print "         not enough data to write extra_info"
            outfile.write( str(i[0]) + " " + \
                           str(i[1]) + " " + \
                           str(i[2]) + " " + \
                           str(i[3]) + " " )
            line= str(i[0]) + " " + \
                  str(i[1]) + " " + \
                  str(i[2]) + " " + \
                  str(i[3]) + " "

            # there are two options here
            #  (1) --> data comes from SELECTOR
            #          in this case all the extra_data will be arranged
            #          as list in i[4]
            #  (2) --> data comes from READING in a file after SELECTOR
            #          processed it already.
            #          in that case all the extra info are arranged as following:
            #          idx, list_acceptors, list_donors
            #
            # ==> this is a mess and was not anticipated, but now too lazy to fix it

            # scenario (1) [OK]
            if isinstance(i[4], list):
              for extra_info in i[4]:
                if isinstance(extra_info, list):
                  for subelement in extra_info:
                    outfile.write(str(subelement) + " ")
                else:
                  outfile.write(str(extra_info) + " ")
            # scenario (2)
            else:
              #print i
              #line += str(i[4]) + " "  #idx
              outfile.write(str(i[4]) + " ")

              if len(i[5:])==2:
                # list_acceptors
                for acceptor in i[5]:
                  if len(acceptor)==3:
                    #line += acceptor[0] + " " + acceptor[1] + " " + acceptor[2] + " "
                    outfile.write(acceptor[0] + " " + acceptor[1] + " " + acceptor[2] + " ")

                # list_donors
                for donor in i[6]:
                  if len(donor)==3:
                    #line += donor[0] + " " + donor[1] + " " + donor[2] + " "
                    outfile.write(donor[0] + " " + donor[1] + " " + donor[2] + " ")
              elif len(i[5:])!=2 and len(i[5:])>0:
                print "[ERROR]. Unknown format..."
                print i
                sys.exit()

            outfile.write("\n")

        # if this frame contains less atoms than the max number of atoms in a frame:
        #  add XX at 0.0 0.0 0.0 so trajectories can be viewed in vmd
        if len(frame) < self.n_atoms:
          # calculate the amount of extra XX needed
          extra_atoms = self.n_atoms-len(frame)
          for extra in range(extra_atoms):
            outfile.write( "XX 0.0 0.0 0.0\n")
    
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            sanity check                                               #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    else:
      print "[ERROR]. write_file (xyz) failed"
      print "         Cannot recognize xyz format, not single frame / trajectory in list_atoms"
      sys.exit()

    # close file  
    outfile.close()

  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            MERGE TRAJECTORIES                                         #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def merge_trajectories ( self, \
                           list_trajectories=[] ):
    """
    this tool is used to merge different trajectories together

    this is important, if the particle number from one traj to another changes
    (happens for almost all selectors)
    --> especially important if we want to analyse dynamical data

    the result is that self.traj contains the concenated trajectory

    list_trajectories ... list of file_format objects, each object representing one traj
    """

    if len(list_trajectories)==0:
      print "[ERROR]. Cannot merge trajectories, because none are specified."
      sys.exit()

    self.title=""

    # figure out the maximum number of atoms that occur in a trajectory
    # also: the number of total frames
    n_frames = 0
    list_n_atoms = []
    for traj in list_trajectories:
      list_n_atoms.append(int(traj.n_atoms))
      n_frames += traj.n_frames

    print list_n_atoms
    self.n_atoms  = int(max(list_n_atoms))
    self.n_frames = n_frames

    # check if the lattice dimensions are the same for all traj
    for idx,traj in enumerate(list_trajectories[1:]):
      if traj.lattice_vector_1 != list_trajectories[idx-1].lattice_vector_1 or \
         traj.lattice_vector_2 != list_trajectories[idx-1].lattice_vector_2 or \
         traj.lattice_vector_3 != list_trajectories[idx-1].lattice_vector_3:
        print "[ERROR]. Cannot merge trajectories, they have different lattice vectors!"
        print "trajectory", idx-1, " --------- trajectory", idx
        print list_trajectories[idx-1].lattice_vector_1, traj.lattice_vector_1
        print list_trajectories[idx-1].lattice_vector_2, traj.lattice_vector_2
        print list_trajectories[idx-1].lattice_vector_3, traj.lattice_vector_3
        sys.exit()

    # assign lattice_vectors
    self.lattice_vector_1 = list_trajectories[0].lattice_vector_1
    self.lattice_vector_2 = list_trajectories[0].lattice_vector_2
    self.lattice_vector_3 = list_trajectories[0].lattice_vector_3

    # merge trajectories
    for traj in list_trajectories:
      self.list_atoms.extend(traj.list_atoms)

    if len(self.list_atoms) != self.n_frames:
      print "[ERROR]. The number of frames is wrong after merging"
      print "         Expected:", self.n_frames
      print "         Found:   ", len(self.list_atoms)
      sys.exit()


    # fill up every frame with "XX" atoms at 0.0 0.0 0.0
    # so that each of them has exactly n_atoms
    for frame in self.list_atoms:
      if len(frame)<self.n_atoms:
        for i in range(self.n_atoms-len(frame)):
          frame.append(["XX", 0.0, 0.0, 0.0, -1])

    # consistency check: is n_atoms right
    for frame in self.list_atoms:
      if len(frame)<self.n_atoms:
        print "[ERROR]. The number of atoms could not be resolved"
        print "         Expected:", self.n_atoms
        print "         Found:   ", len(frame)
        sys.exit()
        



      
##########################################################################################
#                            POSCAR (VASP)                                               # 
##########################################################################################  
class format_POSCAR:
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            __init__                                                   #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def __init__(self):
    self.titel             = ""
    self.lattice_vector_1  = []
    self.lattice_vector_2  = []
    self.lattice_vector_3  = []
    self.scaling           = 1.0
    self.list_atom_names   = [] # atom names  : H O Au
    self.list_atom_numbers = [] # atom numbers: 2 1 10
    self.coord_format      = 0  # CARTESIAN / FRACTIONAL
    self.CARTESIAN         = 1  # makro
    self.FRACTIONAL        = 2  # makro
    self.list_atoms        = [] # [name, index, X, Y, Z, T(F), T(F), T(F)]

    # makros
    self.ATOMNAME = 0
    self.INDEX    = 1
    self.X_COORD  = 2
    self.Y_COORD  = 3
    self.Z_COORD  = 4

  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            read_file                                                  #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def read_file(self, infilename):
    """
    read POSCAR/CONTCAR
    """

    # open POSCAR
    infile = open(infilename,"r")

    n_line = 0
    atom_names = ""
    flag_read_coords=False
    for line in infile:
      n_line += 1
      complete_line = line
      line = line.split() # split words of line into a list
      
      if   (n_line == 1): self.titel = complete_line
      elif (n_line == 2): self.scaling = float(line[0]) 
      elif (n_line == 3): self.lattice_vector_1 = [ float(line[0]), float(line[1]), float(line[2])]
      elif (n_line == 4): self.lattice_vector_2 = [ float(line[0]), float(line[1]), float(line[2])]
      elif (n_line == 5): self.lattice_vector_3 = [ float(line[0]), float(line[1]), float(line[2])]
      elif (n_line == 6): atom_names = line 
      elif (n_line == 7): 
        self.list_atom_numbers = map(int, line)
    
        index = 0
        # assign atom names
        for name in atom_names:
          for _ in range(self.list_atom_numbers[index]): self.list_atom_names.append(name)
          index += 1
    
        if len(self.list_atom_names) != sum(self.list_atom_numbers):
          print "[ERROR]. Could not construct list_atom_names"
          print "         The atom names and numbers don't match up"
          sys.exit()
      elif (n_line > 7) and flag_read_coords==False:
        line = line[0].lower()
        if (line.startswith("s")): continue #n_line == 8 --> Selective Dynamics (not useful here, just ignore it)
        if (line.startswith("c")): self.coord_format = self.CARTESIAN
        if (line.startswith("d")): self.coord_format = self.FRACTIONAL
        flag_read_coords=True
        # set back index to zero to re-iterate over list of names when reading the coordinates
        index = 0
      elif (n_line > 8) and (flag_read_coords==True) and (index < (sum(self.list_atom_numbers))):  
        new_atom = [ self.list_atom_names[index], index, float(line[0]), float(line[1]), float(line[2]) ]
        index += 1
        self.list_atoms.append(new_atom)
    infile.close()

    # [DEBUG] loop through elements
    # for frame in self.list_atoms:
    #   for atom in frame:
    #     print atom
  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            write_file                                                  #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def write_file(self, outfilename):
    """
    write POSCAR  
    """
    outfile  = open(outfilename, 'w')
    
    # title
    outfile.write(self.titel + "\n")
    
    # scaling
    outfile.write(str(self.scaling) + "\n")
    
    # lattice vectors
    line = "%f %f %f\n" % ( self.lattice_vector_1[0], self.lattice_vector_1[1], self.lattice_vector_1[2])
    outfile.write(line)
    line = "%f %f %f\n" % ( self.lattice_vector_2[0], self.lattice_vector_2[1], self.lattice_vector_2[2])
    outfile.write(line)
    line = "%f %f %f\n" % ( self.lattice_vector_3[0], self.lattice_vector_3[1], self.lattice_vector_3[2])
    outfile.write(line)
    
    # figure out list_atom_names
    self.list_atom_names = []
    for atom in self.list_atoms:
      if atom[0] in self.list_atom_names:
        continue
      else:
        self.list_atom_names.append(atom[0])

    # write list_atoms_names
    for item in self.list_atom_names:
      outfile.write(item + "  ")
    outfile.write("\n")
    
    # figure out atom numbers
    self.list_atom_numbers = []
    for atom_name in self.list_atom_names:
      n = 0
      for atom in self.list_atoms:
        if atom[0] == atom_name:
          n += 1
      self.list_atom_numbers.extend([n])

    # write list_atom_numbers
    for item in self.list_atom_numbers:
      outfile.write(str(item) + "  ")
    outfile.write("\n")
    
    # selective dynamics
    outfile.write("Selective Dynamics\n")
    
    # coordinate representation
    if self.coord_format == self.FRACTIONAL:
      outfile.write("Direct\n")
    elif self.coord_format == self.CARTESIAN:
      outfile.write("Cartesian\n")
    
    # --- atoms ---
    # set up local macros, depending on the format of list_atoms
    # list atoms can either be 
    #  (1) [atom_name, x,y,z]
    #  (2) [atom_name, index, x,y,z]
    #  (3) [atom_name, index, x,y,z, T(F), T(F), T(F)]
    ATOM_NAME = 0
    INDEX     = 1
    X_COORD   = 2
    Y_COORD   = 3
    Z_COORD   = 4
    FIX_X     = 5
    FIX_Y     = 6
    FIX_Z     = 7

    if len(self.list_atoms[0]) == 4:
      ATOM_NAME = 0
      X_COORD   = 1
      Y_COORD   = 2
      Z_COORD   = 3
      FIX_X     = -1   # negative value --> not set   
    elif len(self.list_atoms[0]) == 5:
      ATOM_NAME = 0
      X_COORD   = 2
      Y_COORD   = 3
      Z_COORD   = 4
      FIX_X     = -1   # negative value --> not set   
    elif len(self.list_atoms[0]) == 8:
      # this is the default
      pass
    else:
      print "[ERROR]. file_formats (write_file() in format_POSCAR)"
      print "         list_atoms has a format that is not understood"
      print self.list_atoms[0]
      sys.exit()

    for atom_name in self.list_atom_names:
      for atom in self.list_atoms:
        if atom[ATOM_NAME] == atom_name:
          # T(F) values set? NO
          if FIX_X < 0:
            line = "%f %f %f T T T\n" % ( atom[X_COORD], \
                                          atom[Y_COORD], \
                                          atom[Z_COORD]  )
            outfile.write(line)
          # T(F) values set? YES
          else:
            line = "%f %f %f %s %s %s\n" % ( atom[X_COORD], \
                                             atom[Y_COORD], \
                                             atom[Z_COORD], \
                                             atom[FIX_X], \
                                             atom[FIX_Y], \
                                             atom[FIX_Z] )
            outfile.write(line)
    
    outfile.close()
      
##########################################################################################
#                            GRO (GROMACS)                                               # 
##########################################################################################  
class format_GRO:
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            __init__                                                   #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def __init__(self):
    self.titel             = ""
    self.n_atoms           = 0
    self.list_atoms        = [] # [molecule_number, molecule_name, atomtype, index, X, Y, Z]
    self.lattice_vector_1  = []
    self.lattice_vector_2  = []
    self.lattice_vector_3  = []

    # makros to access list_atoms data (for external use)
    self.MOLECULE_NUMBER = 0
    self.MOLECULE_NAME   = 1
    self.ATOMNAME        = 2
    self.INDEX           = 3
    self.X_COORD         = 4
    self.Y_COORD         = 5
    self.Z_COORD         = 6

  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            read_file                                                  #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def read_file(self, infilename, flag_trajectory=False):
    """
    read grofile

    flag_trajectory: is it a grofile with one frame, or a trajectory?
    """

    # open grofile
    infile = open(infilename,"r")

    n_line = 0

    # *** read gro-file with 1 frame ***
    if flag_trajectory == False:
      for line in infile:
        n_line += 1
        line_complete = line
        line = line.split() # split words of line into a list
        
        if   (n_line == 1): self.titel = line_complete
        elif (n_line == 2): self.n_atoms = int(line[0]) 
        else:
          # read box dimensions
          if len(self.list_atoms) == self.n_atoms:
            if len(line) == 3:
              self.lattice_vector_1 = [ float(line[0]), 0.0, 0.0 ]
              self.lattice_vector_2 = [ 0.0, float(line[1]), 0.0 ]
              self.lattice_vector_3 = [ 0.0, 0.0, float(line[2]) ]
            elif len(line) == 9:
              self.lattice_vector_1 = [ float(line[0]), float(line[3]), float(line[4]) ]
              self.lattice_vector_2 = [ float(line[5]), float(line[1]), float(line[6]) ]
              self.lattice_vector_3 = [ float(line[7]), float(line[8]), float(line[2]) ]
            else:
              print "[ERROR]. The boxvector cannot be read"
              print line
              sys.exit()
          else:
            # remove all whitespaces with library re
            # check length of array line
            # if gro file is poorly formatted (not precicely like intended),
            # then we are still able to read it (based on the lenght of line)
            # if this is not possible, then we follow gromacs rules

            # potentially poorly formatted, but number of elements OK
            if len(line) == 6:
              self.list_atoms.append([ int(line_complete[0:5]),                  \
                                       re.sub(r'\s+', '', line_complete[5:10]),  \
                                       re.sub(r'\s+', '', line_complete[10:15]), \
                                       int(line_complete[15:20]),                \
                                       float(line[3]),                           \
                                       float(line[4]),                           \
                                       float(line[5]) ])

            # follow exact GROMACS rules, if one whitespace is off then meaningless data
            else:
              self.list_atoms.append([ int(line_complete[0:5]),                  \
                                       re.sub(r'\s+', '', line_complete[5:10]),  \
                                       re.sub(r'\s+', '', line_complete[10:15]), \
                                       int(line_complete[15:20]),                \
                                       float(line_complete[20:28]),              \
                                       float(line_complete[28:36]),              \
                                       float(line_complete[36:44]) ])
    # *** read multiple frames ***
    else:
      list_atoms_frame = [] # storing all info about the frames
      
      for line in infile:
        n_line += 1
        line_complete = line
        line = line.split() # split words of line into a list
        
        if   (n_line == 1): self.titel = line_complete
        elif (n_line == 2): self.n_atoms = int(line[0]) 
        else:
          # read box dimensions
          if len(list_atoms_frame) == self.n_atoms:
            if len(line) == 3:
              self.lattice_vector_1 = [ float(line[0]), 0.0, 0.0 ]
              self.lattice_vector_2 = [ 0.0, float(line[1]), 0.0 ]
              self.lattice_vector_3 = [ 0.0, 0.0, float(line[2]) ]
              # prepare for new frame
              n_line = 0
              self.list_atoms.append(list_atoms_frame)
              list_atoms_frame = []

            elif len(line) == 9:
              self.lattice_vector_1 = [ float(line[0]), float(line[3]), float(line[4]) ]
              self.lattice_vector_2 = [ float(line[5]), float(line[1]), float(line[6]) ]
              self.lattice_vector_3 = [ float(line[7]), float(line[8]), float(line[2]) ]
              # prepare for new frame
              n_line = 0
              self.list_atoms.append(list_atoms_frame)
              list_atoms_frame = []

            else:
              print "[ERROR]. The boxvector cannot be read ***"
              print line
              sys.exit()
          else:
            # remove all whitespaces with library re
            # check length of array line
            # if gro file is poorly formatted (not precicely like intended),
            # then we are still able to read it (based on the lenght of line)
            # if this is not possible, then we follow gromacs rules

            # potentially poorly formatted, but number of elements OK
            if len(line) == 6:
              list_atoms_frame.append([ int(line_complete[0:5]),                  \
                                        re.sub(r'\s+', '', line_complete[5:10]),  \
                                        re.sub(r'\s+', '', line_complete[10:15]), \
                                        int(line_complete[15:20]),                \
                                        float(line[3]),                           \
                                        float(line[4]),                           \
                                        float(line[5]) ])
            # follow exact GROMACS rules, if one whitespace is off then meaningless data
            else:
              list_atoms_frame.append([ int(line_complete[0:5]),                  \
                                        re.sub(r'\s+', '', line_complete[5:10]),  \
                                        re.sub(r'\s+', '', line_complete[10:15]), \
                                        int(line_complete[15:20]),                \
                                        float(line_complete[20:28]),              \
                                        float(line_complete[28:36]),              \
                                        float(line_complete[36:44]) ])


  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            write_file                                                  #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def write_file(self, outfilename):
    """
    write grofile
    """
    outfile  = open(outfilename, 'w')

    outfile.write(self.titel + "\n")
    outfile.write(str(self.n_atoms) + "\n")

    for atom in self.list_atoms:
      line = "%5d%-5s%5s%5d%8.3f%8.3f%8.3f\n" % ( atom[0],   \
                                                  atom[1],   \
                                                  atom[2],   \
                                                  atom[3],   \
                                                  atom[4], atom[5], atom[6])
      outfile.write(line)

    line = "%.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf %.10lf\n" % (self.lattice_vector_1[0],   \
                                                                                 self.lattice_vector_2[1],   \
                                                                                 self.lattice_vector_3[2],   \
                                                                                 self.lattice_vector_1[1],   \
                                                                                 self.lattice_vector_1[2],   \
                                                                                 self.lattice_vector_2[0],   \
                                                                                 self.lattice_vector_2[2],   \
                                                                                 self.lattice_vector_3[0],   \
                                                                                 self.lattice_vector_3[1])
    outfile.write(line)
    
    # close file  
    outfile.close()

##########################################################################################
#                            NX4A (proton disordered ice)                                # 
##########################################################################################  
class format_nx4a:
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            __init__                                                   #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def __init__(self):
    self.n_atoms     = 0                    # number of atoms
    self.box_x       = []                   # xx, xy, xz
    self.box_y       = []                   # yx, yy, yz
    self.box_z       = []                   # zx, zy, zz
    self.list_atoms  = []                   # [atom_idx, x, y, z]
    self.graph       = []                   # [idx_source, idx_Hbond_O1, idx_Hbond_O2]
  
  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
  #                            read_file()                                                #  
  #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
  def read_file(self, infilename):
    n_line = 0
    section_NPGH = False
    section_NX4A = False
    atom_idx     = 0

    with open(infilename) as f:
      for line in f:
        n_line += 1
        
        line_complete = line
        
        line = line.split() # split words of line into a list

        #----------#
        # read box #
        #----------#
        if n_line==1: # box
          if line[0] != "@BOX3" and line[0] != "@BOX9":
            print "[ERROR]. Could not read file (line 1, @BOX3/@BOX9)", infilename
            print line_complete
            sys.exit()

        if n_line==2:
          if len(line) == 3:
            self.box_x = [float(line[0]), 0.0, 0.0]
            self.box_y = [0.0, float(line[1]), 0.0]
            self.box_z = [0.0, 0.0, float(line[2])]
          elif len(line) == 9:
            self.box_x = [float(line[0]), float(line[1]), float(line[2])]
            self.box_y = [float(line[3]), float(line[4]), float(line[5])]
            self.box_z = [float(line[6]), float(line[7]), float(line[8])]
          else:
            print "[ERROR]. Could not read file (line 2, box-dim)", infilename
            print line_complete
            sys.exit()

        #------------#
        # read graph #
        #------------#
        if n_line == 3:
          if line[0] != "@NGPH":
            print "[ERROR]. Could not read file (line 1, @NGPH)", infilename
            print line_complete
            sys.exit()
          section_NPGH=True
          continue

        if n_line > 3 and section_NPGH==True and len(line)==1:
          self.n_atoms = int(line[0])
          continue
        
        if n_line > 3 and section_NPGH == True and len(line)==2:
          line = [int(x) for x in line]
          
          if line[0] == -1 and line[1] == -1: # graph read fully
            section_NPGH=False
            continue
          else:
            self.graph.append([line[0], line[1]])
            continue

        #----------------#
        # read structure #
        #----------------#
        if n_line > 3 and section_NPGH==False and section_NX4A==False:
          if line[0] != "@NX4A":
            print "[ERROR]. Could not read file (@NX4A)", infilename
            print line_complete
            sys.exit()
          section_NX4A=True
          continue
        
        if n_line > 3 and section_NX4A==True and len(line)==1:
          if self.n_atoms != int(line[0]):
            print "[ERROR]. The number of atomts in @NPGH and @NX4A does not match up."
            sys.exit()
        
        elif n_line > 3 and section_NX4A==True and len(line)!=1:
          line = [float(x) for x in line]
          self.list_atoms.append([atom_idx, line[0], line[1], line[2]])
          atom_idx+=1
          continue

    if len(self.list_atoms) != atom_idx:
      print "[ERROR]. The number of atoms does not match up"
      sys.exit()
