"""
class that contains different ways of analysing a trajectory

the trajectory is handed in as an argument (list)
the format of this list is the one used in file_formats.py

traj(all atoms) --> selector(interesting atoms) --> analysis

analysers:
  (+) density profile in z               [OK]
  (+) 2D (xy) distribution               [OK]
  (+) n_particles                        [OK]
  (+) MSD                                [OK]
  (+) RMSD
  (+) distance distribution
  (+) coordination number
  (+) H-bond length distribution         [OK]
  (+) H-bond angle distribution          [OK]
  (+) H-bond lifetime
  (+) inter-layer exchange
  (+) mobility                           [OK]
  (+) dipole distribution
  (+) atom1-atom2 distance distribution  [OK]

input: trajectory (file_format xyz, gro)
output: depends on analyser
"""

import os.path
import math
import sys
import numpy as np

import file_formats as ff
import transform_coordinates as tt
import pbc_computations as pbc

class analyser:
  ##########################################################################################
  #                            INIT                                                        # 
  ##########################################################################################  
  def __init__(self, traj_format="xyz"):
    """ 
    traj should be set manually
      --> after processing trajectory with selector
    """
    if traj_format=="xyz": 
      self.traj = ff.format_xyz() # initialize trajectory
      self.traj_format="xyz"

    if traj_format=="gro": 
      self.traj = ff.format_GRO() # initialize trajectory
      self.traj_format="gro"
    
  ##########################################################################################
  #                            check if traj is there                                      # 
  ##########################################################################################  
  def check_traj(self):
    """
    check if the data in traj_full makes sense

    exits code if problem
    """

    if self.traj.n_atoms == 0:
      print "[ERROR]. n_atoms == 0"
      sys.exit()
    if self.traj.n_frames == 0:
      print "[ERROR]. n_frames == 0"
      sys.exit()
    if len(self.traj.lattice_vector_1) != 3 or \
       len(self.traj.lattice_vector_2) != 3 or \
       len(self.traj.lattice_vector_3) != 3 :
      print "[ERROR]. lattice_vectors not set"
      sys.exit() 

  ##########################################################################################
  #                            PROFILE IN Z                                                # 
  ##########################################################################################  
  def z_profile(self, reference_z=-1.111, list_reference_z=[], outfilename="height_profile.dat"):
    """
    calculate the density in z direction

    reference_z      ... reference_height (not changing during sim)
                         if != -1.111 --> assumed that it was not set
                         a test on (001) MD traj (20-60ps) revealed that the height changes about 0.4 \AA
                         5th percentile: 18.35, 95th percentile: 18.59 ==> 90 percent of data in 0.3 \AA range
                         this is too much for resolution we are aiming for (histo-binning: 0.1 \AA)
    list_reference_z ... indices of atoms for reference_height
    outfilename      ... filename for output
    
    every atom that is not part of list_reference_z will be used for analysis
    """

    self.check_traj()

    outfile = open(outfilename, "a")

    for idx_frame, frame in enumerate(self.traj.list_atoms):
      # monitor progress
      if int(self.traj.n_frames*0.25)==idx_frame:
        print "         25 % done"
      if int(self.traj.n_frames*0.50)==idx_frame:
        print "         50 % done"
      if int(self.traj.n_frames*0.75)==idx_frame:
        print "         75 % done   "
      
      # get/calculate reference height for this frame
      ref_z = 0.0
      if reference_z != -1.111 and len(list_reference_z)==0:
        ref_z = reference_z
      elif reference_z == -1.111 and len(list_reference_z)>0:
        for idx_reference in list_reference_z:
          ref_z += frame[idx_reference][3]/len(list_reference_z)
        # file_refz=open("ref-z.dat", "a")
        # file_refz.write(str(ref_z)+"\n")
      else:
        print "[ERROR]. Z-PROFILE: Wrong format for ref_z..."
        sys.exit()
      
      for idx,atom in enumerate(frame):
        if idx not in list_reference_z:
          h = atom[3] - ref_z 
          outfile.write(str(h) + "\n")

    outfile.close()

  ##########################################################################################
  #                            XY - 2D probability                                         # 
  ##########################################################################################  
  def probability_2D(self, outfilename="", list_name_select=["O"]):
    """
    Gives the lateral positions of selected atoms

    list_name_select ... contains all atom names that are considered for analysis

    """
    self.check_traj()

    if outfilename == "":
      print "[ERROR]. Please specify the output-filename."
      sys.exit()

    outfile = open(outfilename, "a")

    for idx_frame, frame in enumerate(self.traj.list_atoms):
      # monitor progress
      if int(self.traj.n_frames*0.25)==idx_frame:
        print "         25 % done"
      if int(self.traj.n_frames*0.50)==idx_frame:
        print "         50 % done"
      if int(self.traj.n_frames*0.75)==idx_frame:
        print "         75 % done   "
      

      for atom in frame:
        if atom[0] in list_name_select:
          #                 x-coord              y-coord
          outfile.write(str(atom[1]) + " " + str(atom[2]) + "\n")
      
    outfile.close()

  ##########################################################################################
  #                            n particles                                                 # 
  ##########################################################################################  
  def n_particles(self, outfilename="", list_name_select=["O"]):
    """
    calculate the number of particles
    """
    self.check_traj()
    
    if outfilename == "":
      print "[ERROR]. Please specify the output-filename."
      sys.exit()

    outfile = open(outfilename, "a")

    for idx_frame, frame in enumerate(self.traj.list_atoms):
      # monitor progress
      if int(self.traj.n_frames*0.25)==idx_frame:
        print "         25 % done"
      if int(self.traj.n_frames*0.50)==idx_frame:
        print "         50 % done"
      if int(self.traj.n_frames*0.75)==idx_frame:
        print "         75 % done   "
      
      n_particles=0
      
      for atom in frame:
        if atom[0] in list_name_select:
          n_particles+=1
        
      outfile.write(str(n_particles) + "\n")
      
    outfile.close()


  ##########################################################################################
  #                            create dictionary idx_full -> idx_local                     # 
  ##########################################################################################  
  def create_idx_dict(self):
    """
    the indices specified for bonding information are the indices from full trajectory

    here we generate a dictionary, that links these indices to the actual index (0,1,2,...)
    of an atom in a frame

    the atom indices from frame to frame can change, we therefore need a new dictionary
    for every frame
     --> create list of dictionaries
    """

    # makros
    INDEX      = 4

    self.list_dict_idx = []
    for idx_frame, frame in enumerate(self.traj.list_atoms):
      dict_frame={}
      for idx_atom, atom in enumerate(frame):
        dict_frame[atom[INDEX]] = idx_atom

      self.list_dict_idx.append(dict_frame)


  ##########################################################################################
  #                            d-alpha plot                                                # 
  ##########################################################################################  
  def Hbond_d_alpha(self, outfilename="", list_name_select=["O"]):
    """
    calculate the distribution of bondlenghts and angles
    of selected hydrogen bonds
    """

    self.check_traj()

    # makros
    NAME           = 0
    X              = 1
    Y              = 2
    Z              = 3 
    INDEX          = 4
    LIST_ACCEPTORS = 5
    LIST_DONORS    = 6

    # list_donors / list_acceptors
    PARTNER_TYPE   = 0
    O_PARTNER      = 1
    H_PARTNER      = 2
    
    # check if index-dictionary was created already
    # if not: create it!
    try:
      self.list_dict_idx
    except:  
      self.create_idx_dict()

    if outfilename == "":
      print "[ERROR]. No outputfilename specified!"
      sys.exit()

    outfile = open(outfilename, "a")

    for idx_frame, frame in enumerate(self.traj.list_atoms):
      # monitor progress
      if int(self.traj.n_frames*0.25)==idx_frame:
        print "         25 % done"
      if int(self.traj.n_frames*0.50)==idx_frame:
        print "         50 % done"
      if int(self.traj.n_frames*0.75)==idx_frame:
        print "         75 % done   "
      
      for atom in frame:
        if atom[NAME] in list_name_select:
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
          #                            loop over acceptors                                        #  
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
          for partner in atom[LIST_ACCEPTORS]:
            #=======================================================================================#
            #                            bondlength                                                 #  
            #=======================================================================================# 
            p1 = [ atom[X], atom[Y], atom[Z] ]

            p2 = [ frame[self.list_dict_idx[idx_frame][int(partner[O_PARTNER])]][X] , \
                   frame[self.list_dict_idx[idx_frame][int(partner[O_PARTNER])]][Y] , \
                   frame[self.list_dict_idx[idx_frame][int(partner[O_PARTNER])]][Z] ]
            

            d = pbc.pbc_distance ( p1 = p1,                                   \
                                   p2 = p2,                                   \
                                   lattice_v_1 = self.traj.lattice_vector_1 , \
                                   lattice_v_2 = self.traj.lattice_vector_2 , \
                                   lattice_v_3 = self.traj.lattice_vector_3 )

            #=======================================================================================#
            #                            angle                                                      #  
            #=======================================================================================# 
            # p0 is the atom in the middle (p0 - p1 - p2)
            p0 = [ frame[self.list_dict_idx[idx_frame][int(partner[H_PARTNER])]][X] , \
                   frame[self.list_dict_idx[idx_frame][int(partner[H_PARTNER])]][Y] , \
                   frame[self.list_dict_idx[idx_frame][int(partner[H_PARTNER])]][Z] ]

            angle = pbc.pbc_angle ( p0 = p0,                                   \
                                    p1 = p1,                                   \
                                    p2 = p2,                                   \
                                    lattice_v_1 = self.traj.lattice_vector_1 , \
                                    lattice_v_2 = self.traj.lattice_vector_2 , \
                                    lattice_v_3 = self.traj.lattice_vector_3 )
            
            outfile.write(str(d)+"   "+str(angle)+"\n")

          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
          #                            loop over donors                                           #  
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
          for partner in atom[LIST_DONORS]:
            #=======================================================================================#
            #                            bondlength                                                 #  
            #=======================================================================================# 
            p1 = [ atom[X], atom[Y], atom[Z] ]

            p2 = [ frame[self.list_dict_idx[idx_frame][int(partner[O_PARTNER])]][X] , \
                   frame[self.list_dict_idx[idx_frame][int(partner[O_PARTNER])]][Y] , \
                   frame[self.list_dict_idx[idx_frame][int(partner[O_PARTNER])]][Z] ]
            

            d = pbc.pbc_distance ( p1 = p1,                                   \
                                   p2 = p2,                                   \
                                   lattice_v_1 = self.traj.lattice_vector_1 , \
                                   lattice_v_2 = self.traj.lattice_vector_2 , \
                                   lattice_v_3 = self.traj.lattice_vector_3 )

            #=======================================================================================#
            #                            angle                                                      #  
            #=======================================================================================# 
            # p0 is the atom in the middle (p0 - p1 - p2)
            p0 = [ frame[self.list_dict_idx[idx_frame][int(partner[H_PARTNER])]][X] , \
                   frame[self.list_dict_idx[idx_frame][int(partner[H_PARTNER])]][Y] , \
                   frame[self.list_dict_idx[idx_frame][int(partner[H_PARTNER])]][Z] ]

            angle = pbc.pbc_angle ( p0 = p0,                                   \
                                    p1 = p1,                                   \
                                    p2 = p2,                                   \
                                    lattice_v_1 = self.traj.lattice_vector_1 , \
                                    lattice_v_2 = self.traj.lattice_vector_2 , \
                                    lattice_v_3 = self.traj.lattice_vector_3 )
            
            outfile.write(str(d)+"   "+str(angle)+"\n")

    outfile.close()


  ##########################################################################################
  #                            MOBILITIES                                                  # 
  ##########################################################################################  
  def mobility( self,                  \
                outfilename="",        \
                delta_frame=20,        \
                list_name_select=["O"] ) :
    """
    calculate the mobility of selected particles

    mobility(t) = \frac{1}{N} \sum_{i=1}^{N} [r_i(t+\Delta t) - r_i(t)]^2

    N           ... number of particles (can change from frame to frame)
    delta_frame ... to how far in the future do we compare current position 

    mobility is calculated for each frame

    if particle with same index does not exist at (t) and (t+\Delta t)
      --> don't calculate mobility for that frame and atom

    this analyser requires extra information (index of atom)
    """
    self.check_traj()

    # makros
    NAME           = 0
    X              = 1
    Y              = 2
    Z              = 3 
    INDEX          = 4

    if outfilename == "":
      print "[ERROR]. No outputfilename specified!"
      sys.exit()

    outfile = open(outfilename, "a")

    for idx_frame, frame in enumerate(self.traj.list_atoms):
      # monitor progress
      if int(self.traj.n_frames*0.25)==idx_frame:
        print "         25 % done"
      if int(self.traj.n_frames*0.50)==idx_frame:
        print "         50 % done"
      if int(self.traj.n_frames*0.75)==idx_frame:
        print "         75 % done   "
     
      d_sum_squared = 0
      n_particles   = 0

      for atom in frame:
        if atom[NAME] in list_name_select:
          # store coords r_i(t)
          p1 = [ atom[X], atom[Y], atom[Z] ]
          p2 = []

          # look for same atom in frame+delta_frame
          # if end of traj so that idx_frame+delta_frame does not exist --> continue
          if idx_frame+delta_frame>(self.traj.n_frames-1):
            continue

          for atom_2 in self.traj.list_atoms[idx_frame+delta_frame]:
            if atom_2[NAME] in list_name_select:
              if atom_2[INDEX] == atom[INDEX]:
                p2 = [ atom_2[X], atom_2[Y], atom_2[Z] ]
                break

          # if particle p1 is not there in frame+delta_frame: move on
          if len(p2)==0:
            continue

          d = pbc.pbc_distance ( p1 = p1,                                   \
                                 p2 = p2,                                   \
                                 lattice_v_1 = self.traj.lattice_vector_1 , \
                                 lattice_v_2 = self.traj.lattice_vector_2 , \
                                 lattice_v_3 = self.traj.lattice_vector_3 )

          # collect data
          d_sum_squared += d**2
          n_particles += 1

      # --- ATOM LOOP DONE ---
      if n_particles>0:
        outfile.write(str(d_sum_squared/n_particles) + "\n")
  
    # --- FRAME LOOP DONE ---
    outfile.close()


  ##########################################################################################
  #                            MSD                                                         # 
  ##########################################################################################  
  def MSD ( self,            \
            outfilename="",  \
            skip_frames=50,  \
            timestep=1.0     ) :
    """
    calculates the MSD

    here the number of particles has to be constant 
      --> we can only apply this to all water molecules in trajectory
    and not for example a subset (such as water in first L

    for more info check:
    http://utkstair.org/clausius/docs/che548/pdf/selfD.pdf

    some technical details:
      Np        ... number of particles
      Nt        ... number of time frames
      No        ... = int(Nt/2): number of considererd origins
                    this also specifies for how many future timesteps will be considered
                    
                    this way every delta_t has the same number of datapoints
                    which means that all of them are equally weighted

      e.g. Nt: 1 2 3 4 5 6  --> Nt/2=3
      calculate MSD 1-2, 1-3, 1-4, 2-3, 2-4, 2-5, 3-4, 3-5, 3-6
      each delta_t has Nt/2*Np datapoints --> this is generally used to get statistics

    this analysis needs to be performed on one concenated trajectory

    it is assumed that the order of atoms does not change from frame to frame
      to make sure we don't do anything stupid, this is checked on the run

    skip_frames ... skipping this number of frames each time
    timestep    ... timestep between different frames (in fs)
    """

    self.check_traj()
         
    X     = 1
    Y     = 2
    Z     = 3
    INDEX = 4

    if outfilename == "":
      print "[ERROR]. Please specify the output-filename."
      sys.exit()

    outfile = open(outfilename, "a")

    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            split traj                                                 #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    # split trajectory into sub-trajectory containing each skip_frames-th frame
    traj_split = ff.format_xyz() 

    traj_split.list_atoms=self.traj.list_atoms[::skip_frames]
    traj_split.n_frames=len(traj_split.list_atoms)

    print "\t(+) split trajectory from", self.traj.n_frames, "frames to", len(traj_split.list_atoms), "frames"

    # calculate n_origins
    n_origins = int(traj_split.n_frames/2)

    # initialize list that saves MSD for each delta_t
    # we will have n_origins datapoints for each origin
    # the simple idea:
    #    list_msd = [0.0] * n_origins
    # does not work!
    # because each element in list changes if you change one (weird?)
    # workaround:
    list_msd = []
    for i in range(n_origins):
      list_msd.append([])
    
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            loop over all time origins                                 #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    for idx_origin in range(n_origins):
      # monitor progress
      if int(n_origins*0.10)==idx_origin:
        print "         10 % done"
      if int(n_origins*0.20)==idx_origin:
        print "         20 % done"
      if int(n_origins*0.30)==idx_origin:
        print "         30 % done   "
      if int(n_origins*0.40)==idx_origin:
        print "         40 % done"
      if int(n_origins*0.50)==idx_origin:
        print "         50 % done"
      if int(n_origins*0.60)==idx_origin:
        print "         60 % done   "
      if int(n_origins*0.70)==idx_origin:
        print "         70 % done   "
      if int(n_origins*0.80)==idx_origin:
        print "         80 % done   "
      if int(n_origins*0.90)==idx_origin:
        print "         90 % done   "

      list_reference = []

      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            loop over all relevant timeframes                          #  
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      for idx_frame, frame in enumerate(traj_split.list_atoms[idx_origin:idx_origin+n_origins]):
        
        msd = 0

        for idx_atom,atom in enumerate(frame):
          if len(atom) < 5:
            print "[ERROR]. In MSD procedure. An atom has not the right format"
            print "         Frame: ", idx_frame, "atom: ", idx_atom
            print atom
            sys.exit()

          # store x0, y0, z0 for first frame 
          if idx_frame == 0:
            list_reference.append(atom)
          
          # calculate MSD for this atom
          else:
            # check if trajectory oredered correctly
            if list_reference[idx_atom][INDEX] != atom[INDEX]:
              print "[ERROR]. In MSD procedure: the trajectory is ordered incorrectly"
              sys.exit()

            p0 = [ list_reference[idx_atom][X], list_reference[idx_atom][Y], list_reference[idx_atom][Z] ]
            p1 = [ atom[X], atom[Y], atom[Z] ]

            d = pbc.pbc_distance( p0,                                        \
                                  p1,                                        \
                                  lattice_v_1 = self.traj.lattice_vector_1 , \
                                  lattice_v_2 = self.traj.lattice_vector_2 , \
                                  lattice_v_3 = self.traj.lattice_vector_3 )
            msd += d**2
         
        # normalize MSD to number of atoms
        if int(self.traj.n_atoms) != len(frame):
          print "[ERROR]. In MSD procedure: the number of atoms does not make sense"
          print self.traj.n_atoms, len(frame)
          sys.exit()

        msd = msd/int(self.traj.n_atoms)

        # store MSD for correct delta_t value
        #  --> these values are stored in list_msd
        list_msd[idx_frame].append(msd)
            
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            MSD analysis done --> write output                         #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    # check if list_msd makes sense:
    if len(list_msd) != n_origins:
      print "[ERROR]. In MSD procedure: list_msd has wrong length"
      sys.exit()

    # write column names
    outfile.write("time ")
    for i,datapoint in enumerate(list_msd[0]):
      outfile.write("obs_"+str(i)+ " ")
    outfile.write("\n")

    # write data
    for idx, datapoint in enumerate(list_msd):
      outfile.write(str(idx*timestep*skip_frames)+" ")
      for i in datapoint:
        outfile.write(str(i)+ " ")
      outfile.write("\n")

    outfile.close()

  ##########################################################################################
  #                            LIFETIMES                                                   # 
  ##########################################################################################  
  def lifetimes ( self, \
                  outfilename="", \
                  min_t=10,       \
                  list_name_select=["O"]):
    """
    calculate the lifetime of a particle in a selection

    the lifetime is based wether or not the particle appears/disappears from the selection
     --> to do that we look at indizes of particles in selection

    min_t ... minimum time (specified in frames) the particle has to be gone/appeared in selection
              to be counted
              this is necessary in order to make sure short-term fluctuations don't influence the results

              [NEW ATOM APPEARS]: new atom index pops at at frame 10, min_t=3
              then this atom index has to be there in frame 11, frame 12 and frame 13 to be considered new
              i.e. in total new atom has to be there for at least min_t+1 frames

              [ATOM LEAVES SELECTION]: atom with new index dissapears in frame 10, min_t=3
              then this atom index cannot re-appear in frame 11, 12 or 13 to be considered gone
              i.e. in total the old atom has to be gone for at least min_t+1 frames

    list_name_select ... atom names that are considered for analysis
    """
    self.check_traj()

    # makros
    NAME           = 0
    X              = 1
    Y              = 2
    Z              = 3 
    INDEX          = 4

    if outfilename == "":
      print "[ERROR]. No outputfilename specified!"
      sys.exit()


    list_indices = [] # stores the indices that ever were in selection, 
    #                   + data when this idx entered selection
    #                   + data when this idx left selection

    list_tracked_idx = [] # stores the indices that are currently tracked

    for idx_frame, frame in enumerate(self.traj.list_atoms):
      # monitor progress
      if int(self.traj.n_frames*0.25)==idx_frame:
        print "         25 % done"
      if int(self.traj.n_frames*0.50)==idx_frame:
        print "         50 % done"
      if int(self.traj.n_frames*0.75)==idx_frame:
        print "         75 % done   "
     
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            CHECK IF NEW ATOMS IN SELECTION                            #  
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      for atom in frame:
        if atom[NAME] in list_name_select:

          # first frame: store indizes occuring in this frame
          if idx_frame == 0:
            list_indices.append([atom[INDEX], idx_frame])
            list_tracked_idx.append(atom[INDEX])

          # not first frame
          else:

            if atom[INDEX] in list_tracked_idx:
              continue # if atom was there before already, move on
            else:
              # we found a new atom that wasn't there in previous frame
              #   is atom[INDEX] is gonna remain in selection within the next min_t periods?
              #     if yes --> new lifetime of atom[INDEX] starts
              #     if no  --> move on
              
              flag_new_atom_present=True # if atom is gone from selection only once
              #                            this flag will be set to FALSE


              # check if there are still enough future-frames available:
              # if there are 5 frames [0,1,2,3,4] , and min_t=2, then we need to abort at frame 3 and 4
              if len(self.traj.list_atoms)<=(idx_frame+min_t):
                break

              for future_frame in self.traj.list_atoms[idx_frame+1:idx_frame+min_t+1]:
                flag_new_atom_found=False
                for future_atom in future_frame:
                  if future_atom[INDEX]==atom[INDEX]:
                    flag_new_atom_found=True

                if flag_new_atom_found==False:
                  flag_new_atom_present=False
                  break
             
              # atom DOES NOT disappears within t_min time periods
              #  therefore is counted as birth of a new particle
              if flag_new_atom_present==True:
                list_indices.append([atom[INDEX], idx_frame])  # add new index 
                list_tracked_idx.append(atom[INDEX])           # track this atom too now
        # --- END ATOM LOOP ---

      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            CHECK IF ATOMS LEFT SELECTION                              #  
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      
      # if this is first frame --> move on
      if idx_frame == 0:
        continue

      # make a list with current indizes in selection
      list_current_indices = []
      for atom_current in frame:
        list_current_indices.append(atom_current[INDEX])
      
      # loop over all tracked indices and see if they are still here
      for prev_idx in list_tracked_idx:
        if prev_idx not in list_current_indices:
          
          flag_old_atom_returned=False # if atom returns only once to selection
          #                              this flag will be set to True, meaning that atom did return
        
          # check if there are still enough future-frames available:
          # if there are 5 frames [0,1,2,3,4] , and min_t=2, then we need to abort at frame 3 and 4
          if len(self.traj.list_atoms)<=(idx_frame+min_t):
            break

          for future_frame in self.traj.list_atoms[idx_frame+1:idx_frame+min_t+1]:
            for future_atom in future_frame:
              if future_atom[INDEX]==prev_idx:
                flag_old_atom_returned=True
                break

            if flag_old_atom_returned==True:
              break

          # atom DID NOT return within min_t periods --> lifetime over
          if flag_old_atom_returned == False:
            # add end of lifetime to atom that we lost
            # to do that: find the element of list_indices that has the correct index
            #             AND has only two sub-elements (idx and lifetime_start)

            # stop tracking this atom from now on
            list_tracked_idx.remove(prev_idx)

            # store information of death in list_indices
            flag_found_idx=False
            for i in list_indices:
              if i[0]==prev_idx and len(i)==2:
                i.append(idx_frame-1) # this is the end of ith's lifetime
                #                       -1 because particle "died" in current frame
                flag_found_idx=True
              
            if flag_found_idx==False:
              print "[ERROR]. Could not find atom", prev_idx, "in list_indices"
              print list_indices
              sys.exit()
      # --- END FRAME LOOP ---
    
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            CALCULATE LIFETIMES                                        #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    # for all particles still alive --> there is a value missing (no end of lifetime)
    # --> for these we do not really know how long they would have lived, so we ignore them
    #     and only calculate lifetimes for particles that have completed a whole lifecycle
    list_lifetime = []
    
    for i in list_indices:
      if len(i)==3:
        list_lifetime.append(i[2] - i[1])
    
    # sort list according to index
    list_indices.sort()

    # write outputfile with completed lifetimes
    outfile = open(outfilename, "a")
    outfile.write("     index      birth      death   lifetime\n")

    for i in list_indices:
      if len(i)==3:
        outfile.write(str(i[0]).rjust(10)+" "+str(i[1]).rjust(10)+" "+str(i[2]).rjust(10)+" "+str(i[2]-i[1]).rjust(10)+"\n")

    outfile.close()

    # write outputfile with uncompleted lifetimes
    outfile = open(outfilename+"_INCOMPLETE", "w")
    outfile.write("     index      birth      death   lifetime\n")

    for i in list_indices:
      if len(i)==2:
        outfile.write(str(i[0]).rjust(10)+" "+str(i[1]).rjust(10)+"\n")

    outfile.close()

  ##########################################################################################
  #                            ATOM1 - ATOM2 distance distribution                         # 
  ##########################################################################################  
  def a1_a2_distance_dist(self, selection_type, selection_1, selection_2, outfilename="distance_dist.dat", cutoff=15.0):
    """
    calculate the distance distribution between selection_1 and selection-2

    selection_type: "atom_name"
    cutoff: ignore distances larger than that

    """
    outfile = open(outfilename, "a")

    if selection_type not in ["atom_name"]:
      print "[ERROR] a1_a2_distance_dist(). Selection type", selection_type, "not recognized!"
      sys.exit()

    for idx_frame, frame in enumerate(self.traj.list_atoms):
      # monitor progress
      if int(len(self.traj.list_atoms)*0.25)==idx_frame:
        print "         25 % done"
      if int(len(self.traj.list_atoms)*0.50)==idx_frame:
        print "         50 % done"
      if int(len(self.traj.list_atoms)*0.75)==idx_frame:
        print "         75 % done   "
      
      for idx1, atom1 in enumerate(frame):
        for idx2, atom2 in enumerate(frame):
          if idx2 >= idx1:
            continue

          if atom1[self.traj.ATOMNAME] == selection_1 and atom2[self.traj.ATOMNAME] == selection_2:
            d = pbc.pbc_distance( [atom1[self.traj.X_COORD], atom1[self.traj.Y_COORD], atom1[self.traj.Z_COORD]], \
                                  [atom2[self.traj.X_COORD], atom2[self.traj.Y_COORD], atom2[self.traj.Z_COORD]], \
                                  self.traj.lattice_vector_1, \
                                  self.traj.lattice_vector_2, \
                                  self.traj.lattice_vector_3  )
        
            
            if d*10.0 < cutoff:
              outfile.write(str(d*10.0) + "\n")

    outfile.close()
