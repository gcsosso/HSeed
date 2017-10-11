"""
class that contains different ways of selecting atoms for analysis

traj(all atoms) --> selector(interesting atoms) --> analysis

different selection mechanisms:
  (+) [OK] all
  (+) [OK] index
  (+) geometric criterion ( eg z_low <= z <= z_high)
  (+) H-bonded to a motive (eg SiOH group)

input: trajectory (file_format xyz), len: n_frames
output: selected atoms (file_format xyz), len: n_frames
"""

import os.path
import math
import sys
import numpy as np

import file_formats as ff
import transform_coordinates as tt
import pbc_computations as pbc

class selector:
  ##########################################################################################
  #                            init                                                        # 
  ##########################################################################################  
  def __init__(self):
    self.traj_full = ff.format_xyz() # initialize full trajectory
    self.traj_sel  = ff.format_xyz() # initialize trajectory with selected atoms

  ##########################################################################################
  #                            check if traj is there                                      # 
  ##########################################################################################  
  def check_traj(self):
    """
    check if the data in traj_full makes sense

    exits code if problem
    """

    if self.traj_full.n_atoms == 0:
      print "[ERROR]. n_atoms == 0"
      sys.exit()
    if self.traj_full.n_frames == 0:
      print "[ERROR]. n_frames == 0"
      sys.exit()
    if len(self.traj_full.lattice_vector_1) != 3 or \
       len(self.traj_full.lattice_vector_2) != 3 or \
       len(self.traj_full.lattice_vector_3) != 3 :
      print "[ERROR]. lattice_vectors not set"
      sys.exit() 

  ##########################################################################################
  #                            ALL                                                         # 
  ##########################################################################################  
  def select_all(self):
    self.check_traj()   # check if traj makes sense

    self.traj_sel = self.traj_full

  ##########################################################################################
  #                            INDEX                                                       # 
  ##########################################################################################  
  def select_index(self, index_low=-1, index_high=-1, list_index_select=[], list_name_select=["O", "H"]):
    """
    select atoms which have the correct index 

    the indexing starts with 0 and ends at n_atoms-1! (like vmd)

    option1: specify index_low and index_high, all of the ones in between are selected
             index_low <= select <= index_high
    
    option2: specify all indices explicitly (if index_low = index_high = -1)
             list_index_select contains all the indices explicitly

    option3: specify index_range AND some other indices explicitly

    list_name_select contains all the atom names that are used for analysis
    this is applied only to option1...
    atoms that are indexed explicitly are always included
    """
    self.check_traj()   # check if traj makes sense
  
    # define makros
    INDEX_RANGE    = 0 
    INDEX_EXPLICIT = 1
    INDEX_BOTH     = 2

    # sanity check
    if index_low == -1 and index_high == -1 and list_index_select == []:
      print "[ERROR]. Selection criterion INDEX not set"
      sys.exit()
    elif index_low >=0 and index_high > index_low and len(list_index_select)==0:
      mode = INDEX_RANGE
    elif index_low ==-1 and index_high ==-1 and len(list_index_select)>0:
      mode = INDEX_EXPLICIT
    elif index_low >=0 and index_high > index_low and len(list_index_select)>0:
      mode = INDEX_BOTH
    else:
      print "[ERROR]. Selection criterion INDEX not set"
      sys.exit()
    
    #self.traj_sel.list_atoms = []
    
    # loop over frames
    for idx_frame, frame in enumerate(self.traj_full.list_atoms):
      if int(self.traj_full.n_frames*0.25)==idx_frame:
        print "         25 % done"
      if int(self.traj_full.n_frames*0.50)==idx_frame:
        print "         50 % done"
      if int(self.traj_full.n_frames*0.75)==idx_frame:
        print "         75 % done   "

      list_atom_select = []
      # loop over atoms in frame
      for idx,atom in enumerate(frame):
        if mode == INDEX_RANGE:
          if idx >= index_low and idx <= index_high: # idx in correct range?
            if atom[0] in list_name_select: # should we care about this atom name?
              list_atom_select.append(atom + [[idx]])
        elif mode == INDEX_EXPLICIT:
          if idx in list_index_select:
            list_atom_select.append(atom)
        elif mode == INDEX_BOTH:
          if idx >= index_low and idx <= index_high: # idx in correct range?
            if atom[0] in list_name_select: # should we care about this atom name?
              list_atom_select.append(atom + [[idx]])
          if idx in list_index_select:
            # only add it if it is not included in INDEX_RANGE already
            if idx < index_low or idx > index_high: 
              list_atom_select.append(atom + [[idx]])

      # add all atoms from this frame as a new element in traj_sel
      self.traj_sel.list_atoms.append(list_atom_select)

    # set title
    self.traj_sel.title = ""
    
    # set n_atoms and n_frames
    # n_atoms ... largest number of atoms in a frame
    self.traj_sel.n_atoms  = len(max(self.traj_sel.list_atoms,key=len))
    self.traj_sel.n_frames = len(self.traj_sel.list_atoms)

    if self.traj_sel.n_frames != self.traj_full.n_frames:
      print "[ERROR]. INDEX SELECTOR: Something went wrong, n_frames changed!"
      sys.exit()

    # set lattice_vectors
    self.traj_sel.lattice_vector_1 = self.traj_full.lattice_vector_1
    self.traj_sel.lattice_vector_2 = self.traj_full.lattice_vector_2
    self.traj_sel.lattice_vector_3 = self.traj_full.lattice_vector_3

  ##########################################################################################
  #                            HEIGHT (z)                                                  # 
  ##########################################################################################  
  def select_height( self,                   \
                     index_min,              \
                     index_max,              \
                     h_min,                  \
                     h_max,                  \
                     list_index_ref_z=[],    \
                     reference_z=-1.111,     \
                     list_name_select=["O"], \
                     flag_project_inside=True):
    """
    select atoms which have the correct height (in interval h_min <= h < h_max)

    index_min           ... minimum index to be considered
    index_max           ... maximum index to be considered
    
    option1: specify list_index_ref_z
              the reference_z will be calculated as average over all z of the elements
    
    option2: specify ref_z
              this will be used as ref_z

    list_name_select    ... contains all the atom names that are used for analysis
    flag_project_inside ... decide if atom coords shall be projected inside simulation box
    """
    self.check_traj()   # check if traj makes sense
    
    DYNAMIC_HEIGHT = 0
    STATIC_HEIGHT  = 1

    # figure out mode to calculate reference height
    if len(list_index_ref_z)>0 and reference_z==-1.111:
      mode = DYNAMIC_HEIGHT
    elif len(list_index_ref_z)==0 and reference_z!=-1.111:
      mode = STATIC_HEIGHT
    else:
      print "[ERROR]. Reference height_z not set"
      sys.exit()
    
    # loop over frames
    for idx_frame, frame in enumerate(self.traj_full.list_atoms):
      if int(self.traj_full.n_frames*0.25)==idx_frame:
        print "         25 % done"
      if int(self.traj_full.n_frames*0.50)==idx_frame:
        print "         50 % done"
      if int(self.traj_full.n_frames*0.75)==idx_frame:
        print "         75 % done   "

      # calculate reference height for this frame
      ref_z = 0.0
      if mode == DYNAMIC_HEIGHT:
        for idx_reference in list_index_ref_z:
          ref_z += frame[idx_reference][3]/len(list_index_ref_z)
      elif mode == STATIC_HEIGHT:
        ref_z = reference_z 
      
      list_atom_select = []
      
      # loop over atoms in frame
      for idx,atom in enumerate(frame):
        
        # check if index is OK
        if idx>= index_min and idx<= index_max:
          
          # check if atom type not relevant for analysis:
          if atom[0] not in list_name_select:
            continue

          # check if height is OK
          if atom[3] >= (ref_z+h_min) and atom[3] < (ref_z+h_max):
            
            # project atom into inside of box if necessary
            if flag_project_inside == True:
              #  to do so, we convert x and y into fractional coordinates
              #  if one of them happens to be out of box --> PBC
              coords_frac = tt.cart2frac( self.traj_full.lattice_vector_1 , \
                                          self.traj_full.lattice_vector_2 , \
                                          self.traj_full.lattice_vector_3 , \
                                          [atom[1], atom[2], atom[3]] )

              # pbc
              # the int() makes sure that values end up in interval [0,1]
              # usually: pbc with round!
              coords_frac[0] = coords_frac[0] - int(coords_frac[0])
              coords_frac[1] = coords_frac[1] - int(coords_frac[1])

              if coords_frac[0] < 0:
                coords_frac[0]+=1
              if coords_frac[1] < 0:
                coords_frac[1]+=1

              # convert back to cartesian coords
              coords_cart = tt.frac2cart( self.traj_full.lattice_vector_1 , \
                                          self.traj_full.lattice_vector_2 , \
                                          self.traj_full.lattice_vector_3 , \
                                          coords_frac )
              
              atom[1] = coords_cart[0]
              atom[2] = coords_cart[1]
              atom[3] = coords_cart[2]

              if atom[1]<0 or atom[2]<0 or atom[3]<0:
                print "coords < 0.0"
                print idx_frame, atom
                print coords_frac
                sys.exit()

            list_atom_select.append(atom+[[idx]])

      # --- ATOM LOOP DONE ---
      # add all atoms from this frame as a new element in traj_sel
      self.traj_sel.list_atoms.append(list_atom_select)

    # --- FRAME LOOP DONE ---
    # set title
    self.traj_sel.title = ""
    
    # set n_atoms and n_frames
    # n_atoms ... largest number of atoms in a frame
    self.traj_sel.n_atoms  = len(max(self.traj_sel.list_atoms,key=len))
    self.traj_sel.n_frames = len(self.traj_sel.list_atoms)

    if self.traj_sel.n_frames != self.traj_full.n_frames:
      print "[ERROR]. HEIGHT_Z SELECTOR: Something went wrong, n_frames changed!"
      sys.exit()

    # set lattice_vectors
    self.traj_sel.lattice_vector_1 = self.traj_full.lattice_vector_1
    self.traj_sel.lattice_vector_2 = self.traj_full.lattice_vector_2
    self.traj_sel.lattice_vector_3 = self.traj_full.lattice_vector_3

  ##########################################################################################
  #                            COORDINATING                                                # 
  ##########################################################################################  
  def select_coordination ( self,                        \
                            index_min,                   \
                            index_max,                   \
                            tol_distance,                \
                            list_index_partners=[],      \
                            list_name_select=["O"],      \
                            flag_whole_water=False,      \
                            flag_include_partners=False, \
                            flag_project_inside=True):
    """
    select all atoms that form a coordinative bond with specified motives

    index_min ... minimum index of relevant analyt
    index_max ... maximum index of relevant analyt
    tol_distance ... maximum distance to be considered coordinative bond
    
    list_index_partners ... list with all bonding partners that are considered
    list_name_select    ... atom types which are selected
   
    flag_whole_water      ... determines if whole water molecules is selected or not
    flag_include_partners ... include coordination partners in selection?
    flag_project_inside   ... project atoms inside simulation box?

    """
    self.check_traj()   # check if traj makes sense

    # figure out if arguments OK
    if len(list_index_partners)==0:
      print "[ERROR]. list_index_partners not set"
      sys.exit()
   
    # makros
    X = 1
    Y = 2
    Z = 3

    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            LOOP OVER FRAMES                                           #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    for idx_frame, frame in enumerate(self.traj_full.list_atoms):
      if int(self.traj_full.n_frames*0.25)==idx_frame:
        print "         25 % done"
      if int(self.traj_full.n_frames*0.50)==idx_frame:
        print "         50 % done"
      if int(self.traj_full.n_frames*0.75)==idx_frame:
        print "         75 % done   "
   
      list_atom_select = []

      # include bonding_partners if necessary
      if flag_include_partners == True:
        for partner in list_index_partners:
            tmp_atom=frame[partner]
            tmp_atom[0]+="_bp"
            tmp_atom.append([partner])
            list_atom_select.append(tmp_atom)

      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            LOOP OVER ATOMS IN FRAME                                   #  
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      for idx,atom in enumerate(frame):
        
        # check if index is OK
        if idx>= index_min and idx<= index_max:
          
          # check if atom type not relevant for analysis:
          if atom[0] not in list_name_select:
            continue
          
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
          #                            FIND Ow IN OK DISTANCE FOR COORDINATION                    #  
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
          for partner in list_index_partners:
            
            p1 = [ atom[X], atom[Y], atom[Z] ]
            p2 = [ frame[partner][X] , \
                   frame[partner][Y] , \
                   frame[partner][Z] ]
            
            d = pbc.pbc_distance ( p1 = p1,                                        \
                                   p2 = p2,                                        \
                                   lattice_v_1 = self.traj_full.lattice_vector_1 , \
                                   lattice_v_2 = self.traj_full.lattice_vector_2 , \
                                   lattice_v_3 = self.traj_full.lattice_vector_3 )
            if d<=tol_distance:
              list_atom_select.append(atom+[[idx]])

      # --- ATOM LOOP DONE ---
      # add all atoms from this frame as a new element in traj_sel
      self.traj_sel.list_atoms.append(list_atom_select)

    # --- FRAME LOOP DONE ---
    # set title
    self.traj_sel.title = ""
    
    # set n_atoms and n_frames
    # n_atoms ... largest number of atoms in a frame
    self.traj_sel.n_atoms  = len(max(self.traj_sel.list_atoms,key=len))
    self.traj_sel.n_frames = len(self.traj_sel.list_atoms)

    if self.traj_sel.n_frames != self.traj_full.n_frames:
      print "[ERROR]. COORDINATION SELECTOR: Something went wrong, n_frames changed!"
      sys.exit()

    # set lattice_vectors
    self.traj_sel.lattice_vector_1 = self.traj_full.lattice_vector_1
    self.traj_sel.lattice_vector_2 = self.traj_full.lattice_vector_2
    self.traj_sel.lattice_vector_3 = self.traj_full.lattice_vector_3


  ##########################################################################################
  #                            H-BONDING                                                   # 
  ##########################################################################################  
  def select_Hbond ( self,                        \
                     index_min,                   \
                     index_max,                   \
                     tol_length,                  \
                     tol_angle,                   \
                     list_index_partners=[],      \
                     list_name_select=["O"],      \
                     flag_whole_water=False,      \
                     flag_include_partners=False, \
                     flag_project_inside=True):
    """
    select all atoms that are form hydrogen bonds with specified motives
    
    index_min ... minimum index of relevant analyt
    index_max ... maximum index of relevant analyt
    tol_length ... maximum acceptable Hbond (O-O) length
    tol_angle ... minimum acceptable Hbond (O-H-O) angle
    
    list_index_partners ... list with all bonding partners that are considered
                            the format of this list is the following
                            [[A1], [D1O, D1H1]], [[A2], [D2O, D2H]]
                            in this example: two different potential bonding partners, 1 and 2
                            A ... index for acceptor atom (eg hydroxyl O)
                            D ... D1O donor O, D1H donor H 
                                  (if necessary: add donor H1 and donor H2 in future)
                            if either [A] or [D] are empty ... only acceptor/donor role considered
    list_name_select    ... atom types which are selected
   
    flag_whole_water      ... determines if whole water molecules is selected or not
    flag_include_partners ... include bonding partners in selection?
    flag_project_inside   ... project atoms inside simulation box?
    """
    self.check_traj()   # check if traj makes sense

    # figure out if arguments OK
    if len(list_index_partners)==0:
      print "[ERROR]. list_index_partners not set"
      sys.exit()
    for partner in list_index_partners:
      if len(partner) != 2 or (len(partner[0]) != 1 and len(partner[1]) != 2):
        print "[ERROR]. The partner list does not make sense."
        print "         Length of individual members of list wrong"
        print "         Read the code..."
        print partner[0]
        print partner[1]
        sys.exit()
      if len(partner[1])>0:
        if partner[0][0] != partner[1][0] :
          print "[ERROR]. The partner list does not make sense."
          print "         Indizes wrong"
          print "         Read the code..."
          print partner[0][0], partner[1][0]
          sys.exit()
    
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            LOOP OVER FRAMES                                           #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    for idx_frame, frame in enumerate(self.traj_full.list_atoms):
      if int(self.traj_full.n_frames*0.25)==idx_frame:
        print "         25 % done"
      if int(self.traj_full.n_frames*0.50)==idx_frame:
        print "         50 % done"
      if int(self.traj_full.n_frames*0.75)==idx_frame:
        print "         75 % done   "
   
      list_atom_select = []

      # include bonding_partners if necessary
      if flag_include_partners == True:
        for partner in list_index_partners:
          # find all atoms that belong to partner
          # also mark them as bonding partner for analyser
          # partner as ONLY  Hbond acceptor
          if len(partner[0]) == 1 and len(partner[1])==0:
            tmp_atom=frame[partner[0][0]]
            tmp_atom[0]+="_bp"
            tmp_atom.append([partner[0][0]])
            list_atom_select.append(tmp_atom)
          # partner as Hbond donor (and potentially acceptor too)
          #  the second element is H
          elif len(partner[1]) == 2:
            # append O
            tmp_atom=frame[partner[1][0]]
            tmp_atom[0]+="_bp"
            tmp_atom.append([partner[1][0]])
            list_atom_select.append(tmp_atom)
            
            # append H
            tmp_atom=frame[partner[1][1]]
            tmp_atom[0]+="_bp"
            tmp_atom.append([partner[1][1]])
            list_atom_select.append(tmp_atom)

          else:
            print "[ERROR]. Hbond-selector. The partner-list does not make sense!"
            sys.exit()

      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            LOOP OVER ATOMS IN FRAME                                   #  
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      for idx,atom in enumerate(frame):
        
        # check if index is OK
        if idx>= index_min and idx<= index_max:
          
          # check if atom type not relevant for analysis:
          if atom[0] not in list_name_select:
            continue
          
          list_donors    = []  # complete list of members in list_index_partners that donate Hbond to atom
          list_acceptors = []  # complete list of members in list_index_partners that accept Hbond from atom
          
          flag_Hbond_found = False # keeps track if Hbond was found or not
          
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
          #                            FIND Ow IN OK DISTANCE FOR HBOND                           #  
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
          for partner in list_index_partners:
            coord_partner = []

            # get xyz of partner
            if len(partner[0])==1:
              coord_partner = [ frame[partner[0][0]][1], frame[partner[0][0]][2], frame[partner[0][0]][3] ] 
            else:
              coord_partner = [ frame[partner[1][0]][1], frame[partner[1][0]][2], frame[partner[1][0]][3] ] 
            
            # calculate cartesian distance vector between partner and atom
            dist_v_cart = [ atom[1] - coord_partner[0] , \
                            atom[2] - coord_partner[1] , \
                            atom[3] - coord_partner[2] ]
            
            # convert distance vector to fractional coords
            dist_v_frac = tt.cart2frac( self.traj_full.lattice_vector_1 , \
                                        self.traj_full.lattice_vector_2 , \
                                        self.traj_full.lattice_vector_3 , \
                                        dist_v_cart                     )
                                        
            # apply pbc
            dist_v_frac[:] = [ d-round(d) for d in dist_v_frac ]
            
            # convert back to cartesian coords  
            dist_v_cart = tt.frac2cart( self.traj_full.lattice_vector_1 , \
                                        self.traj_full.lattice_vector_2 , \
                                        self.traj_full.lattice_vector_3 , \
                                        dist_v_frac                     ) 
                                        
            # calculate distance
            d = math.sqrt(dist_v_cart[0]**2 + dist_v_cart[1]**2 + dist_v_cart[2]**2) 
            
            # check if distance is OK to be considererd for Hbond
            if d > tol_length:
              continue
           
            list_donors_tmp    = []  # members in list_index_partners that donate Hbond to atom (expected len: 0 or 1)
            list_acceptors_tmp = []  # members in list_index_partners that accept Hbond from atom (expected len: 0 or 1)
            
            #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
            #                            check if Hbond is accepted by atom                         #  
            #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
            #   partner must have donating ability
            #   i.e. len(partner[1]>1)
            #   for now only one donating group (1 H) supported
            #   can be changed later --> iterate through list of H
            if len(partner[1])==2:
              # get coordinates for O and H from partner
              coord_pO = [ frame[partner[1][0]][1], frame[partner[1][0]][2], frame[partner[1][0]][3] ]
              coord_pH = [ frame[partner[1][1]][1], frame[partner[1][1]][2], frame[partner[1][1]][3] ]
              
              # calculate vector H->O1 and H->O2
              v_HO_1_cart = [ -coord_pH[0] + coord_pO[0], \
                              -coord_pH[1] + coord_pO[1], \
                              -coord_pH[2] + coord_pO[2]  ]
                         
              v_HO_2_cart = [ -coord_pH[0] + atom[1], \
                              -coord_pH[1] + atom[2], \
                              -coord_pH[2] + atom[3]  ] 
                         
              # convert them to fractional coords           
              v_HO_1_frac = tt.cart2frac( self.traj_full.lattice_vector_1 , \
                                          self.traj_full.lattice_vector_2 , \
                                          self.traj_full.lattice_vector_3 , \
                                          v_HO_1_cart                     )      
              v_HO_2_frac = tt.cart2frac( self.traj_full.lattice_vector_1 , \
                                          self.traj_full.lattice_vector_2 , \
                                          self.traj_full.lattice_vector_3 , \
                                          v_HO_2_cart                     ) 
            
              # apply pbc
              v_HO_1_frac[:] = [ d-round(d) for d in v_HO_1_frac ]
              v_HO_2_frac[:] = [ d-round(d) for d in v_HO_2_frac ]
              
              # convert them to back to cartesian coords           
              v_HO_1_cart = tt.frac2cart( self.traj_full.lattice_vector_1 , \
                                          self.traj_full.lattice_vector_2 , \
                                          self.traj_full.lattice_vector_3 , \
                                          v_HO_1_frac                     )      
              v_HO_2_cart = tt.frac2cart( self.traj_full.lattice_vector_1 , \
                                          self.traj_full.lattice_vector_2 , \
                                          self.traj_full.lattice_vector_3 , \
                                          v_HO_2_frac                     )
              
              # calculate angle
              alpha = np.arccos( np.dot (v_HO_1_cart, v_HO_2_cart) / \
                                 (np.linalg.norm(v_HO_1_cart)*np.linalg.norm(v_HO_2_cart)) )
                                 
              # convert to degree
              alpha = alpha*360.0/(2*np.pi)
              
              # found a Hbond?
              if alpha > tol_angle:
                list_acceptors_tmp.append([partner[1][0], partner[1][1]])
                
            #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
            #                            check if Hbond is donated by atom                          #  
            #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
            if len(partner[0])==1:
              # get coordinates for O from partner
              coord_pO = [ frame[partner[0][0]][1], frame[partner[0][0]][2], frame[partner[0][0]][3] ]  
              
              # get coordinates for water molecule
              # assumption: water is stored as O H H in trajectory
              coord_wO  = [ atom[1], atom[2], atom[3] ]
              coord_wH1 = [ frame[idx+1][1], frame[idx+1][2], frame[idx+1][3] ]
              coord_wH2 = [ frame[idx+2][1], frame[idx+2][2], frame[idx+2][3] ] 
              
              # calculate vector H1->O1 and H1->O2
              v_H1O_1_cart = [ -coord_wH1[0] + coord_pO[0], \
                               -coord_wH1[1] + coord_pO[1], \
                               -coord_wH1[2] + coord_pO[2]  ]
                         
              v_H1O_2_cart = [ -coord_wH1[0] + coord_wO[0], \
                               -coord_wH1[1] + coord_wO[1], \
                               -coord_wH1[2] + coord_wO[2]  ]
                               
              # calculate vector H2->O1 and H2->O2
              v_H2O_1_cart = [ -coord_wH2[0] + coord_pO[0], \
                               -coord_wH2[1] + coord_pO[1], \
                               -coord_wH2[2] + coord_pO[2]  ]
                         
              v_H2O_2_cart = [ -coord_wH2[0] + coord_wO[0], \
                               -coord_wH2[1] + coord_wO[1], \
                               -coord_wH2[2] + coord_wO[2]  ]                  
                         
              # convert them to fractional coords           
              v_H1O_1_frac = tt.cart2frac( self.traj_full.lattice_vector_1 , \
                                           self.traj_full.lattice_vector_2 , \
                                           self.traj_full.lattice_vector_3 , \
                                           v_H1O_1_cart                     )      
              v_H1O_2_frac = tt.cart2frac( self.traj_full.lattice_vector_1 , \
                                           self.traj_full.lattice_vector_2 , \
                                           self.traj_full.lattice_vector_3 , \
                                           v_H1O_2_cart                     ) 
              v_H2O_1_frac = tt.cart2frac( self.traj_full.lattice_vector_1 , \
                                           self.traj_full.lattice_vector_2 , \
                                           self.traj_full.lattice_vector_3 , \
                                           v_H2O_1_cart                     )      
              v_H2O_2_frac = tt.cart2frac( self.traj_full.lattice_vector_1 , \
                                           self.traj_full.lattice_vector_2 , \
                                           self.traj_full.lattice_vector_3 , \
                                           v_H2O_2_cart                     )
                                           
              # apply pbc
              v_H1O_1_frac[:] = [ d-round(d) for d in v_H1O_1_frac ]
              v_H1O_2_frac[:] = [ d-round(d) for d in v_H1O_2_frac ]
              v_H2O_1_frac[:] = [ d-round(d) for d in v_H2O_1_frac ]
              v_H2O_2_frac[:] = [ d-round(d) for d in v_H2O_2_frac ]
              
              # convert them to back to cartesian coords           
              v_H1O_1_cart = tt.frac2cart( self.traj_full.lattice_vector_1 , \
                                           self.traj_full.lattice_vector_2 , \
                                           self.traj_full.lattice_vector_3 , \
                                           v_H1O_1_frac                     )      
              v_H1O_2_cart = tt.frac2cart( self.traj_full.lattice_vector_1 , \
                                           self.traj_full.lattice_vector_2 , \
                                           self.traj_full.lattice_vector_3 , \
                                           v_H1O_2_frac                     )
              v_H2O_1_cart = tt.frac2cart( self.traj_full.lattice_vector_1 , \
                                           self.traj_full.lattice_vector_2 , \
                                           self.traj_full.lattice_vector_3 , \
                                           v_H2O_1_frac                     )      
              v_H2O_2_cart = tt.frac2cart( self.traj_full.lattice_vector_1 , \
                                           self.traj_full.lattice_vector_2 , \
                                           self.traj_full.lattice_vector_3 , \
                                           v_H2O_2_frac                     )
                                           
              # calculate angle
              alpha1 = np.arccos( np.dot (v_H1O_1_cart, v_H1O_2_cart) / \
                                  (np.linalg.norm(v_H1O_1_cart)*np.linalg.norm(v_H1O_2_cart)) )
              alpha2 = np.arccos( np.dot (v_H2O_1_cart, v_H2O_2_cart) / \
                                  (np.linalg.norm(v_H2O_1_cart)*np.linalg.norm(v_H2O_2_cart)) )                   
              
              # convert to degree
              alpha1 = alpha1*360.0/(2*np.pi)
              alpha2 = alpha2*360.0/(2*np.pi)
              
              # found a Hbond?
              # atom should only be able to donate one H bond to the same member of list_index_partners!
              #  it can however donate 2 Hbonds in total, but to DIFFERENT members of list_index_partners
              #  only one entry is registered
              if alpha1 > tol_angle and alpha2 > tol_angle:
                print "[ERROR]. Frame: ", idx_frame, "Atom ", idx, " seems to donate two hydrogen bonds to the same atom", partner[0]
                sys.exit()
              elif alpha1 > tol_angle:
                list_donors_tmp.append([partner[0][0], idx+1])
              elif alpha2 > tol_angle:
                list_donors_tmp.append([partner[0][0], idx+2]) 
              
            #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
            #                            consistency check                                          #  
            #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
            #  atom should only either accept OR donate a hydrogen bond to the same member in list_index_partners
            #  it can however be involved in Hbonds with DIFFERENT members in list_index_partners 
            #  the latter case will be handled later
            if len(list_acceptors_tmp) > 0 and len(list_donors_tmp) > 0:
              print "[WARNING]. Frame: ", idx_frame, "Atom ", idx, " seems to be involved in two hydrogen bonds to the same atom ", partner[0]
              print "           list_acceptors_tmp: ", list_acceptors_tmp
              print "           list_donors_tmp   : ", list_donors_tmp
              print "           --> IGNORING THIS CONFIGURATION <--"
              print ""
              list_acceptors_tmp=[]
              list_donors_tmp=[]
              #sys.exit()
            
            # also, no member of list_index_partners should be involved in more than ONLY 1 Hbond
            if len(list_acceptors_tmp) > 1 or len(list_donors_tmp) > 1:
              print "[WARNING]. Frame: ", idx_frame, "Atom ", idx, " seems to be involved in two hydrogen bonds to the same atom ", partner[0]
              print "           list_acceptors_tmp: ", list_acceptors_tmp
              print "           list_donors_tmp   : ", list_donors_tmp
              print "           --> IGNORING THIS CONFIGURATION <--"
              print ""
              list_acceptors_tmp=[]
              list_donors_tmp=[]
              #sys.exit()
              
            # if a Hbond was found, store info in list_donor / list_acceptor
            if len(list_acceptors_tmp) == 1:
              list_acceptors.extend(list_acceptors_tmp)
              flag_Hbond_found = True
            elif len(list_donors_tmp) == 1:
              list_donors.extend(list_donors_tmp)
              flag_Hbond_found = True
          
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
          #                            add atom to selector if Hbond formed                       #  
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
          #  atom can have several H bonds to members in list_index_partners
          #  e.g. accept Hbond from member1 and donate Hbond to member2
          # atom is added to selector in any case
          # but in addition to that, info about idx is stored as well for analyser
          if flag_Hbond_found == True:

            # ***** add oxygen *****
            tmp_atom=atom
            bonding=[]

            # bonding partners from which atom accepts Hbond
            for acceptor in list_acceptors:
              #print acceptor
              #              flag  accept from  H involved
              bonding.append(["A", acceptor[0], acceptor[1]])
            for donor in list_donors:
              #print donor
              #              flag  donate to   H involved
              bonding.append(["D", donor[0],   donor[1]])


            # ***** project atom inside unit cell? *****
            if flag_project_inside == True:
              #  to do so, we convert x and y into fractional coordinates
              #  if one of them happens to be out of box --> PBC
              coords_frac = tt.cart2frac( self.traj_full.lattice_vector_1 , \
                                          self.traj_full.lattice_vector_2 , \
                                          self.traj_full.lattice_vector_3 , \
                                          [atom[1], atom[2], atom[3]] )

              # pbc
              # the int() makes sure that values end up in interval [0,1]
              # usually: pbc with round!
              coords_frac[0] = coords_frac[0] - int(coords_frac[0])
              coords_frac[1] = coords_frac[1] - int(coords_frac[1])

              if coords_frac[0] < 0:
                coords_frac[0]+=1
              if coords_frac[1] < 0:
                coords_frac[1]+=1

              # convert back to cartesian coords
              coords_cart = tt.frac2cart( self.traj_full.lattice_vector_1 , \
                                          self.traj_full.lattice_vector_2 , \
                                          self.traj_full.lattice_vector_3 , \
                                          coords_frac )
              
              atom[1] = coords_cart[0]
              atom[2] = coords_cart[1]
              atom[3] = coords_cart[2]

              if atom[1]<0 or atom[2]<0 or atom[3]<0:
                print "coords < 0.0"
                print idx_frame, atom
                print coords_frac
                sys.exit()
            
            # ***** add extra info *****
            extra_info = [ idx ]
            extra_info.extend(bonding)

            tmp_atom.append(extra_info)
            list_atom_select.append(tmp_atom)

            # ***** add hydrogen *****
            if flag_whole_water == True:
              list_atom_select.append(frame[idx+1] + [[idx+1]])
              list_atom_select.append(frame[idx+2] + [[idx+2]])
          
      # --- ATOM LOOP DONE ---
      # add all atoms from this frame as a new element in traj_sel
      self.traj_sel.list_atoms.append(list_atom_select)

    # --- FRAME LOOP DONE ---
    # set title
    self.traj_sel.title = ""
    
    # set n_atoms and n_frames
    # n_atoms ... biggest number of atoms in a frame
    self.traj_sel.n_atoms  = len(max(self.traj_sel.list_atoms,key=len))
    self.traj_sel.n_frames = len(self.traj_sel.list_atoms)

    if self.traj_sel.n_frames != self.traj_full.n_frames:
      print "[ERROR]. Hbond SELECTOR: Something went wrong, n_frames changed!"
      sys.exit()

    # set lattice_vectors
    self.traj_sel.lattice_vector_1 = self.traj_full.lattice_vector_1
    self.traj_sel.lattice_vector_2 = self.traj_full.lattice_vector_2
    self.traj_sel.lattice_vector_3 = self.traj_full.lattice_vector_3
  
  ##########################################################################################
  #                            NOT_BONDING                                                 # 
  ##########################################################################################  
  def not_bonding( self,                          \
                   index_min,                     \
                   index_max,                     \
                   h_min,                         \
                   h_max,                         \
                   tol_length,                    \
                   tol_angle,                     \
                   tol_distance,                  \
                   list_index_ref_z=[],           \
                   list_index_Hbond_partners=[],  \
                   list_index_coord_partners=[],  \
                   list_name_select=["O"],        \
                   flag_whole_water=False,        \
                   flag_include_partners=False,   \
                   flag_project_inside=True):
    """
    select all atoms that are not Hbonding or coordinating the selected motives
    
    index_min ... minimum index of relevant analyt
    index_max ... maximum index of relevant analyt
    h_min            ... minimum height of water to be considered for analysis
    h_max            ... maximum height of water to be considered for analysis

    list_index_ref_z ... list of reference points to determine slab height
    
    tol_length ... maximum acceptable Hbond (O-O) length
    tol_angle ... minimum acceptable Hbond (O-H-O) angle
    tol_distance ... maximum acceptable coordination bond length
    
    list_index_Hbond_partners ... list with all bonding partners that are considered
                                  the format of this list is the following
                                  [[A1], [D1O, D1H1]], [[A2], [D2O, D2H]]
                                  for more details look at H-BONDING routine
    list_index_coord_partners ... list with all coordination partners that are considered
    
    list_name_select    ... atom types which are selected
   
    flag_whole_water      ... determines if whole water molecules is selected or not
    flag_include_partners ... include bonding partners in selection?
    flag_project_inside   ... project atoms inside simulation box?
    """
    self.check_traj()   # check if traj makes sense

    # makros:
    NAME = 0
    X    = 1
    Y    = 2
    Z    = 3

    # figure out if arguments OK
    if len(list_index_Hbond_partners)==0 and len(list_index_coord_partners)==0:
      print "[ERROR]. No potential bonding partners specified"
      sys.exit()

    if len(list_index_Hbond_partners)>0:
      for partner in list_index_Hbond_partners:
        if len(partner) != 2 or (len(partner[0]) != 1 and len(partner[1]) != 2):
          print "[ERROR]. The partner list does not make sense."
          print "         Length of individual members of list wrong"
          print "         Read the code..."
          print partner[0]
          print partner[1]
          sys.exit()
        if len(partner[1])>0:
          if partner[0][0] != partner[1][0] :
            print "[ERROR]. The partner list does not make sense."
            print "         Indizes wrong"
            print "         Read the code..."
            print partner[0][0], partner[1][0]
            sys.exit()
    
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
    #                            LOOP OVER FRAMES                                           #  
    #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
    for idx_frame, frame in enumerate(self.traj_full.list_atoms):
      if int(self.traj_full.n_frames*0.25)==idx_frame:
        print "         25 % done"
      if int(self.traj_full.n_frames*0.50)==idx_frame:
        print "         50 % done"
      if int(self.traj_full.n_frames*0.75)==idx_frame:
        print "         75 % done   "
   
      list_atom_select = []
      
      # calculate reference height for this frame
      ref_z = 0.0
      for idx_reference in list_index_ref_z:
        ref_z += frame[idx_reference][3]/len(list_index_ref_z)

      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            include bonding partners?                                  #  
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      if flag_include_partners == True:
        if len(list_index_Hbond_partners)>0:
          for partner in list_index_Hbond_partners:
            # find all atoms that belong to partner
            # also mark them as bonding partner for analyser
            # partner as ONLY  Hbond acceptor
            if len(partner[0]) == 1 and len(partner[1])==0:
              tmp_atom=frame[partner[0][0]]
              tmp_atom[0]+="_bp"
              tmp_atom.append([partner[0][0]])
              list_atom_select.append(tmp_atom)
            # partner as Hbond donor (and potentially acceptor too)
            #  the second element is H
            elif len(partner[1]) == 2:
              # append O
              tmp_atom=frame[partner[1][0]]
              tmp_atom[0]+="_bp"
              tmp_atom.append([partner[1][0]])
              list_atom_select.append(tmp_atom)
              
              # append H
              tmp_atom=frame[partner[1][1]]
              tmp_atom[0]+="_bp"
              tmp_atom.append([partner[1][1]])
              list_atom_select.append(tmp_atom)
            else:
              print "[ERROR]. not_bonding-selector. The Hbond partner-list does not make sense!"
              sys.exit()
      
        # include bonding_partners if necessary
        if len(list_index_coord_partners)>1:
          for partner in list_index_coord_partners:
              tmp_atom=frame[partner]
              tmp_atom[0]+="_bp"
              tmp_atom.append([partner])
              list_atom_select.append(tmp_atom)

      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
      #                            LOOP OVER ATOMS IN FRAME                                   #  
      #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
      for idx,atom in enumerate(frame):
        # check if index is OK
        if idx>= index_min and idx<= index_max:

          # check if atom in correct height interval
          if atom[Z] < ref_z+h_min or atom[Z] > ref_z+h_max:
            continue
          
          # check if atom type not relevant for analysis:
          if atom[0] not in list_name_select:
            continue
          
          list_donors       = []  # complete list of Hbond to atom
          list_acceptors    = []  # complete list of Hbond from atom
          list_coordinators = []  # complete list of coordination partners to atom
          
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
          #                            FIND COORDINATION PARTNERS                                 #  
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
          if len(list_index_coord_partners)>0:
            for partner in list_index_coord_partners:
              p1 = [ atom[X], atom[Y], atom[Z] ]
              p2 = [ frame[partner][X] , \
                     frame[partner][Y] , \
                     frame[partner][Z] ]
              
              d = pbc.pbc_distance ( p1 = p1,                                        \
                                     p2 = p2,                                        \
                                     lattice_v_1 = self.traj_full.lattice_vector_1 , \
                                     lattice_v_2 = self.traj_full.lattice_vector_2 , \
                                     lattice_v_3 = self.traj_full.lattice_vector_3 )
              if d<=tol_distance:
                list_coordinators.append(partner)
                break # we only need to find first partner to know that there are coodination partners

          if len(list_coordinators)>0:
            # if we found that atom forms a coordinative bond --> move on to next atom
            continue

          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
          #                            FIND HBONDS                                                #  
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
          if len(list_index_Hbond_partners)>0:
            for partner in list_index_Hbond_partners:
              
              p1 = [ atom[X], atom[Y], atom[Z] ]
              p2 = []
              
              if len(partner[0])==1:
                p2 = [ frame[partner[0][0]][1], frame[partner[0][0]][2], frame[partner[0][0]][3] ] 
              else:
                p2 = [ frame[partner[1][0]][1], frame[partner[1][0]][2], frame[partner[1][0]][3] ] 
              
              d = pbc.pbc_distance ( p1 = p1,                                        \
                                     p2 = p2,                                        \
                                     lattice_v_1 = self.traj_full.lattice_vector_1 , \
                                     lattice_v_2 = self.traj_full.lattice_vector_2 , \
                                     lattice_v_3 = self.traj_full.lattice_vector_3 )
              
              if d<=tol_length:
                #print "potential Hbond found", idx, partner
                #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
                #                            check if Hbond is accepted by atom                         #  
                #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
                if len(partner[1])==2:
                  # get coordinates for O and H from partner
                  coord_pO = [ frame[partner[1][0]][1], frame[partner[1][0]][2], frame[partner[1][0]][3] ]
                  coord_pH = [ frame[partner[1][1]][1], frame[partner[1][1]][2], frame[partner[1][1]][3] ]
                  
                  alpha = pbc.pbc_angle ( p0=coord_pH, \
                                          p1=p1,       \
                                          p2=coord_pO, \
                                          lattice_v_1 = self.traj_full.lattice_vector_1 , \
                                          lattice_v_2 = self.traj_full.lattice_vector_2 , \
                                          lattice_v_3 = self.traj_full.lattice_vector_3 )
                  
                  # found a Hbond?
                  if alpha > tol_angle:
                    list_acceptors.append([partner[1][0], partner[1][1]])
                
                if len(list_acceptors)>0:
                  break

                #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
                #                            check if Hbond is donated from atom                        #  
                #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
                if len(partner[0])==1:
                  # get coordinates for O from partner
                  coord_pO = [ frame[partner[0][0]][1], frame[partner[0][0]][2], frame[partner[0][0]][3] ]  
                  
                  # get coordinates for water molecule
                  # assumption: water is stored as O H H in trajectory
                  coord_wO  = [ atom[1], atom[2], atom[3] ]
                  coord_wH1 = [ frame[idx+1][1], frame[idx+1][2], frame[idx+1][3] ]
                  coord_wH2 = [ frame[idx+2][1], frame[idx+2][2], frame[idx+2][3] ] 
                  
                  alpha1 = pbc.pbc_angle ( p0=coord_wH1, \
                                           p1=coord_wO,  \
                                           p2=coord_pO,  \
                                           lattice_v_1 = self.traj_full.lattice_vector_1 , \
                                           lattice_v_2 = self.traj_full.lattice_vector_2 , \
                                           lattice_v_3 = self.traj_full.lattice_vector_3 )
                  
                  alpha2 = pbc.pbc_angle ( p0=coord_wH2, \
                                           p1=coord_wO,  \
                                           p2=coord_pO,  \
                                           lattice_v_1 = self.traj_full.lattice_vector_1 , \
                                           lattice_v_2 = self.traj_full.lattice_vector_2 , \
                                           lattice_v_3 = self.traj_full.lattice_vector_3 )

                  # found a Hbond?
                  if alpha1 > tol_angle and alpha2 > tol_angle:
                    print "[ERROR]. Frame: ", idx_frame, "Atom ", idx, " seems to donate two hydrogen bonds to the same atom", partner[0]
                    sys.exit()
                  elif alpha1 > tol_angle or alpha2>tol_angle:
                    list_donors.append([partner[0][0], idx+1])

                if len(list_donors)>0:
                  break
            # --- END PARTNER LOOP ---
              
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
          #                            Hbonds found?                                              #  
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
          if len(list_acceptors)>0 or len(list_donors)>0:
            continue

          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#
          #                            no bonding to specified motives found                      #  
          #---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---# 
          # if there was a bond (either coordinative or Hbond), than we would not get this far in loop
          # if we are here however --> we found an atom that does not bind with specified indices

          # ***** project atom inside unit cell? *****
          if flag_project_inside == True:
            #  to do so, we convert x and y into fractional coordinates
            #  if one of them happens to be out of box --> PBC
            coords_frac = tt.cart2frac( self.traj_full.lattice_vector_1 , \
                                        self.traj_full.lattice_vector_2 , \
                                        self.traj_full.lattice_vector_3 , \
                                        [atom[1], atom[2], atom[3]] )

            # pbc
            # the int() makes sure that values end up in interval [0,1]
            # usually: pbc with round!
            coords_frac[0] = coords_frac[0] - int(coords_frac[0])
            coords_frac[1] = coords_frac[1] - int(coords_frac[1])

            if coords_frac[0] < 0:
              coords_frac[0]+=1
            if coords_frac[1] < 0:
              coords_frac[1]+=1

            # convert back to cartesian coords
            coords_cart = tt.frac2cart( self.traj_full.lattice_vector_1 , \
                                        self.traj_full.lattice_vector_2 , \
                                        self.traj_full.lattice_vector_3 , \
                                        coords_frac )
            
            atom[1] = coords_cart[0]
            atom[2] = coords_cart[1]
            atom[3] = coords_cart[2]

            if atom[1]<0 or atom[2]<0 or atom[3]<0:
              print "coords < 0.0"
              print idx_frame, atom
              print coords_frac
              sys.exit()
          
          # ***** add extra info *****
          list_atom_select.append(atom + [[idx]])

          # ***** add hydrogen *****
          if flag_whole_water == True:
            list_atom_select.append(frame[idx+1] + [[idx+1]])
            list_atom_select.append(frame[idx+2] + [[idx+2]])
          
      # --- ATOM LOOP DONE ---
      # add all atoms from this frame as a new element in traj_sel
      self.traj_sel.list_atoms.append(list_atom_select)

    # --- FRAME LOOP DONE ---
    # set title
    self.traj_sel.title = ""
    
    # set n_atoms and n_frames
    # n_atoms ... biggest number of atoms in a frame
    self.traj_sel.n_atoms  = len(max(self.traj_sel.list_atoms,key=len))
    self.traj_sel.n_frames = len(self.traj_sel.list_atoms)

    if self.traj_sel.n_frames != self.traj_full.n_frames:
      print "[ERROR]. Hbond SELECTOR: Something went wrong, n_frames changed!"
      sys.exit()

    # set lattice_vectors
    self.traj_sel.lattice_vector_1 = self.traj_full.lattice_vector_1
    self.traj_sel.lattice_vector_2 = self.traj_full.lattice_vector_2
    self.traj_sel.lattice_vector_3 = self.traj_full.lattice_vector_3
