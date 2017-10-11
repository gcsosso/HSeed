import os.path
import math
import sys
import numpy as np

import transform_coordinates as tt

"""
calculate properties taking into account periodic boundary conditions

(+) pbc_distance
(+) pbc_vector
(+) pbc_angle

"""

##########################################################################################
#                            pbc_distance                                                # 
##########################################################################################  
def pbc_distance ( p1 = [],          \
                   p2 = [],          \
                   lattice_v_1 = [], \
                   lattice_v_2 = [], \
                   lattice_v_3 = [] ):
  # sanity check
  if len(p1) != 3 or len(p2) != 3:
    print "[ERROR]. pbc-distance"
    print "         wrong arguments (particles)"
    sys.exit()
  if len(lattice_v_1) != 3 or len(lattice_v_2) != 3 or len(lattice_v_3) != 3:
    print "[ERROR]. pbc-distance"
    print "         wrong arguments (lattice_v)"
    sys.exit()

  # calculate cartesian distance vector between partner and atom
  dist_v_cart = [ p1[0] - p2[0] , \
                  p1[1] - p2[1] , \
                  p1[2] - p2[2] ]
  
  # convert distance vector to fractional coords
  dist_v_frac = tt.cart2frac( lattice_v_1 , \
                              lattice_v_2 , \
                              lattice_v_3 , \
                              dist_v_cart )
                              
  # apply pbc
  dist_v_frac[:] = [ d-round(d) for d in dist_v_frac ]
  
  # convert back to cartesian coords  
  dist_v_cart = tt.frac2cart( lattice_v_1 , \
                              lattice_v_2 , \
                              lattice_v_3 , \
                              dist_v_frac   ) 
                              
  # calculate distance
  d = math.sqrt(dist_v_cart[0]**2 + dist_v_cart[1]**2 + dist_v_cart[2]**2) 

  return d

##########################################################################################
#                            pbc_vector                                                  # 
##########################################################################################  
def pbc_vector ( p1 = [],          \
                 p2 = [],          \
                 lattice_v_1 = [], \
                 lattice_v_2 = [], \
                 lattice_v_3 = [] ):
  # sanity check
  if len(p1) != 3 or len(p2) != 3:
    print "[ERROR]. pbc-distance"
    print "         wrong arguments (particles)"
    sys.exit()
  if len(lattice_v_1) != 3 or len(lattice_v_2) != 3 or len(lattice_v_3) != 3:
    print "[ERROR]. pbc-distance"
    print "         wrong arguments (lattice_v)"
    sys.exit()

  # calculate cartesian distance vector between partner and atom
  dist_v_cart = [ p1[0] - p2[0] , \
                  p1[1] - p2[1] , \
                  p1[2] - p2[2] ]
  
  # convert distance vector to fractional coords
  dist_v_frac = tt.cart2frac( lattice_v_1 , \
                              lattice_v_2 , \
                              lattice_v_3 , \
                              dist_v_cart )
                              
  # apply pbc
  dist_v_frac[:] = [ d-round(d) for d in dist_v_frac ]
  
  # convert back to cartesian coords  
  dist_v_cart = tt.frac2cart( lattice_v_1 , \
                              lattice_v_2 , \
                              lattice_v_3 , \
                              dist_v_frac   ) 
                              
  return dist_v_cart

##########################################################################################
#                            pbc_angle                                                   # 
##########################################################################################  
def pbc_angle ( p0 = [],          \
                p1 = [],          \
                p2 = [],          \
                lattice_v_1 = [], \
                lattice_v_2 = [], \
                lattice_v_3 = [] ):
  """
  p0 ... the atom in between p1 and p2

  for a Hbond for example, p0 ... H, p1 and p2 ... O
  """
  
  # sanity check
  if len(p1) != 3 or len(p2) != 3:
    print "[ERROR]. pbc-distance"
    print "         wrong arguments (particles)"
    sys.exit()
  if len(lattice_v_1) != 3 or len(lattice_v_2) != 3 or len(lattice_v_3) != 3:
    print "[ERROR]. pbc-distance"
    print "         wrong arguments (lattice_v)"
    sys.exit()

  # calculate vector p0 -> p1 and p0 -> p2
  v_1_cart = [ -p0[0] + p1[0], \
               -p0[1] + p1[1], \
               -p0[2] + p1[2]  ]
             
  v_2_cart = [ -p0[0] + p2[0], \
               -p0[1] + p2[1], \
               -p0[2] + p2[2]  ] 
             
  # convert them to fractional coords           
  v_1_frac = tt.cart2frac( lattice_v_1 , \
                           lattice_v_2 , \
                           lattice_v_3 , \
                           v_1_cart      )      
  v_2_frac = tt.cart2frac( lattice_v_1 , \
                           lattice_v_2 , \
                           lattice_v_3 , \
                           v_2_cart      ) 
  
  # apply pbc
  v_1_frac[:] = [ d-round(d) for d in v_1_frac ]
  v_2_frac[:] = [ d-round(d) for d in v_2_frac ]
  
  # convert them to back to cartesian coords           
  v_1_cart = tt.frac2cart( lattice_v_1 , \
                           lattice_v_2 , \
                           lattice_v_3 , \
                           v_1_frac   )      
  v_2_cart = tt.frac2cart( lattice_v_1 , \
                           lattice_v_2 , \
                           lattice_v_3 , \
                           v_2_frac   )
  
  # calculate angle
  alpha = np.arccos( np.dot (v_1_cart, v_2_cart) / \
                     (np.linalg.norm(v_1_cart)*np.linalg.norm(v_2_cart)) )
                     
  # convert to degree
  alpha = alpha*360.0/(2*np.pi)

  return alpha
