#!/usr/bin/python

import numpy as np

def frac2cart ( lv_1, lv_2, lv_3, point ):
  """
  This function converts fractional to cartesian coordinates

  lattice_vector_1/2/3 ... each one is a list, eg. l_v_1 = [lv1_x, lv1_y, lv1_z ]
  point                ... also a list         eg. point = [ p_x, p_y, p_z ]
  """
  T = np.matrix ([ [ lv_1[0], lv_2[0], lv_3[0] ] ,\
                   [ lv_1[1], lv_2[1], lv_3[1] ] ,\
                   [ lv_1[2], lv_2[2], lv_3[2] ] ])
  v = np.matrix ([ [point[0]], [point[1]], [point[2]] ])

  # multiply them
  w = T*v

  # convert back to list
  list_w = w.tolist()
  list_w = [ list_w[0][0], list_w[1][0], list_w[2][0] ]
  
  return list_w

def cart2frac ( lv_1, lv_2, lv_3, point ):
  """
  This function converts cartesian to fractional coordinates

  lattice_vector_1/2/3 ... each one is a list, eg. l_v_1 = [lv1_x, lv1_y, lv1_z ]
  point                ... also a list         eg. point = [ p_x, p_y, p_z ]
  """
  T = np.matrix ([ [ lv_1[0], lv_2[0], lv_3[0] ] ,\
                   [ lv_1[1], lv_2[1], lv_3[1] ] ,\
                   [ lv_1[2], lv_2[2], lv_3[2] ] ])
  T_inv = T.I
  v = np.matrix ([ [point[0]], [point[1]], [point[2]] ])

  # multiply them
  w = T_inv*v

  # convert back to list
  list_w = w.tolist()
  list_w = [ list_w[0][0], list_w[1][0], list_w[2][0] ]
  
  return list_w 
