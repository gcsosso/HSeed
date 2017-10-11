#!/usr/bin/python

import numpy as np
import math
import copy
import sys
import transform_coordinates as tc

class lattice_building_block:
  """
  Class that contains info about the building blocks available to generate the lattices

  Info contained:
   (*) unit cell of building block (a,b,c)
   (*) lattice positions in building block
  """

  def __init__(self):
    self.a = np.zeros(3) # unit cell vector of building block
    self.b = np.zeros(3) # unit cell vector of building block
    self.c = np.zeros(3) # unit cell vector of building block

    self.list_latticepoints = [] # array containing all the lattice points of building block

  def set_a(self, x, y, z):
    self.a = np.array([x, y, z])

  def set_b(self, x, y, z):
    self.b = np.array([x, y, z])

  def set_c(self, x, y, z):
    self.c = np.array([x, y, z])

  def add_latticepoint(self, x, y, z):
    self.list_latticepoints.append(np.array([x, y, z]))

  def del_latticepoints(self):
    self.list_latticepoints = []

  ########################################################################################
  #                          multiply building blocks                                    #
  ########################################################################################
  def multiply_building_block(self, multiply_x, multiply_y):
    multiply_x = int(multiply_x)
    multiply_y = int(multiply_y)

    self.a = self.a * multiply_x
    self.b = self.b * multiply_y

    # -----------------
    # - multiply in x -
    # -----------------
    list_additional_boxes = []
    for i in range(1, multiply_x):
      list_latticepoints_new_box = copy.deepcopy(self.list_latticepoints)
      for latticepoint in list_latticepoints_new_box:
        latticepoint[0] += i

      list_additional_boxes.append(list_latticepoints_new_box)

    # --- merge boxes ---
    for box in list_additional_boxes:
      self.list_latticepoints += box

    # --- compute correct fractional coords in supercell ---
    for latticepoint in self.list_latticepoints:
      latticepoint[0] = latticepoint[0]/multiply_x
    

    # -----------------
    # - multiply in y -
    # -----------------
    list_additional_boxes = []
    for i in range(1, multiply_y):
      list_latticepoints_new_box = copy.deepcopy(self.list_latticepoints)
      for latticepoint in list_latticepoints_new_box:
        latticepoint[1] += i

      list_additional_boxes.append(list_latticepoints_new_box)
    
    # --- merge boxes
    for box in list_additional_boxes:
      self.list_latticepoints += box

    # --- compute correct fractional coords in supercell ---
    for latticepoint in self.list_latticepoints:
      latticepoint[1] = latticepoint[1]/multiply_y

  ########################################################################################
  #                          convert latticecoords to slabcoords                         #
  ########################################################################################
  def lattice2slab(self, slab_a, slab_b, slab_c, flag_scale_lattice):
    """
    convert fractional lattice points (lattice unit cell) to fractional lattice points (slab unit cell)

    lattice[frac-lattice] --> lattice[cartesian] --> lattice[frac-slab]

    if specified (flag_scale_lattice == True), the lattice is scaled so that it covers the slab uniformly
    """

    # in old version this was not commented out
    # it was necessary for GROMACS RSS to convert lattice to AA (they came in nm)
    # --> now this is done in the main code, because LAMMPS should not do this!
    # [OLD VERSION - FOR GROMACS] slab_a = [ i*10 for i in slab_a]
    # [OLD VERSION - FOR GROMACS] slab_b = [ i*10 for i in slab_b]
    # [OLD VERSION - FOR GROMACS] slab_c = [ i*10 for i in slab_c]

    # --------------------------
    # --- rescale (optional) ---
    # --------------------------
    if (flag_scale_lattice == True):
      """
      for this to work, a should be oriented along x and b along y (at least approximately)
      also: this is intrinsically flawed if the tilt (eg a[1] and b[0] between slab and lattice
      are very different
      and there is nothing which we can do about it if we want to keep the size of the unit cell fixed...


      the best we can do: do the procedure twice, once with scaled and once with unscaled lattice
      """

      quotient_a = self.a[0] / slab_a[0] # scale x of a axis 
      self.a = self.a / quotient_a
      
      quotient_b = self.b[1] / slab_b[1] # scale y of b axis
      self.b = self.b / quotient_b
      
    # ------------------------------
    # --- convert to cart coords ---
    # ------------------------------
    latticepoints_cart = []
    for latticepoint in self.list_latticepoints:
      cart = tc.frac2cart(self.a.tolist(), self.b.tolist(), self.c.tolist(), latticepoint.tolist())
      latticepoints_cart.append(cart)

    # -----------------------------------------------
    # --- convert to frac coords (slab unit cell) ---
    # -----------------------------------------------
    latticepoints_frac = []
    for latticepoint in latticepoints_cart:

      frac = tc.cart2frac(slab_a, slab_b, slab_c, latticepoint)
      latticepoints_frac.append(frac)

    # ---------------------------
    # --- reassign latticepoints ---
    # ---------------------------
    self.set_a(slab_a[0], slab_a[1], slab_a[2])
    self.set_b(slab_b[0], slab_b[1], slab_b[2])
    self.set_c(slab_c[0], slab_c[1], slab_c[2])
    
    self.del_latticepoints()
    for latticepoint in latticepoints_frac:
      self.add_latticepoint(latticepoint[0], latticepoint[1], latticepoint[2])


  def show_building_block(self):
    print self.a
    print self.b
    print self.c
    print "number of latticepoints: ", len(self.list_latticepoints)
    for latticepoint in self.list_latticepoints:
      print "\t", latticepoint
    print "--------------------\n"



def construct_lattice(lattice, multiplier_x, multiplier_y ):
  """
  Construct a lattice in fractional coordinates for the adsorbate

  (1) Choose one of the building blocks for lattice:
         --> rectangular lattice
         --> hexagonal lattice
         --> lines
         --> pentagons
  (2) multiply lattice so that it covers unit cell
  (3) optional: stretch/compress the lattice so the lattice unit cell matches up with the substrate unit cell
  (3a)  compute latticepoints in fractional coordinates relative to substrate unit cell
  """

  lattice.multiply_building_block(multiplier_x, multiplier_y)

  # lattice.show_building_block()
