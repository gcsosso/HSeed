"""
compute properties of dipole distribution
"""

import sys
import numpy as np
from pbc_computations import pbc_vector

def compute_total_dipole(list_atoms, unit_cell=[], water_model="TIP4P/Ice"):
  # makros
  NAME = 0
  X = 1
  Y = 2
  Z = 3
  Q = 4
  
  CONV_eAA_to_D=4.80320 # conversion from eAngstrom to Debye
  
  if len(list_atoms[0]) != 4 or len(list_atoms[-1]) != 4:
    print "[ERROR]. The list_atoms object does not match the expected format [compute_total_dipole()]"
    print "           First element:", list_atoms[0]
    print "           Last element:",  list_atoms[-1]
    sys.exit()

  if water_model != "TIP4P/Ice":
    print "[ERROR]. At the moment only TIP4P/Ice is supported"
    sys.exit()

  total_dipole = 0.0
  for idx, atom in enumerate(list_atoms):
    if atom[NAME] == "MW":
      if list_atoms[idx-2][NAME] != "HW1" and list_atoms[idx-1][NAME] != "HW2":
        print "[ERROR]. The order of list_atoms does not match the expectation (OW, HW1, HW2, MW, ...)"
        sys.exit()

      # calculate MW -> HW1 and MW -> HW2 without pbc
      if len(unit_cell)==0:
        v_1 = np.array( [list_atoms[idx-2][X] - atom[X], \
                         list_atoms[idx-2][Y] - atom[Y], \
                         list_atoms[idx-2][Z] - atom[Z] ])*10

        v_2 = np.array( [list_atoms[idx-1][X] - atom[X], \
                         list_atoms[idx-1][Y] - atom[Y], \
                         list_atoms[idx-1][Z] - atom[Z] ])*10

        d_1 = np.linalg.norm(v_1)
        d_2 = np.linalg.norm(v_2)

        if d_1 > 2.0 or d_2 > 2.0:
          print "[ERROR]. No pbc supported [compute_total_dipole()] if unit cell dimention not specified"
          print "           ", d_1, d_2
          sys.exit()

        total_dipole += (v_1 + v_2)*0.5897
      # with pbc
      else:
        v_1 = np.array(pbc_vector([list_atoms[idx-2][X], list_atoms[idx-2][Y], list_atoms[idx-2][Z]], \
                                  [atom[X], atom[Y], atom[Z]], \
                                  unit_cell[0], \
                                  unit_cell[1], \
                                  unit_cell[2]))*10.0

        v_2 = np.array(pbc_vector([list_atoms[idx-1][X], list_atoms[idx-1][Y], list_atoms[idx-1][Z]], \
                                  [atom[X], atom[Y], atom[Z]], \
                                  unit_cell[0], \
                                  unit_cell[1], \
                                  unit_cell[2]))*10.0

        d_1 = np.linalg.norm(v_1)
        d_2 = np.linalg.norm(v_2)

        if d_1 > 2.0 or d_2 > 2.0:
          print "[ERROR]. pbc does not make sense [compute_total_dipole(), with PBC]. After applying pbc there should be no distances > 2.0 AA"
          print "           ", d_1, d_2
          sys.exit()

        total_dipole += (v_1 + v_2)*0.5897
  return total_dipole, np.linalg.norm(total_dipole)*CONV_eAA_to_D
