#!/usr/bin/python

import ase
from ase import Atoms
from ase.lattice.surface import surface
from ase.io import *


ice = read("0000_Ic_SC_6x4x4.POSCAR")

sf = surface(ice, (-3,8,2), 20)
sf.center(vacuum=10, axis=2)

write("test.POSCAR", sf)
