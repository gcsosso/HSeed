#!/usr/bin/python

import ase
from ase import Atoms
from ase.lattice.surface import surface
from ase.io import *


ice = read("0000_Ih_bulk.POSCAR")

sf = surface(ice, (0,0,1), 1)
sf.center(vacuum=0, axis=2)

write("Ih_SF_001.POSCAR", sf)
