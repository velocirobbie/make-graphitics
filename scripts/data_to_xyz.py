import yaml
from scripts.molecules import *
from scripts import *
import numpy as np
import sys

sim = ReadLammpsData(sys.argv[1])


name = 'out'
output = Writer(sim,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
