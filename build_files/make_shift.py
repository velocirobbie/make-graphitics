import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

#for gap in range(20,21): 
#for i in [1]:
# gap = 3.48
# config['OPLS']['layer_gap']=gap
# config['OPLS']['CC']=1.400

motif = Graphite(config,forcefield)
bulk = Crystal(motif,config,forcefield,[10,6,4])

a = Shifter(bulk)
a.z_shift(1.0,6.0,2.0)

#name = 'layers_'+str(gap)
#output = Writer(bulk,name)
#output.write_xyz(name+'.xyz')
#output.write_lammps(name+'.data')
