import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

for gap in range(20,21): 
 config['OPLS']['layer_gap']=gap

 motif = Graphene(config,forcefield)
 bulk = Crystal(motif,config,forcefield,[7,4,1])

 name = 'layers_'+str(gap)
 output = Writer(bulk,name)
 output.write_xyz(name+'.xyz')
 output.write_lammps(name+'.data')
