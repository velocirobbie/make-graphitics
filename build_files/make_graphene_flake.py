import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

R=25
motif = Hexagon_Graphene(config,forcefield,R)
flake = Crystal(motif,config,forcefield,[1,1,1])
vdw_defs = {1:90, 2:91}

Parameterise(flake,vdw_defs)

name = 'graphene'
output = Writer(flake,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
