import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

motif = Rectangle_Graphene(config,forcefield,50,50)
flake = Crystal(motif,config,forcefield,[1,1,1])
Oxidiser(flake, ratio=2.5, video=20, new_island_freq=3e15, method='rf')

Parameterise(flake,flake.vdw_defs)

name = 'graphene'
output = Writer(flake,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
