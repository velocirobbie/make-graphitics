import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

motif = Rectangle_Graphene(config,forcefield,50,50)
flake = Crystal(motif,config,forcefield,[1,1,1])

oxidiser = Oxidiser(ratio=2.5, video_xyz=20, new_island_freq=1e14, method='rf')
flake = oxidiser.react(flake)

Parameterise(flake,flake.vdw_defs)

name = 'graphene'
output = Writer(flake,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
