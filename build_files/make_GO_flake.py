import yaml
from scripts.molecules import *
from scripts import *
import numpy as np
import json

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'
flake_radius = 25
layout = [1,1,1] # make a 1x1x1 array of flakes

motif = Hexagon_Graphene(config,forcefield,flake_radius)
flake = Crystal(motif,config,forcefield,layout)
Oxidiser(flake, ratio=2.5, video=20, new_island_freq=1e16, method='rf')
name = 'GO_flake'
output = Writer(flake,name)
output.write_xyz(name+'.xyz')

Parameterise(flake,flake.vdw_defs)

name = 'GO_flake'
output = Writer(flake,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
