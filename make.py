import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

motif = Hexagon_Graphene(config,forcefield,15)
layer = Crystal(motif,config,forcefield,[1,1,1])

layer.coords = layer.coords + np.array((
        5 * 2*(3**0.5) * config[forcefield]['CC'],
        24*config[forcefield]['CC'],
        0))#*config[forcefield]['layer_gap']))
#bulk.box_dimensions[2] = bulk.box_dimensions[2] + 20 - config[forcefield]['layer_gap']



name = 'graphene'
output = Writer(layer,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
