import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

graphite = Graphene(config,forcefield)
sim = Crystal(graphite,config,forcefield,[30,20,1])

j = 0
for i in [0,2,3,4,4,5,5,4,4,3,2,0,1,2,3,2,1,0,1,0]:
    j += 1
    motif = Hexagon_Graphene(config,forcefield,5+i*2.44)
    next_layer = Crystal(motif,config,forcefield,[1,1,1])

    next_layer.coords = next_layer.coords + np.array((
        5 * 2*(3**0.5) * config[forcefield]['CC'],
        24*config[forcefield]['CC'],
        j*config[forcefield]['layer_gap']))
    sim  = Combine(sim,next_layer)
#bulk.box_dimensions[2] = bulk.box_dimensions[2] + 20 - config[forcefield]['layer_gap']



name = 'graphene'
output = Writer(sim,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
