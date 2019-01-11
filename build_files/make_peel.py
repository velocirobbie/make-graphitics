import yaml
from scripts.molecules import *
from scripts import *
import numpy as np
import math

config = yaml.load(open('config.yaml'))
forcefield = 'GraFF_5'
graphite = Graphite(config,forcefield)
bulk = Crystal(graphite,config,forcefield,[21,12,1])
bulk.coords = bulk.coords + np.array((
    0,
    0,
    1*config[forcefield]['layer_gap']))


molecule1 = Hexagon_Graphene(config,forcefield,15)
flake1 = Crystal(molecule1,config,forcefield,[1,1,1])

flake1.coords = flake1.coords + np.array((
    10 * 2*math.cos(math.pi/6) * config[forcefield]['CC'],
    6 * 3 *config[forcefield]['CC'],
    3*config[forcefield]['layer_gap']))

sim = Combine(bulk,flake1)
#sim = Combine(sim,top)
output = Writer(sim,'flake on graphite')
output.write_xyz('graphene'+str(1)+'.xyz')
output.write_lammps('data.flake'+str(1))
