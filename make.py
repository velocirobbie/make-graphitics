import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

graphite = Graphite(config,forcefield)
bulk = Crystal(graphite,config,forcefield,[20,12,4])

molecule = Hexagon_Graphene(config,forcefield,10)
flake = Crystal(molecule,config,forcefield,[1,1,1])

flake.coords = flake.coords + np.array((
    5 * 2*(3**0.5) * config[forcefield]['CC'],
    11*config[forcefield]['CC'],
    8*config[forcefield]['layer_gap']))

sim = Combine(bulk,flake)
output = Writer(sim,'flake on graphite')
#output.write_xyz('graphene.xyz')
#output.write_lammps('data.flake')
