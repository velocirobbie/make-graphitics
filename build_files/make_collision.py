import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'GraFF_5'

graphite = Graphene(config,forcefield)
bulk = Crystal(graphite,config,forcefield,[120,37,1])


molecule1 = Hexagon_Graphene(config,forcefield,50)
flake1 = Crystal(molecule1,config,forcefield,[1,1,1])
molecule2 = Hexagon_Graphene(config,forcefield,50)
flake2 = Crystal(molecule2,config,forcefield,[1,1,1])
#make flake carbons different to bulk
#for atom in range(molecule.natoms):
#    if flake.atom_labels[atom] == 1:
#        flake.atom_labels[atom] = 3

flake1.coords = flake1.coords + np.array((
    15 * 2*(3**0.5) * config[forcefield]['CC'],
    54*config[forcefield]['CC'],
    config[forcefield]['layer_gap']+3.7))
flake2.coords = flake2.coords + np.array((
    45* 2*(3**0.5) * config[forcefield]['CC'],
    54*config[forcefield]['CC'],
    config[forcefield]['layer_gap']+3.7))
bulk.coords = bulk.coords + np.array((
    0,
    0,
    3.7))


sim = Combine(bulk,flake1)
sim = Combine(sim,flake2)
sim.box_dimensions[2] = 2*3.5+18
output = Writer(sim,'flake on graphite')
output.write_xyz('graphene.xyz')
output.write_lammps('data.flake')
