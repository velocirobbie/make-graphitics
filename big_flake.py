import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'GraFF_5'

#config[forcefield]['CC'] = config[forcefield]['CC'] * 0.999
for i in [2]:
    graphite = Graphite(config,forcefield)
    molecule = Hexagon_Graphene(config,forcefield,100)
    flake = Crystal(molecule,config,forcefield,[1,1,1])
#make flake carbons different to bulk
#for atom in range(molecule.natoms):
#    if flake.atom_labels[atom] == 1:
#        flake.atom_labels[atom] = 3


    flake.coords = flake.coords + np.array((
        30 * 2*(3**0.5) * config[forcefield]['CC'],
        108*config[forcefield]['CC'],
        (i*2)*config[forcefield]['layer_gap']+3.7))


#output = Shifter(sim,'lammps')
#output.rotate(180,1)
    output = Writer(flake,'flake on graphite')
    output.write_xyz('graphene_bigflake.xyz')
    output.write_lammps('data.bigflake')
