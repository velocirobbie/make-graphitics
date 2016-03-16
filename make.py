import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

graphene = Graphene(config,forcefield)
print graphene.assign_bonds([2,2,1])
print graphene.assign_atom_labels([2,2,1])
bulk = Crystal(graphene,config,forcefield,'layers')

output = Writer(bulk,'graphene')
output.write_xyz('graphene.xyz')
output.write_lammps('data.bulk')
