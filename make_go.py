import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

motif = Hexagon_Graphene(config,forcefield,50)
layer = Crystal(motif,config,forcefield,[1,1,1])

#motif = Graphene(config,forcefield)
#layer = Crystal(motif,config,forcefield,[30,15,1])

a = Oxidiser(layer)

layer.coords = layer.coords + np.array(([50,50,5]))

name = 'graphene'
output = Writer(layer,name)
output.write_xyz(name+'.xyz')

p = Parameterise('hello','scripts/params/oxidise_types.yaml')

output.bond_coeffs = p.match_bonds(a.crystal.bond_types)
output.angle_coeffs = p.match_angles(a.crystal.angle_types)
output.torsion_coeffs = p.match_torsions(a.crystal.torsion_types)
output.improper_coeffs = p.match_impropers(a.crystal.improper_types)
output.pair_coeffs = p.pair_coeffs()

output.write_lammps(name+'.data')
