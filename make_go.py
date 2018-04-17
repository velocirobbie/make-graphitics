import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'
#R=30
#motif = Hexagon_Graphene(config,forcefield,R)
#layer = Crystal(motif,config,forcefield,[1,1,1])

motif = Graphene(config,forcefield)
layer = Crystal(motif,config,forcefield,[30,17,1])

np.random.seed(42)
a = Oxidiser(layer, ratio=1, video=100)
"""
p = Parameterise(layer,a.vdw_defs)

layer.bond_coeffs = p.match_bonds(layer.bond_types)
layer.angle_coeffs = p.match_angles(layer.angle_types)
layer.torsion_coeffs = p.match_torsions(layer.torsion_types)
layer.improper_coeffs = p.match_impropers(layer.improper_types)
layer.pair_coeffs = p.match_pairs()
layer.masses = p.match_masses()
print layer.angle_labels
print layer.masses
print layer.pair_coeffs
print layer.bond_coeffs
"""
#layer.coords = layer.coords + np.array(([R,R,5]))

name = 'graphene'
output = Writer(layer,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
