import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open("config.yaml"))
forcefield = "OPLS"
R = 10
motif = Hexagon_Graphene(config, forcefield, R)
layer = Crystal(motif, config, forcefield, [1, 1, 1])

# motif = Graphene(config,forcefield)
# layer = Crystal(motif,config,forcefield,[30,15,1])

a = Oxidiser(layer)

Parameterise(layer, a.vdw_defs)

# print layer.angle_labels
# print layer.masses
# print layer.pair_coeffs
# print layer.bond_coeffs

layer.coords = layer.coords + np.array(([R + 5, R + 5, 7]))

motif2 = Graphene(config, forcefield)
base = Crystal(motif2, config, forcefield, [45, 25, 1])
print base.box_dimensions
# base.coords = base.coords + np.array(([R,R,0]))

sim = Combine(layer, base)

name = "graphene"
output = Writer(sim, name)
output.write_xyz(name + ".xyz")


output.write_lammps(name + ".data")
