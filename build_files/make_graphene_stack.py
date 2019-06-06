import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

vdw_defs = {1:90, 2:91}
R=25
for i in range(5):
    motif = Hexagon_Graphene(config,forcefield,R)
    new_layer = Crystal(motif,config,forcefield,[1,1,1])

    new_layer.coords = ( new_layer.coords 
                        + np.array((0, 0, i*config[forcefield]['layer_gap'])))
    if i == 0:
        sim = new_layer
    else:
        sim = Combine(sim,new_layer)

Parameterise(sim,vdw_defs)

name = 'graphene_stack'
output = Writer(sim,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
