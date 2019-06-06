import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

vdw_defs = {1:90, 2:91}
R=15
GO_separation = 10 # approx 1 nm in experiment (with water!)
for i in range(5):
    motif = Hexagon_Graphene(config,forcefield,R)
    new_layer = Crystal(motif,config,forcefield,[1,1,1])
    Oxidiser(new_layer, ratio=2.5, video=20, new_island_freq=1e16, method='rf')
    Parameterise(new_layer,new_layer.vdw_defs)

    new_layer.coords = ( new_layer.coords 
                        + np.array((0, 0, i*GO_separation)))
    if i == 0:
        sim = new_layer
    else:
        sim = Combine(sim,new_layer)

name = 'GO_stack'
output = Writer(sim,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
