import yaml
import numpy as np
from scripts.molecules import *
from scripts import *

print '===='
config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'
R=50
motif = Graphene(config,forcefield)
layer = Crystal(motif,config,forcefield,[100,58,1])
#motif = Hexagon_Graphene(config,forcefield,R)
#layer = Crystal(motif,config,forcefield,[1,1,1])

#a = Oxidiser(layer)

name = 'test'
output = Writer(layer,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
#Parameterise(layer, a.vdw_defs)

#layer.coords = layer.coords + np.array(([R+5,R+5,7]))

#motif2 = Graphene(config,forcefield)
#base = Crystal(motif2,config,forcefield,[70,35,1])

#b = Oxidiser(base)

#Parameterise(base, {1:90, 2:91})
#output = Writer(layer,name)
#output.write_xyz(name+'.xyz')

#output.write_lammps(name+'.data')

