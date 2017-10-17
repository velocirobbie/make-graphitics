import yaml
import numpy as np
from scripts.molecules import *
from scripts import *

print '===='
config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'
R=20
motif = Graphene(config,forcefield)
layer = Crystal(motif,config,forcefield,[10,10,1])

a = Oxidiser(layer)

Parameterise(layer, a.vdw_defs)

layer.coords = layer.coords + np.array(([R+5,R+5,7]))

#motif2 = Graphene(config,forcefield)
#base = Crystal(motif2,config,forcefield,[70,35,1])

#b = Oxidiser(base)

#Parameterise(base, {1:90, 2:91})

name = 'test'
output = Writer(layer,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')

