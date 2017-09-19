import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

motif = Hexagon_Graphene(config,forcefield,100)
layer = Crystal(motif,config,forcefield,[1,1,1])

#motif = Graphene(config,forcefield)
#layer = Crystal(motif,config,forcefield,[30,15,1])

a = Oxidiser(layer)

name = 'graphene'
output = Writer(layer,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
