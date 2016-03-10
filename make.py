import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

#graphite = Graphite(config,forcefield)
#bulk = Crystal(graphite,config,forcefield,'vdw')

molecule = Hexagon_Graphene(config,forcefield,10)
sim = Crystal(molecule,config,forcefield,[1,1,1])

#output = Shifter(sim)
#output.rotate(60,20)
output2 = Writer(sim,'factory class test')
output2.write_xyz('graphene.xyz')

#output.write_lammps('data.graphene')
