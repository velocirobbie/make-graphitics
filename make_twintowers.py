import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'GraFF_5'
i=1
#config[forcefield]['CC'] = config[forcefield]['CC'] * 0.999
graphite = Graphene(config,forcefield)
bulk = Crystal(graphite,config,forcefield,[50,20,i])

molecule1 = Hexagon_Graphene(config,forcefield,15)
flake1 = Crystal(molecule1,config,forcefield,[1,1,15])
molecule2 = Hexagon_Graphene(config,forcefield,15)
flake2 = Crystal(molecule2,config,forcefield,[1,1,15])
molecule3 = Hexagon_Graphene(config,forcefield,9)
plane1 = Crystal(molecule3,config,forcefield,[1,1,1])




flake1.coords = flake1.coords + np.array((
    5 * 2*(3**0.5) * config[forcefield]['CC'],
    24*config[forcefield]['CC'],
    (i)*config[forcefield]['layer_gap']+3.7))
flake2.coords = flake2.coords + np.array((
    15* 2*(3**0.5) * config[forcefield]['CC'],
    24*config[forcefield]['CC'],
    (i)*config[forcefield]['layer_gap']+3.7))
bulk.coords = bulk.coords + np.array((
    0,
    0,
    3.7))
plane1.coords = plane1.coords + np.array((
    21* 2*(3**0.5) * config[forcefield]['CC'],
    40*config[forcefield]['CC'],
    10*config[forcefield]['layer_gap']+3.7))


sim = Combine(bulk,flake1)
sim = Combine(sim,flake2)
sim = Combine(sim,plane1)
sim.box_dimensions[2] = i*2*3.5+18
#output = Shifter(sim,'lammps')
#output.rotate(180,1)
output = Writer(sim,'flake on graphite')
output.write_xyz('graphene'+str(i)+'.xyz')
output.write_lammps('data.flake'+str(i))
