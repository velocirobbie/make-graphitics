import yaml
from scripts.molecules import *
from scripts import *
import numpy as np
import json

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'
R=30
motif = Hexagon_Graphene(config,forcefield,R)
layer = Crystal(motif,config,forcefield,[1,1,1])

#motif = Graphene(config,forcefield)
#layer = Crystal(motif,config,forcefield,[30,17,1])
#np.random.seed(12)
a = Oxidiser(layer, ratio=2.5, video=20, new_island_freq=1e16)#, method='rf')
json.dump(a.vdw_defs,open('vdw_defs.json','w'))
p = Parameterise(layer,a.vdw_defs)
#r = Reducer(layer,a.vdw_defs)
#offset for polymerisation step
#qs = ['atom','bond','angle','torsion','improper']
#off =[24, 24,47, 56, 4]
#for i,q_i in enumerate(qs):
#    qlist = q_i+'_labels'
#    q = np.array(getattr(layer, qlist))
#    print type(q)
#    setattr(layer, qlist, q + off[i]) 

#layer.coords = layer.coords + np.array(([R,R,5]))

name = 'graphene'
output = Writer(layer,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
