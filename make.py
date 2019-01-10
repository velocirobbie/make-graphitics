import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

##R=25
#motif = Hexagon_Graphene(config,forcefield,R)
#layer = Crystal(motif,config,forcefield,[1,1,1])
#vdw_defs = {1:90, 2:91}


motif = Graphite(config,forcefield)
layer = Crystal(motif,config,forcefield,[20,10,1])
vdw_defs = {1:90}

Parameterise(layer,vdw_defs)

"""
qs = ['atom','bond','angle','torsion','improper']
off =[28, 22,47, 75, 5]
for i in range(len(qs)):
    qlist = qs[i]+'_labels'
    q = np.array(getattr(layer, qlist))
    print type(q)
    setattr(layer, qlist, q + off[i]) 
"""
name = 'graphene'
output = Writer(layer,name)
output.write_xyz(name+'.xyz')
output.write_lammps(name+'.data')
