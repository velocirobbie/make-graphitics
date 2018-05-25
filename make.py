import yaml
from scripts.molecules import *
from scripts import *
import numpy as np

config = yaml.load(open('config.yaml'))
forcefield = 'OPLS'

motif = Graphene(config,forcefield)
layer = Crystal(motif,config,forcefield,[82,48,1])
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
