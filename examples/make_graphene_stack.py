import yaml
import numpy as np
import makegraphitics as mg

config = yaml.load(open("config.yaml"), Loader=yaml.FullLoader)
forcefield = "OPLS"

vdw_defs = {1: 90, 2: 91}
R = 25
for i in range(5):
    motif = mg.molecules.Hexagon_Graphene(R)
    new_layer = mg.Crystal(motif, [1, 1, 1])

    new_layer.coords = new_layer.coords + np.array(
        (0, 0, i * config[forcefield]["layer_gap"])
    )
    new_layer.vdw_defs = {1: 90, 2: 91}
    if i == 0:
        sim = new_layer
    else:
        sim = mg.Combine(sim, new_layer)

mg.Parameterise(sim)

name = "graphene_stack"
output = mg.Writer(sim, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
