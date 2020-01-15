import yaml
import numpy as np
import makegraphitics as mg

config = yaml.load(open("config.yaml"))
forcefield = "OPLS"

vdw_defs = {1: 90, 2: 91}
R = 15
GO_separation = 10  # approx 1 nm in experiment (with water!)

oxidiser = mg.reactors.Oxidiser(ratio=2.5, video_xyz=20,
                                new_island_freq=1e14, method="rf")

for i in range(3):
    motif = mg.molecules.Hexagon_Graphene(config, forcefield, R)
    new_layer = mg.Crystal(motif, config, forcefield, [1, 1, 1])
    new_layer = oxidiser.react(new_layer)

    new_layer.coords = new_layer.coords + np.array((0, 0, i * GO_separation))
    if i == 0:
        sim = new_layer
    else:
        sim = mg.Combine(sim, new_layer)

mg.Parameterise(sim, new_layer.vdw_defs)

name = "GO_stack"
output = mg.Writer(sim, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
