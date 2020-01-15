import yaml
import numpy as np
import makegraphitics as mg

config = yaml.load(open("config.yaml"))
forcefield = "OPLS"

vdw_defs = {1: 90, 2: 91}

graphite = mg.molecules.Graphene(config, forcefield)
sim = mg.Crystal(graphite, config, forcefield, [40, 30, 1])
sim.vdw_defs = vdw_defs
mg.Parameterise(sim, vdw_defs)

j = 0
for i in [
    3,
    3,
    4,
    4,
    5,
    5,
    0,
    0,
    0,
    15,
    13,
    11,
    9,
    7,
    12,
    10,
    8,
    6,
    9,
    7,
    6,
    5,
    6,
    5,
    4,
    3,
    4,
    2,
    1,
    0,
    -1,
]:
    print j
    j += 1
    motif = mg.molecules.Hexagon_Graphene(config, forcefield, 5 + i * 2.44)
    next_layer = mg.Crystal(motif, config, forcefield, [1, 1, 1])

    next_layer.coords = next_layer.coords + np.array(
        (
            10 * 2 * (3 ** 0.5) * config[forcefield]["CC"],
            30 * config[forcefield]["CC"],
            j * config[forcefield]["layer_gap"],
        )
    )
    mg.Parameterise(next_layer, vdw_defs)
    sim = mg.Combine(sim, next_layer)
"""
for vector in [[4,10,10],[6,65,22],[15,40,10],[11,14,25],[10,70,10],[19,67,13]]:
    motif = Hexagon_Graphene(config,forcefield,8)
    new_flake = Crystal(motif,config,forcefield,[1,1,1])
    new_flake.coords = new_flake.coords + np.array((
        vector[0] * 2*(3**0.5) * config[forcefield]['CC'],
        vector[1] * config[forcefield]['CC'],
        vector[2] * config[forcefield]['layer_gap']))
    Parameterise(new_flake,vdw_defs)
    sim = Combine(sim,new_flake)
"""


name = "graphene"
output = mg.Writer(sim, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
