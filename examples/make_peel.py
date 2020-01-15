import yaml
import numpy as np
import math
import makegraphitics as mg

config = yaml.load(open("config.yaml"))
forcefield = "GraFF_5"
graphite = mg.molecules.Graphite(config, forcefield)
bulk = mg.Crystal(graphite, config, forcefield, [21, 12, 1])
bulk.coords = bulk.coords + np.array((0, 0, 1 * config[forcefield]["layer_gap"]))


molecule1 = mg.molecules.Hexagon_Graphene(config, forcefield, 15)
flake1 = mg.Crystal(molecule1, config, forcefield, [1, 1, 1])

flake1.coords = flake1.coords + np.array(
    (
        10 * 2 * math.cos(math.pi / 6) * config[forcefield]["CC"],
        6 * 3 * config[forcefield]["CC"],
        3 * config[forcefield]["layer_gap"],
    )
)

bulk.vdw_defs = {1: 90}
flake1.vdw_defs = {1: 90, 2: 91}
sim = mg.Combine(bulk, flake1)

output = mg.Writer(sim, "flake on graphite")
output.write_xyz("graphene" + str(1) + ".xyz")
output.write_lammps("flake" + str(1) + ".data")
