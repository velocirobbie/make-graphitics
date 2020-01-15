import yaml
import numpy as np
import makegraphitics as mg

forcefield = "GraFF_5"

config = yaml.load(open("config.yaml"))

graphite = mg.molecules.Graphite(config, forcefield)
bulk = mg.Crystal(graphite, config, forcefield, [122, 71, 2])


molecule1 = mg.molecules.Hexagon_Graphene(config, forcefield, 50)
flake1 = mg.Crystal(molecule1, config, forcefield, [1, 1, 1])
# make flake carbons different to bulk
# for atom in range(molecule.natoms):
#    if flake.atom_labels[atom] == 1:
#        flake.atom_labels[atom] = 3

flake1.coords = flake1.coords + np.array(
    (
        20 * 2 * (3 ** 0.5) * config[forcefield]["CC"],
        72 * config[forcefield]["CC"],
        3.7 + (4) * config[forcefield]["layer_gap"],
    )
)
bulk.coords = bulk.coords + np.array((0, 0, 3.7))


sim = mg.Combine(bulk, flake1)
sim.box_dimensions[2] = 30
# output = Shifter(sim,'lammps')
# output.rotate(180,1)
output = mg.Writer(sim, "flake on graphite")
output.write_xyz("graphene.xyz")
output.write_lammps("data.flake")
