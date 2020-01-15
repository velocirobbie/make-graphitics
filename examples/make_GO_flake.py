import yaml
import makegraphitics as mg

config = yaml.load(open("config.yaml"))
forcefield = "OPLS"
flake_radius = 25
layout = [1, 1, 1]  # make a 1x1x1 array of flakes

motif = mg.molecules.Hexagon_Graphene(config, forcefield, flake_radius)
flake = mg.Crystal(motif, config, forcefield, layout)

oxidiser = mg.reactors.Oxidiser(ratio=2.5, video_xyz=20,
                                new_island_freq=1e14, method="rf")
flake = oxidiser.react(flake)

name = "GO_flake"
output = mg.Writer(flake, name)
output.write_xyz(name + ".xyz")

mg.Parameterise(flake)

output = mg.Writer(flake, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
