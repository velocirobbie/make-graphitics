import makegraphitics as mg

flake_radius = 25
layout = [1, 1, 1]  # make a 1x1x1 array of flakes

motif = mg.molecules.Hexagon_Graphene(flake_radius)
flake = mg.Crystal(motif, layout)

oxidiser = mg.reactors.Oxidiser(
    ratio=2.5, video_xyz=20, new_island_freq=1e14, method="rf"
)
flake = oxidiser.react(flake)

mg.Parameterise(flake)

name = "graphene"
output = mg.Writer(flake, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
