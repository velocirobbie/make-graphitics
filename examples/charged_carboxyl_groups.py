import makegraphitics as mg

flake_radius = 25
layout = [1, 1, 1]  # make a 1x1x1 array of flakes

motif = mg.molecules.Hexagon_Graphene(flake_radius)
flake = mg.Crystal(motif, layout)

oxidiser = mg.reactors.Oxidiser(
    ratio=2.5, new_island_freq=1e14, method="rf",
    carboxyl_charged_ratio=0.5
)
flake = oxidiser.react(flake)

mg.Parameterise(flake)

flake.validate()
