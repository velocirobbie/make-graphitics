import makegraphitics as mg

R = 40
motif = mg.molecules.Hexagon_Graphene(R)
flake = mg.Crystal(motif, [1, 1, 1])
vdw_defs = {1: 90, 2: 91}

mg.Parameterise(flake, vdw_defs)

name = "graphene"
output = mg.Writer(flake, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
