import yaml
import makegraphitics as mg

config = yaml.load(open("config.yaml"))
forcefield = "OPLS"

# makes a graphene flake that is 50x50 Angstroms
motif = mg.molecules.Rectangle_Graphene(config, forcefield, 50, 50)
flake = mg.Crystal(motif, config, forcefield, [1, 1, 1])
vdw_defs = {1: 90, 2: 91}

mg.Parameterise(flake, vdw_defs)

name = "graphene"
output = mg.Writer(flake, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
