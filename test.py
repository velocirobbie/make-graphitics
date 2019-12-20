import yaml
from scripts.molecules import *
from scripts import *
import numpy as np
import json

config = yaml.load(open("config.yaml"))
forcefield = "OPLS"

# simple xyz hexagon
flake_radius = 10
motif = Hexagon_Graphene(config, forcefield, flake_radius)
flake = Crystal(motif, config, forcefield, [1, 1, 1])
name = "test"
output = Writer(flake, name)
output.write_xyz(name)

# GO
flake_radius = 10
motif = Hexagon_Graphene(config, forcefield, flake_radius)
flake = Crystal(motif, config, forcefield, [1, 1, 1])
oxidiser = Oxidiser(ratio=2.5, new_island_freq=1e15, method="rf")
go_flake = oxidiser.react(flake)
Parameterise(go_flake)
name = "test"
output = Writer(go_flake, name)
output.write_xyz(name)
output.write_lammps(name)
output.write_reaxff(name)

# rectangle
motif = Rectangle_Graphene(config, forcefield, 10, 10)
rect_flake = Crystal(motif, config, forcefield, [1, 1, 1])
rect_flake.vdw_defs = {1: 90, 2: 91}
Parameterise(rect_flake)

# combine
sim = Combine(go_flake, rect_flake)
print sim.vdw_defs
Parameterise(sim)
output = Writer(sim, name)
output.write_xyz(name)

# sheet
motif = Graphene(config, forcefield)
sheet = Crystal(motif, config, forcefield, [4, 4, 1])
oxidiser = Oxidiser(ratio=1, new_island_freq=1e15, method="rf")
go_sheet = oxidiser.react(sheet)
Parameterise(go_sheet, sheet.vdw_defs)
output = Writer(go_sheet, name)
output.write_lammps(name)

# graphite
motif = Graphite(config, forcefield)
graphite = Crystal(motif, config, forcefield, [4, 4, 2])
vdw_defs = {1: 90}
Parameterise(graphite, vdw_defs)
output.write_reaxff(name)
