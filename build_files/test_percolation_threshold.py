import yaml
from scripts.molecules import *
from scripts import *
import numpy as np
import json
from math import pi, cos
import sys

np.random.seed(6)
config = yaml.load(open("config.yaml"))
forcefield = "OPLS"
x_length = 200
y_length = 200

# calculate array of unit cells to make sheet
# unit cell is the orthorombic unit cell of graphene
unit_cell_x = 2.0 * config[forcefield]["CC"] * cos(pi / 6.0)
unit_cell_y = 3.0 * config[forcefield]["CC"]
x_cells = int(x_length / unit_cell_x)
y_cells = int(y_length / unit_cell_y)
layout = [x_cells, y_cells, 1]  # make an array of unit cells with this dimension

for rate_exp in [16, 17, 18, 19, 20]:
    motif = Graphene(config, forcefield)
    sheet = Crystal(motif, config, forcefield, layout)

    Oxidiser(sheet, ratio=2.5, video=20, new_island_freq=10 ** rate_exp, method="rf")
    Parameterise(sheet, sheet.vdw_defs)

    name = "GO_sheet_" + str(rate_exp)
    output = Writer(sheet, name)
    output.write_xyz(name + ".xyz")
    output.write_lammps(name + ".data")
