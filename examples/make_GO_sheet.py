import yaml
from math import pi, cos
import makegraphitics as mg

config = yaml.load(open("config.yaml"), Loader=yaml.FullLoader)
forcefield = "OPLS"
x_length = 20
y_length = 20

# calculate array of unit cells to make sheet
# unit cell is the orthorombic unit cell of graphene
unit_cell_x = 2.0 * config[forcefield]["CC"] * cos(pi / 6.0)
unit_cell_y = 3.0 * config[forcefield]["CC"]
x_cells = int(x_length / unit_cell_x)
y_cells = int(y_length / unit_cell_y)
layout = [x_cells, y_cells, 1]  # make an array of unit cells with this dimension

motif = mg.molecules.Graphene(forcefield=forcefield)
sheet = mg.Crystal(motif, layout)

oxidiser = mg.reactors.Oxidiser(ratio=2.5, video_xyz=20,
                                new_island_freq=1e14, method="rf")
sheet = oxidiser.react(sheet)

mg.Parameterise(sheet, sheet.vdw_defs)

name = "GO_sheet"
output = mg.Writer(sheet, name)
output.write_xyz(name + ".xyz")
output.write_lammps(name + ".data")
