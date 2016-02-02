from scripts.lattice import Lattice
from scripts.graphene_cell import Graphene
from scripts.write_coords import Writer
import yaml

config = yaml.load(open('config.yaml'))
CC = config['crystal']['CC']
layer_gap = config['crystal']['layer_gap']
vdw_cutoff = config['system']['vdw_cutoff']
N_layers = 1

graphene = Graphene(CC)
sheet = Lattice(graphene.orthorhombic_graphene_shape())

cell_coords = graphene.orthorhombic_graphene_coords()
lattice_dimensions = sheet.lattice_size_layers(vdw_cutoff,N_layers)
lattice_points = sheet.create_lattice_points(lattice_dimensions)
coords = sheet.cell_onto_lattice(cell_coords,lattice_points)

print 'Cell dimensions:'
print graphene.orthorhombic_graphene_shape()
print 'Lattice structure:'
print sheet.lattice_size_layers(vdw_cutoff,N_layers)
print 'System size:'
print sheet.system_size(lattice_dimensions)

output = Writer(coords)
output.write_xyz('graphene.xyz','graphene')
output.write_lammps(sheet.system_size(lattice_dimensions),'data.graphene')
