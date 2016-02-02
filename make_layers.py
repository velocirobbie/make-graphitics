from scripts.lattice import Lattice
from scripts.graphite_cell import Graphite
from scripts.write_coords import Writer
import yaml

config = yaml.load(open('config.yaml'))
CC = config['crystal']['CC']
layer_gap = config['crystal']['layer_gap']
vdw_cutoff = config['system']['vdw_cutoff']
N_layers = config['system']['N_layers']

graphite = Graphite(CC,layer_gap)
layers = Lattice(graphite.orthorhombic_ABgraphite_shape())

cell_coords = graphite.orthorhombic_ABgraphite_coords()
lattice_dimensions = layers.lattice_size_layers(vdw_cutoff,N_layers)
lattice_points = layers.create_lattice_points(lattice_dimensions)
coords = layers.cell_onto_lattice(cell_coords,lattice_points)

print 'Cell dimensions:'
print graphite.orthorhombic_ABgraphite_shape()
print 'Lattice structure:'
print lattice_dimensions
print 'System size:'
print layers.system_size(lattice_dimensions)

output = Writer(coords)
output.write_xyz('layers.xyz','layers of graphite')
output.write_lammps(layers.system_size(lattice_dimensions),'data.layers')
