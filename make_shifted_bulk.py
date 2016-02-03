from scripts.lattice import Lattice
from scripts.graphite_cell import Graphite
from scripts.write_coords import Writer
from scripts.shifty import Shifter
import yaml

config = yaml.load(open('config.yaml'))
CC = config['crystal']['CC']
layer_gap = config['crystal']['layer_gap']
vdw_cutoff = config['system']['vdw_cutoff']
N_layers = config['system']['N_layers']

graphite = Graphite(CC,layer_gap)
bulk = Lattice(graphite.cell_shape())

cell_coords = graphite.cell_coords()
lattice_dimensions = bulk.lattice_size_vdw(vdw_cutoff)
lattice_points = bulk.create_lattice_points(lattice_dimensions)
coords = bulk.cell_onto_lattice(cell_coords,lattice_points)
molecule_labels = graphite.assign_molecules(lattice_dimensions)

system_size = bulk.system_size(lattice_dimensions)

print 'Cell dimensions:'
print graphite.cell_shape()
print 'Lattice structure:'
print lattice_dimensions
print 'System size:'
print system_size

system_size[2] = 70

output = Shifter(coords,molecule_labels,8,'lammps',system_size)
output.z_shift(12,13,1)
