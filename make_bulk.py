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
bulk = Lattice(graphite.cell_shape())

cell_coords = graphite.cell_coords()
lattice_dimensions = bulk.lattice_size_vdw(vdw_cutoff)
lattice_points = bulk.create_lattice_points(lattice_dimensions)
coords = bulk.cell_onto_lattice(cell_coords,lattice_points)
molecule_labels = graphite.assign_molecules(lattice_dimensions)

print 'Cell dimensions:'
print graphite.cell_shape()
print 'Lattice structure:'
print bulk.lattice_size_vdw(vdw_cutoff)
print 'System size:'
print bulk.system_size(lattice_dimensions)

output = Writer(coords,molecule_labels)
output.write_xyz('bulk.xyz','bulk graphite')
output.write_lammps(bulk.system_size(lattice_dimensions),'data.bulk')
