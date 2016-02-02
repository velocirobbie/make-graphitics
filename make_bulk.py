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
bulk = Lattice(graphite.orthorhombic_ABgraphite_shape())

cell_coords = graphite.orthorhombic_ABgraphite_coords()
lattice_dimensions = bulk.lattice_size_vdw(vdw_cutoff)
lattice_points = bulk.create_lattice_points(lattice_dimensions)
coords = bulk.cell_onto_lattice(cell_coords,lattice_points)

print 'Cell dimensions:'
print graphite.orthorhombic_ABgraphite_shape()
print 'Lattice structure:'
print bulk.lattice_size_vdw(vdw_cutoff)
print 'System size:'
print bulk.system_size(lattice_dimensions)

output = Writer(coords)
output.write_xyz('bulk.xyz','bulk graphite')
output.write_lammps(bulk.system_size(lattice_dimensions),'in.bulk')
