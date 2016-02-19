from scripts.lattice import Lattice
from scripts.graphite_cell import Graphite
from scripts.connector  import Connector
import yaml

def process(params,layer_or_bulk):
    config = yaml.load(open('config.yaml'))
    CC = config[params]['CC']
    layer_gap = config[params]['layer_gap']
    vdw_cutoff = config['system']['vdw_cutoff']
    N_layers = config['system']['N_layers']

    graphite = Graphite(CC,layer_gap)
    cell_coords = graphite.cell_coords()
    bulk = Lattice(graphite.cell_shape())
    connect = Connector()

    if layer_or_bulk == 'bulk':
        lattice_dimensions = bulk.lattice_size_vdw(vdw_cutoff)
    elif layer_or_bulk == 'layer':
        lattice_dimensions = bulk.lattice_size_layers(
                vdw_cutoff,N_layers)
    else:
        raise ValueError(layers_or_bulk,'is not a valid choice')

    lattice_points = bulk.create_lattice_points(lattice_dimensions)

    coords = bulk.cell_onto_lattice(cell_coords,lattice_points)
    molecule_labels = graphite.assign_molecules(lattice_dimensions)
    bonds = connect.graphite_bonds(lattice_dimensions)
    angles = connect.angles(bonds)
    torsions = connect.torsions(bonds)
    box_dimensions = bulk.system_size(lattice_dimensions)
    
    return [coords, molecule_labels, bonds, angles, torsions, box_dimensions]

