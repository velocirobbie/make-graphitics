from scripts.lattice import Lattice
from scripts.coronene_cell import Coronene
from scripts.connector  import Connector
import yaml

def process(params):
    config = yaml.load(open('config.yaml'))
    CC = config[params]['CC']
    layer_gap = config[params]['layer_gap']
    vdw_cutoff = config['system']['vdw_cutoff']
    N_layers = config['system']['N_layers']
    CH = config[params]['CH']
    dq = config[params]['dq']

    coronene = Coronene(CC,CH,layer_gap)
    cell_coords = coronene.cell_coords()
    stack = Lattice(coronene.cell_shape())
    connect = Connector()

    lattice_dimensions = [1,1,2]
    lattice_points = stack.create_lattice_points(lattice_dimensions)

    coords = stack.cell_onto_lattice(cell_coords,lattice_points)
    molecule_labels = coronene.assign_molecules(lattice_dimensions)
    atom_labels = coronene.assign_atom_labels(lattice_dimensions)
    atom_charges = coronene.assign_atom_charges(lattice_dimensions,dq)
    bond_types =    [[1,1],
                     [1,2]]
    angle_types =   [[1,1,1],
                     [1,1,2]]
    torsion_types = [[1,1,1,1],
                     [1,1,1,2],
                     [2,1,1,2]]
    bonds = connect.coronene_bonds(lattice_dimensions)
    bond_labels =connect.bond_labels(atom_labels,bonds,bond_types)
    angles = connect.angles(bonds)
    angle_labels =connect.angle_labels(atom_labels,angles,angle_types) 
    torsions = connect.torsions(bonds)
    torsion_labels =connect.torsion_labels(atom_labels,torsions,torsion_types)
    box_dimensions = stack.system_size(lattice_dimensions)
    
    return [coords, atom_labels, atom_charges, molecule_labels, 
            bonds, bond_labels,
            angles, angle_labels,
            torsions, torsion_labels,
            box_dimensions]

