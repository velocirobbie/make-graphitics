from lattice import Lattice
from connector import Connector

class Crystal(object):
    def __init__(self,molecule,config,
            forcefield=[],lattice_dimensions=[]):

        self.molecule = molecule
        self.config = config
        if forcefield: self.forcefield = forcefield
        else: self.forcefield = 'OPLS'
        self.lattice = self.init_lattice()
        self.lattice_dimensions = self.determine_lattice(
                lattice_dimensions)
        self.generate_structure()
        self.generate_connections()

    def init_lattice(self):
        self.cell_coords = self.molecule.cell_coords()
        return Lattice(self.molecule.cell_shape())

    def determine_lattice(self,lattice_dimensions):
        if not lattice_dimensions:
            return [1,1,1]
        elif lattice_dimensions == 'vdw':
            return self.lattice.lattice_size_vdw(
                    self.config['system']['vdw_cutoff'])
        elif lattice_dimensions == 'layers':
            return self.lattice.lattice_size_layers(
                    self.config['system']['vdw_cutoff'],
                    self.config['system']['N_layers'])
        else: return lattice_dimensions
    
    def generate_structure(self):
        self.lattice_points = self.lattice.create_lattice_points(
                self.lattice_dimensions)
        self.coords = self.lattice.cell_onto_lattice(
                self.cell_coords,self.lattice_points)
        self.box_dimensions = self.lattice.system_size(
                self.lattice_dimensions)
        self.molecule_labels = self.molecule.assign_molecules(
                self.lattice_dimensions)
        self.atom_labels = self.molecule.assign_atom_labels(
                self.lattice_dimensions)
        self.atom_charges = self.molecule.assign_atom_charges(
                self.lattice_dimensions,self.config[self.forcefield]['dq'])
    def generate_connections(self):
        connect = Connector()
        bond_types, angle_types, torsion_types = self.molecule.connection_types()
        self.bonds = self.molecule.assign_bonds(
                self.lattice_dimensions)
        self.bond_labels = connect.bond_labels(
                self.atom_labels,self.bonds,bond_types)

        self.angles = connect.angles(self.bonds)
        self.angle_labels = connect.angle_labels(
                self.atom_labels,self.angles,angle_types)

        self.torsions = connect.torsions(self.bonds)
        self.torsion_labels = connect.torsion_labels(
                self.atom_labels,self.torsions,torsion_types)

                
