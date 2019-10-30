from lattice import Lattice
from connector import Connector
import numpy as np

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
        self.bond_types, self.angle_types, self.dihedral_types, self.improper_types = self.molecule.connection_types()
        self.bonds = self.molecule.assign_bonds(
                self.lattice_dimensions)
        self.bond_labels = connect.bond_labels(
                self.atom_labels,self.bonds,self.bond_types)

        self.bond_graph = self.generate_bond_graph(self.bonds)

        self.angles = connect.angles(self.bonds, self.bond_graph)
        self.angle_labels = connect.angle_labels(
                self.atom_labels,self.angles,self.angle_types)

        self.dihedrals = connect.dihedrals(self.bonds, self.bond_graph)
        self.dihedral_labels = connect.dihedral_labels(
                self.atom_labels,self.dihedrals,self.dihedral_types)

        self.impropers = connect.impropers(self.bonds, self.bond_graph)
        self.improper_labels = connect.improper_labels(
                self.atom_labels,self.impropers,self.improper_types)
        
    def generate_bond_graph(self, bonds):
        N = int(np.amax(bonds)) # Number of atoms
        bond_graph = dict()
        for i in xrange(N):
            bond_graph[i] = set()

        for bond in bonds:
            bond_graph[bond[0]-1].add(bond[1]-1)
            bond_graph[bond[1]-1].add(bond[0]-1)
        return bond_graph
            

