from lattice import Lattice
from connector import Connector
import numpy as np
import os
import yaml


class Sim(object):
    def generate_connections(self):
        try:
            self.atom_labels, self.bonds
        except AttributeError:
            raise Exception("Simulation has not been assigned atom types and bonds yet")

        connect = Connector()
        self.bond_types = connect.find_bond_types(self.atom_labels, self.bonds)
        self.bond_labels = connect.bond_labels(
            self.atom_labels, self.bonds, self.bond_types
        )

        self.bond_graph = self.generate_bond_graph(self.bonds)

        self.angles = connect.angles(self.bonds, self.bond_graph)
        self.angle_types = connect.find_angle_types(self.atom_labels, self.angles)
        self.angle_labels = connect.angle_labels(
            self.atom_labels, self.angles, self.angle_types
        )

        self.dihedrals = connect.dihedrals(self.bonds, self.bond_graph)
        self.dihedral_types = connect.find_dihedral_types(
            self.atom_labels, self.dihedrals
        )
        self.dihedral_labels = connect.dihedral_labels(
            self.atom_labels, self.dihedrals, self.dihedral_types
        )

        self.impropers = connect.impropers(self.bonds, self.bond_graph)
        self.improper_types = connect.find_improper_types(
            self.atom_labels, self.impropers
        )
        self.improper_labels = connect.improper_labels(
            self.atom_labels, self.impropers, self.improper_types
        )

    def generate_bond_graph(self, bonds):
        N = int(np.amax(bonds))  # Number of atoms
        bond_graph = dict()
        for i in xrange(N):
            bond_graph[i] = set()

        for bond in bonds:
            bond_graph[bond[0] - 1].add(bond[1] - 1)
            bond_graph[bond[1] - 1].add(bond[0] - 1)
        return bond_graph

    def bonded_to(self, centre):
        return list(self.bond_graph[centre])

    def crystal_params(self):
        path = os.path.dirname(__file__) + "/params/"
        return yaml.load(open(path + "config.yaml"), Loader=yaml.FullLoader)
