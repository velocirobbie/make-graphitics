import numpy as np
from math import pi, sin, cos
import yaml
from molecule import Molecule


class GraphiteStrip(Molecule):
    # Handling AB graphite with orthorhombic unit cell (8 atom)
    # zig zag edges
    def __init__(self, config, forcefield, length):
        self.CC = config[forcefield]["CC"]
        self.CH = config[forcefield]["CH"]
        self.layer_gap = config[forcefield]["layer_gap"]
        self.length = length

    # orthorhombic unitcell lattice parameters
    def cell_shape(self):
        a = 2.0 * self.CC * cos(pi / 6.0)
        c = 2.0 * self.layer_gap

        b = 2.0 * self.CC + 2.0 * self.CH
        units_per_cell = 1
        while b < self.length:
            b += 3.0 * self.CC
            units_per_cell += 1

        if units_per_cell < 3:
            raise ValueError("make it bigger")
        self.units_per_cell = units_per_cell

        cell_dimensions = [a, b, c]
        return cell_dimensions

    def cell_coords(self):
        CC = self.CC
        CH = self.CH
        layer_gap = self.layer_gap
        cos_CC = cos(pi / 6.0) * CC
        sin_CC = 0.5 * CC

        def simple_unit():
            C1 = [0, 0, 0]
            C2 = [0, CC, 0]
            C3 = [cos_CC, CC + sin_CC, 0]
            C4 = [cos_CC, 2 * CC + sin_CC, 0]
            C5 = [0, CC, layer_gap]
            C6 = [0, CC * 2, layer_gap]
            C7 = [cos_CC, 2 * CC + sin_CC, layer_gap]
            C8 = [cos_CC, sin_CC, layer_gap]
            unit_coords = np.array([C1, C2, C3, C4, C5, C6, C7, C8])
            return unit_coords

        def bottom_unit():
            H1 = [CC - CH, 0, 0]
            C2 = [0, CC, 0]
            C3 = [cos_CC, CC + sin_CC, 0]
            C4 = [cos_CC, 2 * CC + sin_CC, 0]
            H5 = [0, 2 * CC - CH, layer_gap]
            C6 = [0, CC * 2, layer_gap]
            C7 = [cos_CC, 2 * CC + sin_CC, layer_gap]
            unit_coords = np.array([H1, C2, C3, C4, H5, C6, C7])
            return unit_coords

        def top_unit():
            C1 = [0, 0, 0]
            C2 = [0, CC, 0]
            C3 = [cos_CC, CC + sin_CC, 0]
            H4 = [cos_CC, CH + CC + sin_CC, 0]
            C5 = [0, CC, layer_gap]
            H6 = [0, CC + CH, layer_gap]
            C8 = [cos_CC, sin_CC, layer_gap]
            unit_coords = np.array([C1, C2, C3, H4, C5, H6, C7])
            return unit_coords

        cell_coords = bottom_unit()
        for i in range(self.units_per_cell - 2):
            unit = simple_unit()
            for j in unit:
                j[2] += (i + 1) * 3 * CC
            cell_coords = np.vstack(cell_cords, unit)
        cell_coords = np.vstack(cell_coords, top_unit())
        return cell_coords

    def assign_molecules(self, lattice_dimensions):
        l = [1, 1, 1, 1, 2, 2, 2]  # bottom
        for i in range(self.units_per_cell - 2):
            l += [1, 1, 1, 1, 2, 2, 2, 2]  # simple
        l += [1, 1, 1, 1, 2, 2, 2]  # top
        unit_cell_molecule_label = l

        molecule_labels = []
        for x in range(lattice_dimensions[0]):
            for y in range(lattice_dimensions[1]):
                for z in range(lattice_dimensions[2]):
                    labels = np.array(unit_cell_molecule_label) + (z * 2)
                    molecule_labels.extend(list(labels))
        return molecule_labels

    def assign_atom_labels(self, lattice_dimensions):
        l = [2, 1, 1, 1, 2, 1, 1]  # bottom
        for i in range(self.units_per_cell - 2):
            l += [1, 1, 1, 1, 1, 1, 1, 1]  # simple
        l += [1, 1, 1, 2, 1, 1, 2]  # top
        cell_label = l

        atom_labels = []
        for x in range(lattice_dimensions[0]):
            for y in range(lattice_dimensions[1]):
                for z in range(lattice_dimensions[2]):
                    atom_labels.extend(list(cell_labels))
        return atom_labels

    def assign_atom_charges(self, lattice_dimensions, q):
        l = [2, 1, 1, 1, 2, 1, 1]  # bottom
        for i in range(self.units_per_cell - 2):
            l += [1, 1, 1, 1, 1, 1, 1, 1]  # simple
        l += [1, 1, 1, 2, 1, 1, 2]  # top
        cell_label = l

        atom_charges = []
        cell_charges = [0, 0, 0, 0, 0, 0, 0, 0]
        for x in range(lattice_dimensions[0]):
            for y in range(lattice_dimensions[1]):
                for z in range(lattice_dimensions[2]):
                    atom_charges.extend(list(cell_charges))
        return atom_charges

    def assign_bonds(self, lattice_dimensions):
        internal_bonds = np.array([[1, 2], [2, 3], [3, 4], [5, 6], [6, 7], [8, 5]])
        bonds = np.empty((0, 2), dtype=int)
        # loop through all cells
        for x in range(lattice_dimensions[0]):
            for y in range(lattice_dimensions[1]):
                for z in range(lattice_dimensions[2]):
                    cell_position = [x, y, z]
                    (
                        xcell_position,
                        ycell_position,
                        xycell_position,
                    ) = find_adjacent_cells(cell_position, lattice_dimensions)

                    # Add bonds within cell
                    i = self.index_cell([x, y, z], lattice_dimensions, 8)
                    bonds = np.vstack((bonds, internal_bonds + i))
                    # Add bonds that cross cell boundries
                    bonds = self.add_cross_bond(
                        lattice_dimensions, xcell_position, [3, 2], bonds, i
                    )
                    bonds = self.add_cross_bond(
                        lattice_dimensions, xycell_position, [4, 1], bonds, i
                    )
                    bonds = self.add_cross_bond(
                        lattice_dimensions, ycell_position, [4, 1], bonds, i
                    )

                    bonds = self.add_cross_bond(
                        lattice_dimensions, xcell_position, [7, 6], bonds, i
                    )
                    bonds = self.add_cross_bond(
                        lattice_dimensions, xcell_position, [8, 5], bonds, i
                    )
                    bonds = self.add_cross_bond(
                        lattice_dimensions, ycell_position, [7, 8], bonds, i
                    )

        return bonds

    def index_cell(self, cell_position, lattice_dimensions, atoms_per_cell):
        N = atoms_per_cell
        total = cell_position[2] * N
        total += cell_position[1] * lattice_dimensions[2] * N
        total += cell_position[0] * lattice_dimensions[2] * lattice_dimensions[1] * N
        return total

    def add_cross_bond(self, lattice_dimensions, cell_position, atoms, bonds, i):
        # cell_position : the adjoining cell
        # atoms [i,j]: ith atom in cell, jth atom in adjoining
        atoms[0] += i
        atoms[1] += self.index_cell(cell_position, lattice_dimensions, 8)
        return np.vstack((bonds, atoms))

    def connection_types(self):
        bond_types = [[1, 1]]
        angle_types = [[1, 1, 1]]
        dihedral_types = [[1, 1, 1, 1]]
        improper_types = [[1, 1, 1, 1]]
        return bond_types, angle_types, dihedral_types, improper_types


def find_connections(bonds, centre):
    connections = np.where(bonds == centre)
    connections = np.vstack((connections[0], connections[1]))
    return connections.transpose()


def find_neighbours(bonds, centre):
    connections = find_connections(bonds, centre)
    neighbours = []
    for connection in connections:
        # Find atom connected to centre
        neighbour = bonds[connection[0]][connection[1] - 1]
        neighbours.append(neighbour)
    return neighbours


def find_adjacent_cells(cell_position, lattice_dimensions):
    x, y, z = cell_position
    # Account for periodic system
    if x == lattice_dimensions[0] - 1:
        xcell_position = [0, y, z]
    else:
        xcell_position = [x + 1, y, z]

    if y == lattice_dimensions[1] - 1:
        ycell_position = [x, 0, z]
    else:
        ycell_position = [x, y + 1, z]

    xycell_position = [xcell_position[0], ycell_position[1], z]

    return xcell_position, ycell_position, xycell_position
