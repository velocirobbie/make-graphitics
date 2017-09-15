import numpy as np
from math import pi,sin,cos
import yaml
from molecule import Molecule

class Graphene(Molecule):
    def __init__(self,config,forcefield):
        self.CC = config[forcefield]['CC']
        self.layer_gap = config[forcefield]['layer_gap']

    def cell_shape(self):
        a = 2.0 * self.CC * cos(pi/6.0)
        b = 3.0 * self.CC
        c = self.layer_gap
        cell_dimensions = [a,b,c]
        return cell_dimensions

    def cell_coords(self):
        CC = self.CC
        cos_CC = cos(pi/6.0) * CC
        sin_CC = 0.5 * CC
        C1 = [0,0,0]
        C2 = [0,CC,0]
        C3 = [cos_CC, CC + sin_CC,0]
        C4 = [cos_CC, 2*CC + sin_CC,0]
        cell_coords = np.array([C1,C2,C3,C4])
        return cell_coords

    def assign_molecules(self,lattice_dimensions):
        unit_cell_molecule_label = [1,1,1,1]
        molecule_labels = []
        for x in range(lattice_dimensions[0]):
             for y in range(lattice_dimensions[1]):
                 for z in range(lattice_dimensions[2]):
                     labels = np.array(unit_cell_molecule_label)+(z*2)
                     molecule_labels.extend(list(labels))
        return molecule_labels

    def assign_atom_labels(self,lattice_dimensions):
        atom_labels = []
        cell_labels = [1,1,1,1]
        for x in range(lattice_dimensions[0]):
             for y in range(lattice_dimensions[1]):
                 for z in range(lattice_dimensions[2]):
                     atom_labels.extend(list(cell_labels))
        return atom_labels

    def assign_atom_charges(self,lattice_dimensions,q):
        atom_charges = []
        cell_charges = [0,0,0,0]
        for x in range(lattice_dimensions[0]):
             for y in range(lattice_dimensions[1]):
                 for z in range(lattice_dimensions[2]):
                    atom_charges.extend(list(cell_charges))
        return atom_charges

    def assign_bonds(self,lattice_dimensions):
        internal_bonds = np.array([[1,2],[2,3],[3,4]])
        bonds = np.empty((0,2),dtype=int)
        # loop through all cells
        for x in range(lattice_dimensions[0]):
         for y in range(lattice_dimensions[1]):
          for z in range(lattice_dimensions[2]):
            cell_position = [x,y,z]
            xcell_position,ycell_position,xycell_position = (
             find_adjacent_cells(cell_position,lattice_dimensions))
            
            # Add bonds within cell
            i = self.index_cell([x,y,z],lattice_dimensions,4)
            bonds = np.vstack((bonds,internal_bonds+i))
            # Add bonds that cross cell boundries
            bonds = self.add_cross_bond(lattice_dimensions,
                    xcell_position,[3,2],bonds,i)
            bonds = self.add_cross_bond(lattice_dimensions,
                    xycell_position,[4,1],bonds,i)
            bonds = self.add_cross_bond(lattice_dimensions,
                    ycell_position,[4,1],bonds,i)
            
        return bonds

    def index_cell(self,cell_position,lattice_dimensions,
                   atoms_per_cell):
        N = atoms_per_cell
        total = cell_position[2] * N
        total += cell_position[1] * lattice_dimensions[2] * N
        total += (cell_position[0] * lattice_dimensions[2] 
                * lattice_dimensions[1] * N)
        return total

    def add_cross_bond(self,lattice_dimensions,
                   cell_position,atoms,bonds,i):
        # cell_position : the adjoining cell
        # atoms [i,j]: ith atom in cell, jth atom in adjoining
        atoms[0] += i
        atoms[1] += self.index_cell(cell_position,
                                   lattice_dimensions,4)
        return np.vstack((bonds,atoms))

    def connection_types(self):
        bond_types = [[1,1]]
        angle_types = [[1,1,1]]
        torsion_types = [[1,1,1,1]]
        improper_types = [[1,1,1,1]]
        return bond_types, angle_types, torsion_types, improper_types

def find_connections(bonds,centre):
    connections = np.where(bonds==centre)
    connections =np.vstack((connections[0],connections[1]))
    return connections.transpose()
            
def find_neighbours(bonds,centre):
    connections = find_connections(bonds,centre)
    neighbours = []
    for connection in connections:
        # Find atom connected to centre
        neighbour = bonds[connection[0]][connection[1]-1]
        neighbours.append(neighbour)
    return neighbours

def find_adjacent_cells(cell_position,lattice_dimensions):
    x,y,z = cell_position
    # Account for periodic system 
    if x == lattice_dimensions[0]-1: 
        xcell_position=[0,y,z]
    else: xcell_position = [x+1,y,z]
    
    if y == lattice_dimensions[1]-1: 
        ycell_position=[x,0,z]
    else: ycell_position = [x,y+1,z]          
                                                    
    xycell_position = [xcell_position[0],   
                       ycell_position[1],z]

    return xcell_position, ycell_position, xycell_position
 
