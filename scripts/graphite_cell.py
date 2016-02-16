import numpy as np
from math import pi,sin,cos
import yaml

class Graphite(object):
    # Handling AB graphite with orthorhombic unit cell (8 atom)
    def __init__(self,CC,layer_gap):
        self.CC = CC
        self.layer_gap = layer_gap
        
    #orthorhombic unitcell lattice parameters
    def cell_shape(self):
        a = 2.0 * self.CC * cos(pi/6.0)
        b = 3.0 * self.CC
        c = 2.0 * self.layer_gap
        cell_dimensions = [a,b,c]
        return cell_dimensions

    def cell_coords(self):
        CC = self.CC
        layer_gap = self.layer_gap
        cos_CC = cos(pi/6.0) * CC
        sin_CC = 0.5 * CC
        C1 = [0,0,0]
        C2 = [0,CC,0]
        C3 = [cos_CC, CC + sin_CC,0]
        C4 = [cos_CC, 2*CC + sin_CC,0]
        C5 = [0,CC, layer_gap]
        C6 = [0,CC*2,layer_gap]
        C7 = [cos_CC, 2*CC + sin_CC, layer_gap]
        C8 = [cos_CC, sin_CC, layer_gap]
        cell_coords = np.array([C1,C2,C3,C4,C5,C6,C7,C8])
        return cell_coords

    def bonds(self,lattice_dimensions):
        internal_bonds = np.array([[1,2],[2,3],[3,4],
                                   [5,6],[6,7],[7,8]])
        bonds = np.empty((0,2))
        
        def add_cross_bond(lattice_dimensions,
                           cell_position,atoms,bonds,i):
            # cell_position : the adjoining cell
            # atoms [i,j]: ith atom in cell, jth atom in adjoining
            atoms[0] += i
            atoms[1] += self.index_cell(cell_position,
                                        lattice_dimensions)
            return np.vstack((bonds,atoms))

        for x in range(lattice_dimensions[0]):
         for y in range(lattice_dimensions[1]):
          for z in range(lattice_dimensions[2]):
            cell_position = [x,y,z]
            
            #Find adjacent cells ===================#
            if x == lattice_dimensions[0]-1:        #
                xcell_position = [0,y,z]            #
            else:                                   #
                xcell_position = [x+1,y,z]          #
                                                    #
            if y == lattice_dimensions[1]-1:        #
                ycell_position = [x,0,z]            #
            else:                                   #
                ycell_position = [x,y+1,z]          #
                                                    #
            xycell_position = [xcell_position[0],   #
                               ycell_position[1],z] #
            #=======================================#
            i = self.index_cell([x,y,z],lattice_dimensions)
            bonds = np.vstack((bonds,internal_bonds+i))

            bonds = add_cross_bond(lattice_dimensions,
                    xcell_position,[3,2],bonds,i)
            bonds = add_cross_bond(lattice_dimensions,
                    xycell_position,[4,1],bonds,i)
            bonds = add_cross_bond(lattice_dimensions,
                    ycell_position,[4,1],bonds,i)
            
            bonds = add_cross_bond(lattice_dimensions,
                    xcell_position,[7,6],bonds,i)
            bonds = add_cross_bond(lattice_dimensions,
                    xcell_position,[8,5],bonds,i)
            bonds = add_cross_bond(lattice_dimensions,
                    ycell_position,[8,5],bonds,i)

        return bonds

    def index_cell(sel,cell_position,lattice_dimensions):
        total = cell_position[2] * 8
        total += cell_position[1] * lattice_dimensions[2] * 8
        total += (cell_position[0] * lattice_dimensions[2] 
                * lattice_dimensions[1] * 8)
        return total

    def assign_molecules(self,lattice_dimensions):
        unit_cell_molecule_label = [1,1,1,1,2,2,2,2]
        molecule_labels = []
        for x in range(lattice_dimensions[0]):
             for y in range(lattice_dimensions[1]):
                 for z in range(lattice_dimensions[2]):
                     labels = np.array(unit_cell_molecule_label)+(z*2)
                     molecule_labels.extend(list(labels))
        return molecule_labels

