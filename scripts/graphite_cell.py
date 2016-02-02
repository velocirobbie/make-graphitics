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

    def assign_molecules(self,lattice_dimensions):
        unit_cell_molecule_label = [1,1,1,1,2,2,2,2]
        molecule_labels = []
        for x in range(lattice_dimensions[0]):
             for y in range(lattice_dimensions[1]):
                 for z in range(lattice_dimensions[2]):
                     labels = np.array(unit_cell_molecule_label)+(z*2)
                     molecule_labels.extend(list(labels))
        return molecule_labels
