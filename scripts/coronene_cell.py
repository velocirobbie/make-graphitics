from math import pi,cos,sin
import numpy as np

class Coronene(object):
    def __init__(self,CC,CH,layer_gap):
        self.CC = CC
        self.CH = CH
        self.layer_gap = layer_gap

    def cell_shape(self):
        a = 100
        b = 100
        c = self.layer_gap
        return [a,b,c]

    def cell_coords(self):
        CC = self.CC
        CH = self.CH
        cos_CC = cos(pi/6.0) * CC
        sin_CC = 0.5 * CC
        C1 = [0,CC,0]
        C2 = [0,2 * CC,0]
        C3 = [-cos_CC,2.5 * CC,0]
        C4 = [+cos_CC,2.5 * CC,0]
        H1 = [-cos_CC,2.5 * CC + CH,0]
        H2 = [+cos_CC,2.5 * CC + CH,0]
        section = np.array((C1,C2,C3,C4,H1,H2))
        coords = np.empty((0,3))
        
        def rotate_vector(r,theta):
            sint= sin(theta); cost = cos(theta)
            a = cost * r[0] - sint * r[1]
            b = sint * r[0] + cost * r[1]
            c = r[2]
            return [a,b,c]

        for theta in range(6):
            for atom in section:
                coord = rotate_vector(atom,theta*pi/3.0)
                coords = np.vstack((coords,coord))
            
        return coords

    def assign_molecules(self,lattice_dimensions):
        molecule_labels = []
        for z in range(lattice_dimensions[2]):
            labels = np.ones((36),dtype=int)
            molecule_labels.extend(list(labels+z))
        return molecule_labels

    def assign_atom_labels(self,lattice_dimensions):
        atom_labels = []
        segment_lables = [1,1,1,1,2,2]
        for z in range(lattice_dimensions[2]):
            for i in range(6):
                atom_labels.extend(list(segment_lables))
        return atom_labels

    def assign_atom_charges(self,lattice_dimensions,q):
        atom_charges = []
        segment_charges = [0,0,-q,-q,q,q]
        for z in range(lattice_dimensions[2]):
            for i in range(6):
                atom_charges.extend(list(segment_charges))
        return atom_charges


