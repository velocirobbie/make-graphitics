import numpy as np
from math import pi,sin,cos
import yaml

class Graphene(object):
    def __init__(self,CC,layer_gap=3.354):
        self.CC = CC
        # layer gap is a placeholder, code needs 3 dimensions
        self.layer_gap = layer_gap
        
    #orthorhombic unitcell lattice parameters
    def orthorhombic_graphene_shape(self):
        a = 2.0 * self.CC * cos(pi/6.0)
        b = 3.0 * self.CC
        c = self.layer_gap
        cell_dimensions = [a,b,c]
        return cell_dimensions

    def orthorhombic_graphene_coords(self):
        CC = self.CC
        cos_CC = cos(pi/6.0) * CC
        sin_CC = 0.5 * CC
        C1 = [0,0,0]
        C2 = [0,CC,0]
        C3 = [cos_CC, CC + sin_CC,0]
        C4 = [cos_CC, 2*CC + sin_CC,0]
        cell_coords = np.array([C1,C2,C3,C4])
        return cell_coords
