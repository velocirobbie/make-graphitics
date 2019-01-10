from math import pi,cos,sin,sqrt
import numpy as np
from molecule import Molecule

class Hexagon_Graphene(Molecule):
    def __init__(self,config,forcefield,radius):
        self.CC = config[forcefield]['CC']
        self.CH = config[forcefield]['CH']
        self.layer_gap = config[forcefield]['layer_gap']
        self.order = int((radius + (sqrt(3)/2)*
            (self.CC - self.CH)) / (self.CC * sqrt(3)))
        self.radius = self.order * self.CC *sqrt(3)
        self.natoms = 0
        for ring in range(self.order):
            self.natoms += ((ring + 1) * 2 - 1) * 6
        self.natoms += self.order * 6 # Hydrogens
        
    def cell_shape(self):
        a = 100 + 2 * self.radius
        b = 100 + 2 * self.radius
        c = self.layer_gap
        return [a,b,c]

    def cell_coords(self):
        CC = self.CC
        CH = self.CH
        cos_CC = cos(pi/6.0) * CC
        sin_CC = 0.5 * CC
        
        def add_atom(array,new):
            return np.vstack((array,new))

        section = np.empty((0,3))
        section = add_atom(section,[0,CC,0])
        
        for ring in range(1,self.order):
            section = add_atom(section,
                [-ring * cos_CC,
                 CC + ring * (CC + sin_CC),
                 0])
            for branch in range(ring):
                section = add_atom(section,
                    [-(ring-1) * cos_CC + branch * 2 * cos_CC,
                     sin_CC + ring * (CC + sin_CC),
                     0])
                section = add_atom(section,
                    [-(ring-2) * cos_CC + branch * 2 * cos_CC,
                     CC + ring * (CC + sin_CC),
                     0])
        
        # Add Hydrogens round edge
        for branch in range(self.order):
            section = add_atom(section,
                    [-(self.order-1) * cos_CC + branch * 2 *cos_CC,
                     CC + CH + (self.order-1) * (sin_CC + CC),
                     0])
        
        def rotate_vector(r,theta):
            sint= sin(-theta); cost = cos(theta)
            a = cost * r[0] - sint * r[1]
            b = sint * r[0] + cost * r[1]
            c = r[2]
            return [a,b,c]
        
        coords = np.empty((0,3))
        for theta in range(6):
            for atom in section:
                coord = rotate_vector(atom,theta*pi/3.0)
                coords = add_atom(coords,coord)

        return coords

    def assign_molecules(self,lattice_dimensions):
        molecule_labels = []
        for z in range(lattice_dimensions[2]):
            labels = np.ones(self.natoms,dtype=int)
            molecule_labels.extend(list(labels+z))
        return molecule_labels

    def assign_atom_labels(self,lattice_dimensions):
        atom_labels = []
        segment_carbons = [1] * ((self.natoms - (self.order*6))/6)
        segment_hydrogens = [2] * (self.order)
        for z in range(lattice_dimensions[2]):
            for i in range(6):
               atom_labels += segment_carbons + segment_hydrogens
        return atom_labels

    def assign_atom_charges(self,lattice_dimensions,q):
        atom_charges = []
        number_of_inner_carbons = 0
        for ring in range(self.order -1):
            number_of_inner_carbons += ((ring + 1) * 2 - 1)
        segment_inner_carbons = [0] * number_of_inner_carbons
        segment_outer_carbons = [-q] + [0,-q] * (self.order-1)
        segment_hydrogens = [q] * self.order
        segment_charges = (segment_inner_carbons +
                           segment_outer_carbons +
                           segment_hydrogens)
        for z in range(lattice_dimensions[2]):
            for i in range(6):
                atom_charges += segment_charges
        return atom_charges

    def assign_bonds(self,lattice_dimensions):
        segment_bonds = []
        natoms_section = self.natoms/6
        natoms_passed = 0
        for ring in range(1,self.order):
            natoms_ring = (ring+1) * 2 - 1
            for branch in range(ring):
                b1 = natoms_passed + 1 + branch * 2
                b2 = b1 + natoms_ring -1
                segment_bonds += [[b1,b2]]
                segment_bonds += [[b2,b2-1]]
                segment_bonds += [[b2,b2+1]]
            natoms_passed += ring * 2 - 1
 
        for i in range(self.order):
            H = natoms_section - i
            C = natoms_section - self.order - (i * 2)
            segment_bonds += [[C,H]]
        
        last_on_ring = 0
        for i in range(self.order):
            last_on_ring += i * 2 + 1
            first_on_ring = last_on_ring - (i * 2)
            segment_bonds += [[last_on_ring,
                first_on_ring+natoms_section]]
       
        segment_bonds = np.array(segment_bonds)
        molecule_bonds = np.empty((0,2),dtype=int)
        for i in range(6):
            molecule_bonds = np.vstack((molecule_bonds,
                segment_bonds + natoms_section * i))
        
        for connection in range(self.order):
            molecule_bonds[-(connection+1),1] -= self.natoms 
        
        bonds = np.empty((0,2),dtype=int)
        for z in range(lattice_dimensions[2]):
            bonds = np.vstack((bonds,
                molecule_bonds+(z*self.natoms)))

        return bonds

    def connection_types(self):
        bond_types = [[1,1],[1,2]]
        angle_types = [[1,1,1],[1,1,2]]
        dihedral_types = [[1,1,1,1],[1,1,1,2],[2,1,1,2]]
        improper_types = [[1,1,1,1],[1,1,1,2]]
        return bond_types, angle_types, dihedral_types, improper_types





