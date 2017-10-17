import numpy as np
from connector import Connector

class Oxidiser(object):

    """ 
    Atom types
    1 = Cg, graphitic (aromatic)
    2 = Hg, graphitic edge

    3 = Ct, tertiary C-OH
    4 = Oa, C-OH
    5 = Ha, C-OH

    3 = Ct, epoxy C-O-C 
    6 = Oe, epoxy 
    
        Edges
    11 = Cb, Benzyl 
    7 = Oa, C-OH
    5 = Ha, C-OH

    11 = Cb, Benzyl carbon
    8 = Cc, Carboxylic carbon
    9 = Oc, Ketone oxygen
    10 = Oa, alcohol
    5   = Ha, alcohol
    """
    
    def __init__(self, crystal, ratio = 2.6,
                 affinity = 0.999,
                 surface_OHratio = 0.5,
                 edge_OHratio = 0.25,
                 edge_carboxyl_ratio = 0.25):
        self.crystal = crystal
        self.molecule = crystal.molecule
        
        self.affinity = affinity
    
        # C/O ratio 
        self.ratio = ratio
        self.surface_OHratio = surface_OHratio
        self.edge_OHratio = edge_OHratio
        self.edge_carboxyl_ratio = edge_carboxyl_ratio
        
        self.Ncarbons = self.calc_Ncarbons(crystal)
        self.Nhydrogens = len(crystal.atom_labels) - self.Ncarbons

        self.NO = int(self.Ncarbons / self.ratio)
   
        self.oxidise(crystal, self.NO)
#        crystal.bond_labels = [1] * len(crystal.bonds)
        self.generate_connections(crystal)
        self.vdw_defs = {1: 90, # Cg, graphitic (aromatic)
                         2: 91, # Hg, graphitic edge
                         3: 101,# Ct, tertiary C-OH
                         4: 96, # Oa, C-OH
                         5: 97, # Ha, C-OH
                         # 3 = Ct, epoxy C-O-C 
                         6: 122,# Oe, epoxy 
                         11: 108,# Cb, Benzyl 
                         7: 109,# Oa, C-OH
                         #5 = Ha, C-OH
                         #11 = Cb, Benzyl carbon
                         8: 209, # Cc, Carboxylic carbon
                         9: 210, # Oc, Ketone oxygen
                         10: 211 # Oa, alcohol
                         #5   = Ha, alcohol
                         }


    def generate_connections(self,crystal):
        connect = Connector()
        crystal.bond_types = connect.find_bond_types(crystal.atom_labels,
                                                        crystal.bonds)
        crystal.bond_labels = connect.bond_labels(
                crystal.atom_labels,crystal.bonds,
                crystal.bond_types)
        #
        crystal.angles = connect.angles(crystal.bonds)
        crystal.angle_types = connect.find_angle_types(crystal.atom_labels,
                                                        crystal.angles)
        crystal.angle_labels = connect.angle_labels(
                crystal.atom_labels,crystal.angles,
                crystal.angle_types)
        #
        crystal.torsions = connect.torsions(crystal.bonds) 
        crystal.torsion_types = connect.find_torsion_types(crystal.atom_labels,
                                                           crystal.torsions)
        crystal.torsion_labels = connect.torsion_labels(
                crystal.atom_labels,crystal.torsions,
                crystal.torsion_types)
        #
        crystal.impropers = connect.impropers(crystal.bonds) 
        crystal.improper_types = connect.find_improper_types(crystal.atom_labels,
                                                           crystal.impropers)
        crystal.improper_labels = connect.improper_labels(
                crystal.atom_labels,crystal.impropers,
                crystal.improper_types)

    def find_12_neighbours(self, crystal, i, j):
        first_neighbours = self.bonded_to(crystal.bonds, i)
        first_neighbours += self.bonded_to(crystal.bonds, j)
        first_neighbours = set(first_neighbours) - {i, j}
        if len(first_neighbours) != 4:
            raise ValueError ('Not enough first neighbours',i,j, first_neighbours)
        second_neighbous = {}
        for n in first_neighbours:
            second_neighbours  += set(self.bonded_to(crystal.bonds, n))
        second_neighbours = second_neighbours - first_neighbours
        if len(second_neighbours) != 8:
            raise ValueError ('Not enough second neighbours',i,j,second_neighbours)



    def bonded_to(self, bonds, centre):
        ibonds = self.find_connections(bonds,centre)
        bonded_to = []
        for x in ibonds:
            bonded_to += [ bonds[x[0]][x[1]-1]]
        return bonded_to
        

    def find_connections(self,bonds,centre):
        connections = np.where(bonds==centre)
        connections = np.vstack((connections[0],
                                 connections[1]))
        return connections.transpose()  

    def neighbour_matrix(self, crystal):
        Nbonds = len(crystal.bonds)
        neighbours = np.empty(Nbonds,14)
        for i in range(Nbonds):
            c1 = crystal.bonds[i][0]-1
            c2 = crystal.bonds[i][1]-1
            label1 = crystal.atom_labels[ c1 ]
            label2 = crystal.atom_labels[ c2 ]

            if label1 == 1 and label2 == 1:
                self.find_14_neighbours(crystal.bonds, c1, c2)
 
    def affinity_matrix(self, crystal):
        Nbonds = len(crystal.bonds)
        affinity = np.zeros(Nbonds)
        for i in range(Nbonds):
            c1 = crystal.bonds[i][0]-1
            c2 = crystal.bonds[i][1]-1
            label1 = crystal.atom_labels[ c1 ]
            label2 = crystal.atom_labels[ c2 ]
            
            # if i is a graphitic bond
            if label1 == 1 and label2 == 1:
                # is it near oxidised sections
                bonded_to = (  self.bonded_to(crystal.bonds,c1+1) 
                             + self.bonded_to(crystal.bonds,c2+1))
                mc = np.random.random()
                sp3_C = [3]
                edge = [2,8]
                edge_flag =  False
                for atom in bonded_to:
                    if crystal.atom_labels[atom -1] in sp3_C:
                        mc += 0.25
                    if crystal.atom_labels[atom -1] in edge:
                        edge_flag = True
                        break
                if mc > self.affinity and not edge_flag:
                    self.add_epoxy(crystal,c1,c2)
                    searching = False
        return attempts



    def oxidise(self, crystal, Ntotal):
        # edges first
        N = 0
        i = 0
        edge_OH = 0
        carboxyl = 0
        for i in range(len(crystal.atom_labels)):
            if crystal.atom_labels[i] == 2:
                r = np.random.random()
                if r < self.edge_OHratio:
                    self.add_edge_OH(crystal,i)
                    N += 1
                    edge_OH += 1
                elif r > 1 - self.edge_carboxyl_ratio:
                    self.add_carboxyl(crystal,i)
                    N += 2
                    carboxyl += 1
        
        print 'added ',edge_OH,' OH and ',carboxyl,' COOH at edges'
        
        OH_added = 0
        epoxy_added = 0
        OH_attempts = 0
        epoxy_attempts = 0 
        while N < Ntotal:
            r = np.random.random()
            if r < self.surface_OHratio:
                OH_attempts += self.find_OH_site(crystal)
                OH_added += 1
            else:
                epoxy_attempts += self.find_epoxy_site(crystal)
                epoxy_added += 1
            N += 1
        print 'after ',OH_attempts,'\tattempts, ',\
                OH_added,'\tOH were added'
        print 'after ',epoxy_attempts,'\tattempts, ',\
                epoxy_added,'\tepoxy were added'

        e = 0; o = 0
        o,e = self.oxidise_islands(crystal)
        print e+o,' island removed, with ',e,' epoxys and ',o,' OH'
        print '=========='
        print 'C/O = ',float(self.Ncarbons+carboxyl)/(N+e+o)
        print 'OH/epoxy = ',float(OH_added+o)/(epoxy_added+e)
        print 'Carboxy : OH : epoxy ratio = 1:',\
              #float(OH_added+o)/carboxyl,':',\
              #float(epoxy_added+e)/carboxyl
        
    def oxidise_islands(self, crystal):
        removed = 1 # not 0
        sp3 = [3]
        epoxy_added = 0
        OH_added = 0
        while removed != 0:
            epoxy_added_cycle = 0
            OH_added_cycle = 0
            removed = 0
            for bond in crystal.bonds:
                c1 = bond[0]-1
                c2 = bond[1]-1
                label1 = crystal.atom_labels[ c1 ]
                label2 = crystal.atom_labels[ c2 ]
                # if i is a graphitic bond
                if label1 == 1 and label2 == 1:
                  # is it near oxidised sections
                  bonded_to = ( self.bonded_to(crystal.bonds,c1+1) 
                               +self.bonded_to(crystal.bonds,c2+1))
                  sp3_flag = 0
                  for atom in bonded_to:
                      if crystal.atom_labels[atom -1] in sp3:
                          sp3_flag += 1
                  
                  #if bond is surrounded by sp3 carbons
                  if sp3_flag  >= 3:
                      self.add_epoxy(crystal, c1, c2)
                      removed += 1
                      epoxy_added_cycle += 1
        
            for i in range(len(crystal.atom_labels)):
                if crystal.atom_labels[i] == 1:
                    bonded_to = self.bonded_to(crystal.bonds,i+1)
                    sp3_flag = 0
                    for atom in bonded_to:
                      if crystal.atom_labels[atom -1] in sp3:
                           sp3_flag += 1
                    # if bond is surrounded by sp3 carbons
                    if sp3_flag >= 3:
                      self.add_OH(crystal, i)
                      removed += 1
                      OH_added_cycle += 1

            print 'cycle: islands removed', removed,' with ',\
                   epoxy_added_cycle,' epoxy and ',\
                   OH_added_cycle,' OH'
            OH_added += OH_added_cycle
            epoxy_added += epoxy_added_cycle

        return OH_added, epoxy_added
 
    def find_epoxy_site(self, crystal):
        Nbonds = len(crystal.bonds)
        searching = True
        attempts = 0
        while searching:
            attempts += 1
            # pick random atom
            i = np.random.randint(Nbonds)
            c1 = crystal.bonds[i][0]-1
            c2 = crystal.bonds[i][1]-1
            label1 = crystal.atom_labels[ c1 ]
            label2 = crystal.atom_labels[ c2 ]
            # if i is a graphitic bond
            if label1 == 1 and label2 == 1:
                # is it near oxidised sections
                bonded_to = (  self.bonded_to(crystal.bonds,c1+1) 
                             + self.bonded_to(crystal.bonds,c2+1))
                mc = np.random.random()
                sp3_C = [3]
                edge = [2,8]
                edge_flag =  False
                for atom in bonded_to:
                    if crystal.atom_labels[atom -1] in sp3_C:
                        mc += 0.25
                    if crystal.atom_labels[atom -1] in edge:
                        edge_flag = True
                        break
                if mc > self.affinity and not edge_flag:
                    self.add_epoxy(crystal,c1,c2)
                    searching = False
        return attempts

    def find_OH_site(self, crystal):
        Natoms = len(crystal.atom_labels)
        searching = True
        attempts = 0
        while searching:
            attempts += 1
            # pick random atom
            i = np.random.randint(Natoms) 
            # if i is a graphitic carbon
            if crystal.atom_labels[i] == 1:
                # what is i bonded to
                bonded_to = self.bonded_to(crystal.bonds,i+1)

                mc = np.random.random()
                sp3_C = [3]
                edge = [2,8]
                edge_flag = False
                for atom in bonded_to:
                    if crystal.atom_labels[atom -1] in edge:
                        edge_flag = True 
                        break
                    if crystal.atom_labels[atom -1] in sp3_C:
                        mc += 0.34
                if mc > self.affinity and not edge_flag:
                    searching = False 
                    self.add_OH(crystal,i)
        return attempts
       
    def add_edge_OH(self,crystal, H_at):
        bonded_to = self.bonded_to(crystal.bonds,H_at+1)
        C_at = bonded_to[0] - 1
        if len(bonded_to) != 1: raise ValueError
        
        C_coord = crystal.coords[C_at]
        H_coord = crystal.coords[H_at]
        CO = 1.4 
        OH = 0.7 
        bond_vector = H_coord - C_coord
        bond_vector = bond_vector / np.linalg.norm(bond_vector)
        o_coord = C_coord + bond_vector * CO
        above = np.random.randint(2)        
        h_coord = (  o_coord 
                   + bond_vector * OH 
                   + np.array(([0,0,OH * (-1)**above ] )) )
        
        molecule = crystal.molecule_labels[C_at]
        O_at = H_at # H becomes O to preserve bond already there
        H_at = len(crystal.atom_labels)
        crystal.atom_labels[C_at] = 11 #3 is a C-OH carbon
        crystal.atom_charges[C_at] = 0.15
        crystal.coords[O_at] = o_coord
        crystal.atom_labels[O_at] = 7
        crystal.atom_charges[O_at] = -0.585
 
        crystal.coords = np.vstack((crystal.coords,h_coord))
        crystal.atom_labels += [5] # H
        crystal.atom_charges += [0.435]
        crystal.molecule_labels += [molecule]

#        self.change_bond_label(crystal,C_at,O_at,8)
        new_bond = np.array(([O_at+1, H_at+1]))
        self.crystal.bonds = np.vstack((crystal.bonds,new_bond))
#        crystal.bond_labels += [5]
    
    def add_carboxyl(self, crystal, H_at):
        bonded_to = self.bonded_to(crystal.bonds,H_at+1)
        C_at = bonded_to[0] - 1
        if len(bonded_to) != 1: raise ValueError
        
        C_coord = crystal.coords[C_at]
        H_coord = crystal.coords[H_at]
        CC = 1.4
        CO = 1.4 
        OH = 1.1 
        angle = np.pi / 3
        sangle = np.sin(angle) * CO
        cangle = np.cos(angle) * CO
        bond_vector = H_coord - C_coord
        bond_vector = bond_vector / np.linalg.norm(bond_vector)
        
        C1_coord = C_coord + bond_vector * CC
        above = (-1)**( np.random.randint(2) )
        O1_coord = (  C1_coord 
                    + bond_vector * cangle
                    + np.array([0,0,sangle*above]) )
        O2_coord = O1_coord * np.array([1,1,-1])
        H_coord  = O2_coord + bond_vector * OH
        
        molecule = crystal.molecule_labels[C_at]

        C1_at = H_at
        O1_at = len(crystal.atom_labels)
        O2_at = O1_at + 1
        H_at = O1_at + 2

        crystal.atom_labels[C_at] = 11 
        crystal.atom_charges[C_at] = -0.115

        crystal.coords[C1_at] = C1_coord
        crystal.atom_labels[C1_at] = 8
        crystal.atom_charges[C1_at] = 0.635
        crystal.coords = np.vstack((crystal.coords,O1_coord))
        crystal.coords = np.vstack((crystal.coords,O2_coord))
        crystal.coords = np.vstack((crystal.coords,H_coord))
        crystal.atom_labels += [9,10,5]
        crystal.atom_charges += [-0.44,-0.53,0.45]
        crystal.molecule_labels += [molecule]*3 

        new_bonds = np.array(([C1_at+1,O1_at+1],
                              [C1_at+1,O2_at+1],
                              [O2_at+1,H_at+1]))
        crystal.bonds = np.vstack((crystal.bonds,new_bonds))
#        crystal.bond_labels += [10,11,5]
#        self.change_bond_label(crystal,C_at,C1_at,9)
    
    def add_OH(self,crystal,at):
        crystal.atom_labels[at] = 3 #3 is a C-OH carbon
        crystal.atom_charges[at] = 0.265
        molecule = crystal.molecule_labels[at]

        above = np.random.randint(2)
        CO = 1.4 * (-1)**above
        OH = 1.0 * (-1)**above
        o_coord = crystal.coords[at] + np.array([0,0,CO])
        h_coord = o_coord + np.array([0,0,OH])
        crystal.coords = np.vstack((crystal.coords,o_coord))
        crystal.coords = np.vstack((crystal.coords,h_coord))
        crystal.atom_labels += [4,5]
        crystal.atom_charges += [-0.683,0.418]
        crystal.molecule_labels += [molecule,molecule]
        
        hid = len(crystal.atom_labels)
        oid = hid -1
        new_bonds = np.array(([at+1, oid],
                              [oid, hid]))
        crystal.bonds = np.vstack((crystal.bonds,new_bonds))
#        crystal.bond_labels += [4,5]
#        self.remove_graphitic_bonds(crystal, at)
 


    def add_epoxy(self,crystal,c1,c2):
        crystal.atom_labels[c1] = 3 #3 is epoxy carbon
        crystal.atom_labels[c2] = 3 #3 is epoxy carbon
        crystal.atom_charges[c1] = 0.2
        crystal.atom_charges[c2] = 0.2
        
        molecule = crystal.molecule_labels[c1]

        above = np.random.randint(2)
        c1c2 = crystal.coords[c2]-crystal.coords[c1]
        CO = 0.9 * (-1)**above
        o_coord = ( crystal.coords[c1] 
                    + np.array([0,0,CO])
                    + 0.5 * c1c2 
                    )
                                        
        crystal.coords = np.vstack((crystal.coords,o_coord))
        crystal.atom_labels += [6]
        crystal.atom_charges += [-0.4]
        crystal.molecule_labels += [molecule]
        oid = len(crystal.atom_labels)
        new_bonds = np.array(([c1+1, oid],
                              [c2+1, oid]))
        crystal.bonds = np.vstack((crystal.bonds,new_bonds))
#        crystal.bond_labels += [6,6]
#        self.change_bond_label(crystal, c1, c2, 7)
#        self.remove_graphitic_bonds(crystal, c1)
#        self.remove_graphitic_bonds(crystal, c2)
        
    def remove_graphitic_bonds(self, crystal, a):
        connections = self.find_connections(crystal.bonds,a+1)
        for i in range(len(connections)):
            bond = connections[i][0]
            if crystal.bond_labels[bond] == 1:
                crystal.bond_labels[bond] = 3
#            elif crystal.bond_labels[bond] =

    def change_bond_label(self, crystal, a1, a2, label):
        connections = self.find_connections(crystal.bonds,a1+1)
        for i in range(len(connections)):
            bond = connections[i][0]
            if a2+1 in crystal.bonds[bond]:
                crystal.bond_labels[bond] = label

    def calc_Ncarbons(self,crystal):
        N = 0
        for atom_type in crystal.atom_labels:
            if atom_type == 1:
                N += 1
        return N

    def bonded_to(self, bonds, centre):
        ibonds = self.find_connections(bonds,centre)
        bonded_to = []
        for x in ibonds:
            bonded_to += [ bonds[x[0]][x[1]-1]]
        return bonded_to
        

    def find_connections(self,bonds,centre):
        connections = np.where(bonds==centre)
        connections = np.vstack((connections[0],
                                 connections[1]))
        return connections.transpose()  

