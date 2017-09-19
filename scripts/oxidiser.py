import numpy as np

class Oxidiser(object):

    """ 
    Atom types
    1 = C, graphitic (aromatic)
    2 = H, graphitic edge

    3 = C, graphitic C-OH
    4 = O, C-OH
    5 = H, C-OH

    6 = C, epoxy C-O-C 
    7 = O, epoxy 
    """
    
    def __init__(self, crystal, ratio = 2.7, 
                 surface_OHratio = 0.7,
                 edge_OHratio = 0.3,
                 edge_carboxyl_ratio = 0.3):
        self.crystal = crystal
        self.molecule = crystal.molecule
        
        # #C/O ratio 
        self.ratio = ratio
        self.surface_OHratio = surface_OHratio
        self.edge_OHratio = edge_OHratio
        self.edge_carboxyl_ratio = edge_carboxyl_ratio
        
        self.Ncarbons = self.calc_Ncarbons(crystal)
        self.Nhydrogens = len(crystal.atom_labels) - self.Ncarbons

        self.NO = int(self.Ncarbons / self.ratio)
   
        self.oxidise(crystal, self.NO)

        if 2 in self.crystal.atom_labels:
            self.carboxylic = 0.3

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

    def oxidise(self, crystal, Ntotal):
        #edges first
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
        print 'after ',OH_attempts,'\tattempts, ',OH_added,'\tOH were added'
        print 'after ',epoxy_attempts,'\tattempts, ',epoxy_added,'\tepoxy were added'

        e = 0; o = 0
        o,e = self.oxidise_islands(crystal)
        print e+o,'\tisland removed, with ',e,'\tepoxys and ',o,'\tOH'
        print '=========='
        print 'C/O = ',float(self.Ncarbons+carboxyl)/(N+e+o)
        print 'OH/epoxy = ',float(OH_added+o)/(epoxy_added+e)
        print 'Carboxy : OH : epoxy ratio = 1:',float(OH_added+o)/carboxyl,':',float(epoxy_added+e)/carboxyl
        
    def oxidise_islands(self, crystal):
        removed = 100000000
        sp3 = [3,6]
        epoxy_added = 0
        OH_added = 0
        while removed > 0:
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
                  bonded_to = (  self.bonded_to(crystal.bonds,c1+1) 
                               + self.bonded_to(crystal.bonds,c2+1))
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
                    if sp3_flag == 3:
                      self.add_OH(crystal, i)
                      removed += 1
                      OH_added_cycle += 1
            print 'cycle: islands removed', removed,' with ', epoxy_added_cycle,' epoxy and ', OH_added_cycle,' OH'
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
                sp3_C = [3,6]
                edge = [2,8]
                edge_flag =  False
                for atom in bonded_to:
                    if crystal.atom_labels[atom -1] in sp3_C:
                        mc += 0.25
                    if crystal.atom_labels[atom -1] in edge:
                        edge_flag = True
                        break
                if mc > 0.99 and not edge_flag:
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
                sp3_C = [3,6]
                edge = [2,8]
                edge_flag = False
                for atom in bonded_to:
                    if crystal.atom_labels[atom -1] in edge:
                        edge_flag = True 
                        break
                    if crystal.atom_labels[atom -1] in sp3_C:
                        mc += 0.34
                if mc > 0.99 and not edge_flag:
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
        #new_h_angle = np.random.random() * 2 * np.pi
        #new_h_coord = o_coord + np.array([np.sin(h_angle) * 0.7,
        #                              np.cos(h_angle) * 0.7,
        #                              OH])
        h_coord = (  o_coord 
                   + bond_vector * OH 
                   + np.array(([0,0,OH * (-1)**above ] )) )
        
        molecule = crystal.molecule_labels[C_at]
        
        crystal.atom_labels[C_at] = 8 #3 is a C-OH carbon
        crystal.atom_charges[C_at] = 0.2

        crystal.coords = np.vstack((crystal.coords,o_coord))
        crystal.coords[H_at] = h_coord
        crystal.atom_labels += [4]
        crystal.atom_charges += [-0.4]
        crystal.molecule_labels += [molecule]
    
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
        
        crystal.atom_labels[C_at] = 8 #3 is a C-OH carbon
        crystal.atom_charges[C_at] = 0.2

        crystal.coords = np.vstack((crystal.coords,C1_coord))
        crystal.coords = np.vstack((crystal.coords,O1_coord))
        crystal.coords = np.vstack((crystal.coords,O2_coord))
        crystal.coords[H_at] = H_coord
        crystal.atom_labels += [6,4,4]
        crystal.atom_charges += [-0.4,-0.4,0.6]
        crystal.molecule_labels += [molecule]*3 

    def add_OH(self,crystal,at):
        crystal.atom_labels[at] = 3 #3 is a C-OH carbon
        crystal.atom_charges[at] = 0.2
        
        molecule = crystal.molecule_labels[at]

        above = np.random.randint(2)
        CO = 1.4 * (-1)**above
        OH = 0.7 * (-1)**above
        o_coord = crystal.coords[at] + np.array([0,0,CO])
        h_angle = np.random.random() * 2 * np.pi
        h_coord = o_coord + np.array([np.sin(h_angle) * 0.7,
                                      np.cos(h_angle) * 0.7,
                                      OH])
        crystal.coords = np.vstack((crystal.coords,o_coord))
        crystal.coords = np.vstack((crystal.coords,h_coord))
        crystal.atom_labels += [4,5]
        crystal.atom_charges += [-0.4,0.2]
        crystal.molecule_labels += [molecule,molecule]

    def add_epoxy(self,crystal,c1,c2):
        crystal.atom_labels[c1] = 6 #6 is epoxy carbon
        crystal.atom_labels[c2] = 6 #6 is epoxy carbon
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
        crystal.atom_labels += [7]
        crystal.atom_charges += [-0.4]
        crystal.molecule_labels += [molecule]





        
