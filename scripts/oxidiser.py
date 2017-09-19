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
    
    def __init__(self, crystal, ratio = 2.5, OHratio = 0.66):
        self.crystal = crystal
        self.molecule = crystal.molecule
        
        # #C/O ratio 
        self.ratio = ratio
        self.OHratio = OHratio
        
        self.Ncarbons = self.calc_Ncarbons(crystal)
        self.Nhydrogens = len(crystal.atom_labels) - self.Ncarbons

        self.NO = int(self.Ncarbons / self.ratio)
   
        self.oxidise(crystal, self.NO)
#        self.add_NOH(crystal, self.NOH)
#        self.add_Nepoxy(crystal, self.Nepoxy)
#        self.crystal.natom_types += 2    


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

    def oxidise(self, crystal, NO):
        Ntotal = NO
        N = 0
        a = 0
        b = 0
        OH_attempts = 0
        epoxy_attempts = 0
        while N < Ntotal:
            r = np.random.random()
            if r < self.OHratio:
                OH_attempts += self.find_OH_site(crystal)
                a += 1
            else:
                epoxy_attempts += self.find_epoxy_site(crystal)
                b += 1
            N += 1
        print 'after ',OH_attempts,' attempts, ',a,' OH were added'
        print 'after ',epoxy_attempts,' attempts, ',b,' epoxy were added'


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
                bonded_to = (  self.bonded_to(crystal.bonds,c1) 
                             + self.bonded_to(crystal.bonds,c2))
                mc = np.random.random()
                oxidised_C = [3,6]
                for atom in bonded_to:
                    if crystal.atom_labels[atom -1] in oxidised_C:
                        mc += 0.25
                if mc > 0.99:
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

                hydrogen_flag = False
                near_oxidised = False
                mc = np.random.random()
                oxidised_C = [3,6]
                for atom in bonded_to:
                    if crystal.atom_labels[atom -1] == 2:
                        hydrogen_flag = True
                    if crystal.atom_labels[atom -1] in oxidised_C:
                        mc += 0.34
                if (not hydrogen_flag) and mc > 0.99:
                    self.add_OH(crystal,i)
                    searching = False 
        return attempts
        
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





        
