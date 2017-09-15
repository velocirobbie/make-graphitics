import numpy as np

class Oxidiser(object):
    def __init__(self, crystal, ratio = 6, OH = 0.75):
        self.crystal = crystal
        self.molecule = crystal.molecule
        
        # #C/O ratio 
        self.ratio = ratio
        self.OH = OH
        self.epoxy = 1 - OH
        
        self.Ncarbons = self.calc_Ncarbons(crystal)
        self.Nhydrogens = len(crystal.atom_labels) - self.Ncarbons

        self.NOH = int(self.Ncarbons / self.ratio * self.OH)
        self.Nepoxy = int(self.Ncarbons / self.ratio * self.epoxy)
    
        self.add_NOH(crystal, self.NOH)
        self.add_Nepoxy(crystal, self.Nepoxy)
#        self.crystal.natom_types += 2    


        if 2 in self.crystal.atom_labels:
            self.carboxylic = 0.3

    def calc_Ncarbons(self,crystal):
        N = 0
        for atom_type in crystal.atom_labels:
            if atom_type == 1:
                N += 1
        return N

    def find_connections(self,bonds,centre):
        connections = np.where(bonds==centre)
        connections = np.vstack((connections[0],
                                 connections[1]))
        return connections.transpose()
 
    def add_Nepoxy(self, crystal, N):
        Nbonds = len(crystal.bonds)
        n = 0
        a = 0
        while n < N:
            a +=1
            # pick random atom
            i = np.random.randint(Nbonds)
            c1 = crystal.bonds[i][0]-1
            c2 = crystal.bonds[i][1]-1
            label1 = crystal.atom_labels[ c1 ]
            label2 = crystal.atom_labels[ c2 ]
            # if i is a graphitic bond
            if label1 == 1 and label2 == 1:
                self.add_epoxy(crystal,c1,c2)
                n +=1
        print 'attempts = ',a,' epoxy added = ',n
 
    def add_NOH(self, crystal, N):
       
        Natoms = len(crystal.atom_labels)
        n = 0
        a = 0
        while n < N:
            a +=1
            # pick random atom
            i = np.random.randint(Natoms) 
            # if i is a graphitic carbon
            if crystal.atom_labels[i] == 1:
                # what is i bonded to
                ibonds = self.find_connections(crystal.bonds,i+1)
                bonded_to = []
                for x in ibonds:
                    bonded_to += [ crystal.bonds[x[0]][x[1]-1] ]
                
                hydrogen_flag = False
                hydroxy_flag = False
                for atom in bonded_to:
                    if crystal.atom_labels[atom-1] == 2:
                        hydrogen_flag = True
#                    if crystal.atom_labels[atom-1] == 3:
#                        hydroxy_flag = True
                if (not hydrogen_flag) and (not hydroxy_flag):
                    self.add_OH(crystal,i)
                    n +=1
        print 'attempts = ',a,' OH added = ',n
        
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





        
