import numpy as np
from connector import Connector
from write_coords import Writer
from oxidise_rf import init_random_forest
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
    
    def __init__(self, crystal, 
                 ratio = 2.5,               # Target overall C:O ratio
                 surface_OHratio = 0.5,     # Surface OH/epoxy fraction
                 edge_OHratio = 0.25,       # edge H:OH:carboxyl ratio
                 edge_carboxyl_ratio = 0.25,# edge H:OH:carboxyl ratio
                 method='rf',               # which method to calculate a site's affinity 
                                            #empirical / rf (random forest)
                 new_island_freq=0,         # Freq s-1 attempt to add new island 
                 video=False):
        self.crystal = crystal
        self.molecule = crystal.molecule
        
        self.method = method
        self.new_island_freq = new_island_freq
        self.video = video

        self.ratio = ratio
        self.surface_OHratio = surface_OHratio
        self.edge_OHratio = edge_OHratio
        self.edge_carboxyl_ratio = edge_carboxyl_ratio
        
        self.Ncarbons = self.calc_Ncarbons(crystal)
        self.Nhydrogens = len(crystal.atom_labels) - self.Ncarbons

        self.NO = int(self.Ncarbons / self.ratio)
        
        self.CCbonds, self.neighbours = self.neighbour_matrix(self.crystal)
        self.NCCbonds = len(self.CCbonds)
        self.affinities_above, self.affinities_below = self.affinity_matrix(
                                                                self.crystal)
 
        if self.method == 'rf':
            # Read in data from Yang2014, generate random forest regressor
            self.rf = init_random_forest()

        self.atom_states = np.zeros(len(crystal.atom_labels))
        for i in range(len(crystal.atom_labels)):
            if crystal.atom_labels[i] == 2:
                self.atom_states[i] = 3
        
        # lists to record oxidisation process
        self.affinity_order = [0]
        self.time_order = []
        self.time_elapsed_list = []
        self.node_order = []

        with open('affinity.dat','w') as f:
            f.write('#affinity, dt, time_since_island,'+ 
                    'poison_mean, new_islands, available_CC_bonds\n')

        self.oxidise(crystal, self.NO)

        self.crystal.generate_connections()
        self.vdw_defs = {1: 90, # Cg, graphitic (aromatic)
                         2: 91, # Hg, graphitic edge
                         3: 101,# Ct, tertiary C-OH
                         4: 96, # Oa, C-OH
                         5: 97, # Ha, C-OH
                         6: 122,# Oe, epoxy 
                         11: 108,# Cb, Benzyl 
                         7: 109,# Oa, C-OH
                         8: 209, # Cc, Carboxylic carbon
                         9: 210, # Oc, Ketone oxygen
                         10: 211 # Oa, alcohol
                        } # OPLS definitions 
        crystal.vdw_defs = self.vdw_defs

    def oxidise(self, crystal, Ntotal):
        # edges first
        N = 0
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
                else:
                    pass # leave as H

        print 'added ',edge_OH,' OH and ',carboxyl,' COOH at edges'
        
        OH_added = 0
        epoxy_added = 0
        time_elapsed = 0 # since last new island
	dt = 0
        nodes = 0
        new_island = 1
        while N < Ntotal:
            available_CC_bonds = np.sum(np.array(self.affinities_above != 0))
            if not new_island:
                new_island = np.random.poisson(  float(dt) 
                                               * self.new_island_freq
                                               * available_CC_bonds)
            self.node_order += [new_island]
            self.time_elapsed_list += [time_elapsed]

            # choose site
            if new_island:
                dt = 0
                new_island -= 1 
                time_elapsed = 0
                site, above = self.find_new_island()
                atom1,atom2 = self.CCbonds[site] - 1
                state1 = self.atom_states[atom1]
                state2 = self.atom_states[atom2]
                nodes += 1
                print 'new_island accepted,',nodes,'nodes (',new_island,')'
            else:
                site, above, dt = self.find_site()
                time_elapsed += dt
            if above == 0:
                print 'Could not reach C/O ratio:',self.ratio
                print N,' oxygens added'
                break

            # oxygenate at site,above 
            r = np.random.random() # between 0,1
            if r < self.surface_OHratio:
                # add OH
                r2 = np.random.randint(2)
                atom1 = self.CCbonds[site][r2] - 1
                self.add_OH(crystal, above, atom1)
                self.atom_states[atom1] = 1 * above
                self.update_affinity(atom1+1)
                OH_added += 1
            else:
                # add epoxy
                atom1, atom2 = self.CCbonds[site]
                atom1, atom2 = atom1 -1, atom2-1
                self.add_epoxy(crystal, above, atom1, atom2)
                self.atom_states[atom1] = 2 * above
                self.atom_states[atom2] = 2 * above
                self.update_affinity(atom1+1)
                self.update_affinity(atom2+1)
                epoxy_added += 1
            N += 1
            
            # outputs
            if not N % 20:
                print N,'/',Ntotal,'\toxygens added\t',nodes,'nodes'

            if self.video and not N % self.video:
                out = Writer(crystal)
                out.write_xyz(option='a')
                
            with open('affinity.dat','a') as f:
                f.write(str(self.affinity_order[-1])+'\t'+
                        str(dt)+'\t'+
                        str(time_elapsed)+'\t'+
                        str(dt*self.new_island_freq*self.NCCbonds)+'\t'+
                        str(self.node_order[-1])+'\t'+
                        str(available_CC_bonds)+'\n')

        print OH_added,'\tOH were added'
        print epoxy_added,'\tepoxy were added'
        print nodes,'nodes'

        with open('output.dat','w') as f:
            f.write('Cs \t'     +str(self.Ncarbons)+'\n'+
                    'Os \t'     +str(N)+'\n'+
                    'islands \t'+str(nodes)+'\n')

        print '=========='
        print 'C/O = ',float(self.Ncarbons+carboxyl)/(N)
        if epoxy_added != 0:
            print 'OH/epoxy = ',float(OH_added)/(epoxy_added)
        else:
            print 'OH/epoxy = inf'

    def find_12_neighbours(self, crystal, i, j):
        expected_first_neighbours = 4
        expected_second_neighbours = 8

        first_neighbours = self.bonded_to(crystal.bonds, i-1)
        first_neighbours += self.bonded_to(crystal.bonds, j-1)
        first_neighbours = [n+1 for n in first_neighbours]
        first_neighbours = set(first_neighbours) - {i, j}
        if len(first_neighbours) != expected_first_neighbours:
            raise ValueError ('Not enough first neighbours',i,j, first_neighbours)
        
        second_neighbours = set()
        for atom in first_neighbours:
            if crystal.atom_labels[atom-1] == 2:
                expected_second_neighbours -= 2
        for n in first_neighbours:
            second_neighbours  = (second_neighbours | 
                            set(self.bonded_to(crystal.bonds, n-1)) )
        second_neighbours = {n+1 for n in second_neighbours}
        second_neighbours = second_neighbours - first_neighbours - {i,j}
        if len(second_neighbours) != expected_second_neighbours:
            raise ValueError ('Not enough second neighbours',i,j,second_neighbours)
        
        return list(first_neighbours)+list(second_neighbours)

    def affinity_matrix(self, crystal):
        affinities_above = np.ones(self.NCCbonds) 
        affinities_below = np.ones(self.NCCbonds)
        for i in range(self.NCCbonds):
            first_neighbours = self.neighbours[i][0:4]
            for atom in first_neighbours:
                if crystal.atom_labels[atom-1] == 2:
                    affinities_above[i] = 0
                    affinities_below[i] = 0
        return affinities_above, affinities_below

    def neighbour_matrix(self, crystal):
        Nbonds = len(crystal.bonds)
        CCbonds = []
        neighbours = []
        
        for i in range(Nbonds):
            c1 = crystal.bonds[i][0]
            c2 = crystal.bonds[i][1]
            label1 = crystal.atom_labels[ c1-1 ]
            label2 = crystal.atom_labels[ c2-1 ]
            
            if label1 == 1 and label2 == 1:
                CCbonds += [ [c1, c2] ]
                neighbours += [self.find_12_neighbours(crystal, c1, c2) ]
        
        return np.array(CCbonds), neighbours
    
    def update_affinity(self, atom):
        for i in range(self.NCCbonds):
            if atom in self.neighbours[i] and self.affinities_above[i] != 0:
                self.calc_affinities(i)
            elif atom in self.CCbonds[i]:
                self.affinities_above[i] = 0
                self.affinities_below[i] = 0

    def calc_affinities(self, site):
        calc_affinity = getattr(self, 'calc_affinity_' + self.method)
        n = []
        for i in self.neighbours[site]:
            n += [ self.atom_states[i-1] ]
        first = n[0:5]
        second = n [5:]
        above = calc_affinity(first,second) 
        self.affinities_above[site] = above
        first = -np.array(first)
        second = -np.array(second)
        below = calc_affinity(first,second)
        self.affinities_below[site] = below
       
    def calc_affinity_rf(self, first, second):
        edge = False
        X = [0] * 8
        for state in first:
            if state == 1: X[0] += 1
            if state ==-1: X[1] += 1
            if state == 2: X[2] += 1
            if state ==-2: X[3] += 1
        for state in second:
            if state == 1: X[4] += 1
            if state ==-1: X[5] += 1
            if state == 2: X[6] += 1
            if state ==-2: X[7] += 1
            if state == 3: edge = True
        if edge: rate = 1
        else: 
            exponent = self.rf.predict([X])
            rate = 10 ** exponent[0]
        return rate

    def calc_affinity_empirical(self, first, second):
        steric = 0
        polar = 0
        hbond = 0
        edge = 0
        m = [ -3.867, 0.185,
              23.169, -5.138,
              11.648, -4.413]
        for state in first:
            if state == 1: steric += 1
            elif state == 2: steric += 1
            if abs(state) == 1: polar += 1
            if abs(state) == 2: polar += 0.633
            
            #if state == 1: hbond =
        
        for state in second:
            if state == 1: hbond += 1
            if state == 3: edge = 1
        
        steric = (  m[0]    * steric 
                  + m[1] * steric * steric)
        polar  = (  m[2] * polar
                  + m[3] * polar * polar)
        hbond  = (  m[4] * hbond
                  + m[5] * hbond * hbond) 

        if edge:
            rate = 1
        else:
            rate = 10 ** (steric + polar + hbond)
        return rate

    def find_new_island(self):
        # number of sites that are not CH
        total = ( sum(np.array(self.affinities_above != 0)) 
                + sum(np.array(self.affinities_below != 0)) ) 
        r = np.random.random() * total
        R = 0
        above = 0
        for i in range(self.NCCbonds):
            R += bool(self.affinities_above[i])
            if R > r:
                above = 1
                break
        if not above:
            for i in range(self.NCCbonds):
                R += bool(self.affinities_below[i])
                if R > r:
                 above = -1
                 break
        if above == 0:
            #no possible oxidation sites
            raise Exception('Couldnt find a new island site')
            pass 
        first_neighbours = self.neighbours[i][0:4]
        for atom in first_neighbours:
            if self.crystal.atom_labels[atom-1] == 2:
                raise Exception("i've picked an unallowed oxidation site...")
        return i, above

    def find_site(self):
        total = sum(self.affinities_above) + sum(self.affinities_below)
        r = np.random.random() * total
        R = 0
        above = 0
        for i in range(self.NCCbonds):
            R += self.affinities_above[i]
            if R > r:
                above = 1
                self.affinity_order += [self.affinities_above[i]]
                break
        if not above:
            for i in range(self.NCCbonds):
                R += self.affinities_below[i]
                if R > r:
                 self.affinity_order += [self.affinities_below[i]]
                 above = -1
                 break
        if above == 0:
            #no possible oxidation sites
            #raise Exception('Couldnt find a site')
            pass 
        
        time = 1/( total )
        self.time_order += [time]
        return i, above, time


      
    def add_edge_OH(self,crystal, H_at):
        bonded_to = self.bonded_to(crystal.bonds,H_at)
        C_at = bonded_to[0] 
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

        new_bond = np.array(([O_at+1, H_at+1]))
        self.crystal.bonds = np.vstack((crystal.bonds,new_bond))
    
    def add_carboxyl(self, crystal, H_at):
        bonded_to = self.bonded_to(crystal.bonds,H_at)
        C_at = bonded_to[0] 
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
        O2_coord = (  O1_coord + np.array([0,0,-2*sangle*above]) )

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
    
    def add_OH(self,crystal,above, at):
        crystal.atom_labels[at] = 3 #3 is a C-OH carbon
        crystal.atom_charges[at] = 0.265
        molecule = crystal.molecule_labels[at]

        CO = 1.4 * above
        OH = 1.0 * above
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
 


    def add_epoxy(self,crystal,above, c1, c2):
        crystal.atom_labels[c1] = 3 #3 is epoxy carbon
        crystal.atom_labels[c2] = 3 #3 is epoxy carbon
        crystal.atom_charges[c1] = 0.2
        crystal.atom_charges[c2] = 0.2
        
        molecule = crystal.molecule_labels[c1]

        c1c2 = crystal.coords[c2]-crystal.coords[c1]
        if c1c2[0] > 2:
            c1c2[0] += crystal.box_dimensions[0,1]
        elif c1c2[0] < -2:
            c1c2[0] += crystal.box_dimensions[0,1]
        if c1c2[1] > 2:
            c1c2[1] += crystal.box_dimensions[1,1]
        elif c1c2[1] < -2:
            c1c2[1] += crystal.box_dimensions[1,1]
        CO = 0.9 * above
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
        
    def remove_graphitic_bonds(self, crystal, a):
        connections = self.find_connections(crystal.bonds,a+1)
        for i in range(len(connections)):
            bond = connections[i][0]
            if crystal.bond_labels[bond] == 1:
                crystal.bond_labels[bond] = 3

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
        ibonds = self.find_connections(bonds,centre+1)
        bonded_to = []
        for x in ibonds:
            bonded_to += [ bonds[x[0]][x[1]-1] - 1]
        return bonded_to
        

    def find_connections(self,bonds,centre):
        connections = np.where(bonds==centre)
        connections = np.vstack((connections[0],
                                 connections[1]))
        return connections.transpose()  

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
                  bonded_to = ( self.bonded_to(crystal.bonds,c1) 
                               +self.bonded_to(crystal.bonds,c2))
                  sp3_flag = 0
                  for atom in bonded_to:
                      if crystal.atom_labels[atom] in sp3:
                          sp3_flag += 1
                  
                  #if bond is surrounded by sp3 carbons
                  if sp3_flag  >= 3:
                      self.add_epoxy(crystal, c1, c2)
                      removed += 1
                      epoxy_added_cycle += 1
        
            for i in range(len(crystal.atom_labels)):
                if crystal.atom_labels[i] == 1:
                    bonded_to = self.bonded_to(crystal.bonds,i)
                    sp3_flag = 0
                    for atom in bonded_to:
                      if crystal.atom_labels[atom] in sp3:
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

