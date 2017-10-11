import yaml
from opls_reader import OPLS_Reader

class Parameterise(object):
    def __init__(self, crystal, vdw_defs, forcefield = 'OPLS'):
        
        self.vdw_defs = vdw_defs
        
        if forcefield == 'OPLS':
            paramfile = 'params/oplsaa.prm'
            data = OPLS_Reader(paramfile)
        else:
            raise TypeError(forcefield,'forcefield not implemented')
        self.vdw_type = data.vdw_type
        self.bond_data = data.bond
        self.angle_data = data.angle
        self.torsion_data = data.torsion
        self.improper_data = data.improper
        self.pair_data = data.pair
        self.mass_data = data.mass
        self.charge_data = data.charge
        
        self.type_defs = {}
        for label in self.vdw_defs:
            vdw = self.vdw_defs[label]
            self.type_defs[label] = self.vdw_type['type'][ 
                                            self.vdw_type['vdw'].index(vdw) ]
        print 'Atom label -> OPLS vdw definitions: \t',  self.vdw_defs
        print 'Atom label -> OPLS type definitions: \t', self.type_defs

        crystal.bond_coeffs = self.match_bonds(crystal.bond_types)
        crystal.angle_coeffs = self.match_angles(crystal.angle_types)
        crystal.torsion_coeffs = self.match_torsions(crystal.torsion_types)
        crystal.improper_coeffs = self.match_impropers(crystal.improper_types)
        crystal.pair_coeffs = self.match_pairs()
        crystal.masses = self.match_masses()

    def match_charges(self):
        charge_data = self.charge_data
        charge_coeffs = {'a':[], 'q':[]}
        for label in self.vdw_defs:
            found = 0
            for j in range(len(charge_data['a'])):
                atom_data = charge_data['a'][j] 
                if self.vdw_defs[label] == atom_data:
                    found += 1
                    charge_coeffs['a'] += [label]
                    charge_coeffs['q']    += [ charge_data['q'][j] ]
            if found != 1:
                raise ValueError('WRONG',label,'\t found ',found,' entries')
        return charge_coeffs



    def match_masses(self):
        mass_data = self.mass_data
        mass_coeffs = {'a':[], 'm':[]}
        for label in self.vdw_defs:
            found = 0
            for j in range(len(mass_data['a'])):
                atom_data = mass_data['a'][j] 
                if self.vdw_defs[label] == atom_data:
                    found += 1
                    mass_coeffs['a'] += [label]
                    mass_coeffs['m']    += [ mass_data['m'][j] ]
            if found != 1:
                raise ValueError('WRONG',label,'\t found ',found,' entries')
        return mass_coeffs

        
    def match_pairs(self):
        pair_data = self.pair_data
        pair_coeffs = {'a':[], 'e':[], 's':[]}
        for label in self.vdw_defs:
            found = 0
            for j in range(len(pair_data['a'])):
                atom_data = pair_data['a'][j] 
                if self.vdw_defs[label] == atom_data:
                    found += 1
                    pair_coeffs['a'] += [label]
                    pair_coeffs['e']    += [ pair_data['e'][j] ]
                    pair_coeffs['s']    += [ pair_data['s'][j] ]
            if found != 1:
                raise ValueError('WRONG',label,'\t found ',found,' entries')
        return pair_coeffs



    def match_bonds(self, bond_types):
        bond_data = self.bond_data
        bond_coeffs = {'type':[], 'k':[], 'r':[]}
        for i in range(len(bond_types)):
            a1 = self.type_defs[ bond_types[i][0] ]
            a2 = self.type_defs[ bond_types[i][1] ]
            atoms = [a1,a2]
            found = 0
            for j in range(len(bond_data['k'])):
                atom_data =[ bond_data['a1'][j], bond_data['a2'][j] ]
                flag1 = atoms==atom_data
                flag2 = atoms==list(reversed(atom_data))
                if flag1 or flag2:
                    found += 1
                    bond_coeffs['type'] += [i+1]
                    bond_coeffs['k']    += [ bond_data['k'][j] ]
                    bond_coeffs['r']    += [ bond_data['r'][j] ]
            if found != 1:
                raise ValueError('WRONG',atoms,'\t found ',found,' entries')
        return bond_coeffs

    def match_angles(self, angle_types):

        def search_angles(angle_data, atoms, angle_coeffs):
            found = 0
            for j in range(len(angle_data['k'])):
                atom_data =[ angle_data['a1'][j], 
                             angle_data['a2'][j],
                             angle_data['a3'][j] ]
                flag1 = atoms==atom_data
                flag2 = atoms==list(reversed(atom_data))
                if flag1 or flag2:
                    found += 1
                    angle_coeffs['type'] += [i+1]
                    angle_coeffs['k']    += [ angle_data['k'][j] ]
                    angle_coeffs['r']    += [ angle_data['r'][j] ]
            return angle_coeffs, found

        angle_data = self.angle_data
        angle_coeffs = {'type':[], 'k':[], 'r':[]}
        questionable_substitutions = 0
        for i in range(len(angle_types)):
            a1 = self.type_defs[ angle_types[i][0] ]
            a2 = self.type_defs[ angle_types[i][1] ]
            a3 = self.type_defs[ angle_types[i][2] ]
            atoms = [a1,a2,a3]
            #print atoms, angle_types[i]
            angle_coeffs, found = search_angles(angle_data, atoms, angle_coeffs)

            if found == 0:
                for j in range(3):
                    if atoms[j] == 48: 
                        atoms[j] = 47 # Try alkene carbon instead of aromatic
                    if atoms[j] == 49: 
                        atoms[j] = 46 # Try alkene H
                angle_coeffs, found = search_angles(angle_data, atoms, 
                                                    angle_coeffs)
                if found == 0:
                    if atoms[1] == 47:
                        if atoms[0] == 13: atoms[0] = 47
                        elif atoms[2] == 13: atoms[0] = 47
                        angle_coeffs, found = search_angles(angle_data, atoms, 
                                                            angle_coeffs)
                        if found != 0:
                            print atoms
                            questionable_substitutions +=1
            if found != 1:
                 #raise ValueError('WRONG',atoms,'\t found ',found,' entries')
                 print 'WRONG',atoms,'\t found ',found,' entries'
        if questionable_substitutions != 0:
            print 'made ',questionable_substitutions,' questionable angle subs'
        return angle_coeffs

    def match_torsions(self, torsion_types):
        def add_coeff(torsion_data, torsion_coeffs, j):
            N = len(torsion_coeffs['type'])
            torsion_coeffs['type'] +=  [N+1]
            torsion_coeffs['k1']    += [ torsion_data['k1'][j] ]
            torsion_coeffs['k2']    += [ torsion_data['k2'][j] ]
            torsion_coeffs['k3']    += [ torsion_data['k3'][j] ]
            torsion_coeffs['k4']    += [ torsion_data['k4'][j] ]
            return torsion_coeffs
 
        def check_wildcards(atoms, torsion_data, torsion_coeffs):
            found = 0
            t = len(torsion_data['k1'])
            for i in range(t):
                a1 = torsion_data['a1'][i]
                a2 = torsion_data['a2'][i]
                a3 = torsion_data['a3'][i]
                a4 = torsion_data['a4'][i]
                if a1 == 0:
                    flag1 = atoms[1:4] == [a2,a3,a4]
                    flag2 = atoms[0:3] == [a4,a3,a2]
                    if flag1 or flag2:
                        found += 1
                        torsion_coeffs = add_coeff(torsion_data,torsion_coeffs,i)
                        break
                if a4 == 0:
                    flag1 = atoms[0:3] == [a1,a2,a3]
                    flag2 = atoms[1:4] == [a4,a2,a1]
                    if (flag1 or flag2):
                        foudn += 1
                        torsion_coeffs = add_coeff(torsion_data,torsion_coeffs,i)
                        break
                if a1 == 0 and a4 == 0:
                    flag1 = [atoms[1],atoms[2]] == [a2,a3]
                    flag2 = [atoms[2],atoms[1]] == [a3,a2]
                    if flag1 or flag2:
                        found += 1
                        torsion_coeffs = add_coeff(torsion_data,torsion_coeffs,i)
                        break
            return torsion_coeffs, found

        def search_torsions(torsion_data, atoms, torsion_coeffs):
            found = 0
            for j in range(len(torsion_data['k1'])):
                atom_data =[ torsion_data['a1'][j], 
                             torsion_data['a2'][j],
                             torsion_data['a3'][j],
                             torsion_data['a4'][j] ]
                flag1 = atoms==atom_data
                flag2 = atoms==list(reversed(atom_data))
                if flag1 or flag2:
                    found += 1
                    torsion_coeffs = add_coeff(torsion_data,torsion_coeffs,j)
            if found == 0:
                torsion_coeffs, found = check_wildcards(atoms,torsion_data, 
                                                        torsion_coeffs)
            return torsion_coeffs, found


        torsion_data = self.torsion_data
        torsion_coeffs = {'type':[], 'k1':[], 'k2':[], 'k3':[], 'k4':[]}
        questionable = 0
        for i in range(len(torsion_types)):
            a1 = self.type_defs[ torsion_types[i][0] ]
            a2 = self.type_defs[ torsion_types[i][1] ]
            a3 = self.type_defs[ torsion_types[i][2] ]
            a4 = self.type_defs[ torsion_types[i][3] ]
            atoms = [a1,a2,a3,a4]
            #print atoms
            
            torsion_coeffs,found = search_torsions(torsion_data, atoms, 
                                                   torsion_coeffs)
            
            if found == 0:
                for j in range(4):
                    if atoms[j] == 48: 
                        atoms[j] = 47 # Try alkene carbon instead of aromatic
                    if atoms[j] == 49: 
                        atoms[j] = 46 # Try alkene H
                torsion_coeffs,found = search_torsions(torsion_data, atoms, 
                                                   torsion_coeffs)
            if found == 0:
                change = False
                if atoms[0] == 5:
                    atoms[0] = 46
                    change = True
                if atoms[3] == 5:
                    atoms[3] = 46
                torsion_coeffs,found = search_torsions(torsion_data, atoms, 
                                                       torsion_coeffs)
                if found != 0:
                    print found, [a1,a2,a3,a4], atoms, 'made 5 -> 46 sub'

            if found == 0:
                if atoms[1:4] == [13,47,3] or atoms[0:3] == [3,47,13]:
                    atoms = [0,13,47,47]
                    torsion_coeffs,found = search_torsions(torsion_data, atoms,
                                                           torsion_coeffs)
                if found != 0:
                    print found, [a1,a2,a3,a4], atoms, 'made 13,47,47 sub'


            if found == 0:
                if atoms[1] == 13 and atoms[2] == 13:
                    if atoms[0] == 47: atoms[0] = 3
                    elif atoms[3] == 47: atoms[0] = 3
                    torsion_coeffs,found = search_torsions(torsion_data, atoms, 
                                                   torsion_coeffs)
                    if found != 0:
                        print found, [a1,a2,a3,a4], atoms, 'sub 47 -> 3'

            if found == 0:
                if atoms == [13,47,5,7] or atoms == [7,5,47,13]:
                    atoms = [7,5,47,47]
                    torsion_coeffs,found = search_torsions(torsion_data, atoms, 
                                                   torsion_coeffs)
                    if found != 0:
                        print found, [a1,a2,a3,a4], atoms, 'made phenol sub'
            
            if found != 1:
                raise ValueError('Torsion',[a1,a2,a3,a4],
                                 'found ',found,' entries')
                #print 'WRONG',[a1,a2,a3,a4],'\t found ',found,' entries'
        if questionable != 0:
            print 'made ',questionable,' questionable torsion substitutions'            
        return torsion_coeffs

    def match_impropers(self, improper_types):
        def search_impropers(improper_data, centre, neighbours, improper_coeffs):
            found = 0
            for j in range(len(improper_data['k'])):
                data_neighbours =[ improper_data['a1'][j], 
                                   improper_data['a2'][j],
                                   improper_data['a3'][j] ]
                data_centre     = improper_data['centre'][j] 
                flag1 = data_centre == centre
                # Test if other improper atom in connected to centre
                flag2 = set(neighbours) >= {data_neighbours[2]}
                flag3 = data_neighbours == [0,0,0]
                if flag1 and (flag2 or flag3):
                    found += 1
                    improper_coeffs['type'] += [i+1]
                    improper_coeffs['k']    += [ improper_data['k'][j] ]
                    improper_coeffs['r']    += [ improper_data['r'][j] ]
            return improper_coeffs, found

        improper_data = self.improper_data
        improper_coeffs = {'type':[], 'k':[], 'r':[]}
        questionable_substitutions = 0
        for i in range(len(improper_types)):
            centre = self.type_defs[ improper_types[i][0] ]
            a1 = self.type_defs[ improper_types[i][1] ]
            a2 = self.type_defs[ improper_types[i][2] ]
            a3 = self.type_defs[ improper_types[i][3] ]
            neighbours = [a1,a2,a3]
            
            improper_coeffs, found = search_impropers(improper_data, centre,
                                                      neighbours, improper_coeffs)
            if found != 1:
                 raise ValueError('Improper',centre,neighbours,
                                  'found ',found,' entries')
                 #print 'WRONG',centre,neighbours,'\t found ',found,' entries'
        return improper_coeffs 


