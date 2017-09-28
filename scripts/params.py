import yaml

class Parameterise(object):
    def __init__(self, params, definitions):
        self.params = params
        self.definitions = yaml.load(open(definitions))
        self.bond_data = self.read_bond_data()
        self.angle_data = self.read_angle_data()
        self.torsion_data = self.read_torsion_data()
        self.improper_data = self.read_improper_data()

    def pair_coeffs(self, pairfile='scripts/params/opls_pair.dat'):
        pair_data = {'a':[],'s':[],'e':[]}
        with open(pairfile) as f:
            for line in f:
                line = line.split()
                pair_data['a'] += [int(line[0])]
                pair_data['e'] += [float(line[1])]
                pair_data['s'] += [float(line[2])]
        return pair_data


    def read_bond_data(self,bondfile='scripts/params/opls_bond.dat'):
        bond_data = {'a1':[], 'a2':[], 'k':[], 'r':[]}
        with open(bondfile) as f:
            for line in f:
                line = line.split()
                bond_data['a1'] += [ int(line[1]) ]
                bond_data['a2'] += [ int(line[2]) ]
                bond_data['k']  += [ float(line[3]) ]
                bond_data['r']  += [ float(line[4]) ]
        return bond_data

    def match_bonds(self, bond_types):
        bond_data = self.bond_data
        bond_coeffs = {'type':[], 'k':[], 'r':[]}
        for i in range(len(bond_types)):
            a1 = self.definitions[ bond_types[i][0] ]
            a2 = self.definitions[ bond_types[i][1] ]
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

    def read_angle_data(self,anglefile='scripts/params/opls_angle.dat'):
        angle_data = {'a1':[], 'a2':[], 'a3':[], 'k':[], 'r':[]}
        with open(anglefile) as f:
            for line in f:
                line = line.split()
                angle_data['a1'] += [ int(line[1]) ]
                angle_data['a2'] += [ int(line[2]) ]
                angle_data['a3'] += [ int(line[3]) ]
                angle_data['k']  += [ float(line[4]) ]
                angle_data['r']  += [ float(line[5]) ]
        return angle_data

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
            a1 = self.definitions[ angle_types[i][0] ]
            a2 = self.definitions[ angle_types[i][1] ]
            a3 = self.definitions[ angle_types[i][2] ]
            atoms = [a1,a2,a3]
            
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

    def read_torsion_data(self,torsionfile='scripts/params/opls_torsion.dat'):
        torsion_data = {'a1':[], 'a2':[], 'a3':[], 'a4':[],
                        'k1':[], 'k2':[], 'k3':[], 'k4':[]}
        with open(torsionfile) as f:
            for line in f:
                line = line.split()
                torsion_data['a1'] += [ int(line[1]) ]
                torsion_data['a2'] += [ int(line[2]) ]
                torsion_data['a3'] += [ int(line[3]) ]
                torsion_data['a4'] += [ int(line[4]) ]
                torsion_data['k1']  += [ float(line[5]) ]
                torsion_data['k2']  += [ float(line[8]) ]
                torsion_data['k3']  += [ float(line[11]) ]
                if len(line) == 17:
                    torsion_data['k4']  += [ float(line[14]) ]
                else:
                    torsion_data['k4']  += [ 0.0 ]
        return torsion_data

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
            a1 = self.definitions[ torsion_types[i][0] ]
            a2 = self.definitions[ torsion_types[i][1] ]
            a3 = self.definitions[ torsion_types[i][2] ]
            a4 = self.definitions[ torsion_types[i][3] ]
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
                #raise ValueError('WRONG',atoms,'\t found ',found,' entries')
                print 'WRONG',[a1,a2,a3,a4],'\t found ',found,' entries'
        if questionable != 0:
            print 'made ',questionable,' questionable torsion substitutions'            
        return torsion_coeffs

    def match_impropers(self, improper_types):
        def search_impropers(improper_data, atoms, improper_coeffs):
            found = 0
            for j in range(len(improper_data['k'])):
                atom_data =[ improper_data['centre'][j], 
                             improper_data['a1'][j],
                             improper_data['a2'][j],
                             improper_data['a3'][j] ]
                flag1 = atom_data[0] == atoms[0]
                # Test if other improper atom in connected to centre
                flag2 = set(atoms[1:]) >= {atom_data[3]}
                flag3 = atom_data[1:] == [0,0,0]
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
            centre = self.definitions[ improper_types[i][0] ]
            a1 = self.definitions[ improper_types[i][1] ]
            a2 = self.definitions[ improper_types[i][2] ]
            a3 = self.definitions[ improper_types[i][3] ]
            atoms = [centre,a1,a2,a3]
            
            improper_coeffs, found = search_impropers(improper_data, 
                                                      atoms, improper_coeffs)
            if found != 1:
                 #raise ValueError('WRONG',atoms,'\t found ',found,' entries')
                 print 'WRONG',atoms,'\t found ',found,' entries'
        return improper_coeffs 

    def read_improper_data(self,improperfile='scripts/params/opls_improper.dat'):
        improper_data = {'a1':[], 'a2':[], 'centre':[], 'a3':[], 'k':[], 'r':[]}
        with open(improperfile) as f:
            for line in f:
                line = line.split()
                improper_data['a1'] += [ int(line[1]) ]
                improper_data['a2'] += [ int(line[2]) ]
                improper_data['centre'] += [ int(line[3]) ]
                improper_data['a3'] += [ int(line[4]) ]
                improper_data['k']  += [ float(line[5]) ]
                improper_data['r']  += [ float(line[6]) ]
        return improper_data


