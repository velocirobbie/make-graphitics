import numpy as np

class ReadLammpsData(object):
    def __init__(self,filename):
        """ Given a lammps data file this script will read in all 
        the available data. 
        Data that can be read in are listed in self.attributes.
        Coefficients will not be read and should be in a separate
        input file for lammps.
        Copy attributes to another class with:
        ReadLammpsData.__dict__ = obj.__dict__.copy()
        """
        self.filename = filename
        self.attributes = {'masses', 'coords', 'molecule_labels', 
                      'atom_charges', 'atom_labels','atom_ids',
                      'box_dimensions',
                      'bonds', 'bond_labels','Nbond_types',
                      'angles', 'angle_labels','Nangle_types',
                      'dihedrals', 'dihedral_labels',
                      'Ndihedral_types',
                      'impropers', 'improper_labels',
                      'Nimproper_types',
                      'pair_coeffs','bond_coeffs','angle_coeffs',
                      'dihedral_coeffs','improper_coeffs'}

        with open(self.filename) as datafile:
            number_of_lines = 0
            for line in datafile:
                number_of_lines += 1
            datafile.seek(0)
        
            self.count = 0
            while self.count < number_of_lines:
                line = self.read(datafile)
                self.analyse(line,datafile)
        self.box_dimensions = np.array([[self.xlo,self.xhi],
                                        [self.ylo,self.yhi],
                                        [self.zlo,self.zhi]])
        self.validate()
        extra_attributes = set(self.__dict__.keys())-self.attributes
        for attribute in extra_attributes:
            delattr(self,attribute)


    def read(self,datafile):
        self.count += 1
        return datafile.readline().split()

    def analyse(self,line,datafile):
        def could_not_read(unknown):
            print 'Could not decipher: ',unknown,' on line ',self.count
        
        def is_number(s):
            try:
                float(s)
                return True
            except ValueError:
                return False
        Ncoeff = {
                'Pair'      : 'Natom_types',
                'Bond'      : 'Nbond_types',
                'Angle'     : 'Nangle_types',
                'Dihedral'  : 'Ndihedral_types',
                'Improper'  : 'Nimproper_types'
                }

        main = {
                'Masses'    : self.read_masses,
                'Atoms'     : self.read_atoms,
                'Velocities': self.read_velocities,
                'Bonds'     : self.read_bonds,
                'Angles'    : self.read_angles,
                'Dihedrals' : self.read_dihedrals,
                'Impropers' : self.read_impropers
                }
        
        header_numbers = {
                'atoms'     : 'Natoms',
                'bonds'     : 'Nbonds',
                'angles'    : 'Nangles',
                'dihedrals' : 'Ndihedrals',
                'impropers' : 'Nimpropers'
                }

        header_types = {
                'atom'     : 'Natom_types',
                'bond'     : 'Nbond_types',
                'angle'    : 'Nangle_types',
                'dihedral' : 'Ndihedral_types',
                'improper' : 'Nimproper_types'
                }

        box = [['xlo','xhi'],['ylo','yhi'],['zlo','zhi']]

        if not line: return # blank line
        l = line[0]
        if l == '#': return # catch comment line
        elif l[0] == '#': return # catch comment line
        if is_number(l):
            if len(line) == 2:
                attribute = header_numbers.get(line[1],False)
                if not attribute: could_not_read(line)
                else: setattr(self,attribute,int(l))
            
            if len(line) ==3:
                attribute = header_types.get(line[1],False)
                if not attribute: could_not_read(line)
                else: setattr(self,attribute,int(l))

            if len(line) == 4:
                if line[2:4] in box:
                    setattr(self,line[2],float(line[0]))
                    setattr(self,line[3],float(line[1]))
                else: could_not_read(line)
        elif len(line) == 2:
            print 'reading '+l+' Coeffs'
            self.read(datafile)
            name = line[0].lower()+'_'+line[1].lower()
            N = getattr(self,Ncoeff.get(l))
            self.read_coeffs(datafile,name,N)

        else: 
            func = main.get(l,False)
            if not func: could_not_read(line)
            else: 
                print 'reading '+l
                self.read(datafile)
                func(datafile)

    def read_masses(self,datafile):
        self.masses = {}
        for i in range(self.Natom_types):
            line = self.read(datafile)
            self.masses[int(line[0])] = float(line[1])

    def read_atoms(self,datafile):
        self.atom_ids = np.empty(self.Natoms,dtype=int)
        self.coords = np.zeros((self.Natoms,3))
        self.atom_charges = np.zeros(self.Natoms)
        self.molecule_labels = np.zeros(self.Natoms,dtype=int)
        self.atom_labels = np.zeros(self.Natoms,dtype=int)
        for i in range(self.Natoms):
            line = self.read(datafile)
            self.atom_ids[i] = line[0]
            self.coords[i] = line[4:7]
            self.molecule_labels[i] = line[1]
            self.atom_labels[i] = line[2]
            self.atom_charges[i] = line[3]

    def read_velocities(self,datafile):
        print '--- Ignoring Velocities --- noone cares'
        for i in range(self.Natoms):
            line = self.read(datafile)

    def read_data_line(self,datafile,atom_labels,atoms):
        line = self.read(datafile)
        index = int(line[0]) - 1 # index line up with numpy array
        atom_labels[index] = int(line[1])
        atoms[index] = [int(atom) for atom in line[2:]]

    def read_bonds(self,datafile):
        self.bonds = np.zeros((self.Nbonds,2),int)
        self.bond_labels = np.zeros((self.Nbonds),dtype=int)
        for i in range(self.Nbonds):
            self.read_data_line(datafile,
                                self.bond_labels,self.bonds)

    def read_angles(self,datafile):
        self.angles = np.zeros((self.Nangles,3),int)
        self.angle_labels = np.zeros(self.Nangles,dtype=int)
        for i in range(self.Nangles):
            self.read_data_line(datafile,
                                self.angle_labels,self.angles)

    def read_dihedrals(self,datafile):
        self.dihedrals = np.zeros((self.Ndihedrals,4),int)
        self.dihedral_labels = np.zeros(self.Ndihedrals,dtype=int)
        for i in range(self.Ndihedrals):
            self.read_data_line(datafile,
                                self.dihedral_labels,self.dihedrals)

    def read_impropers(self,datafile):
        self.impropers = np.zeros((self.Nimpropers,4),int)
        self.improper_labels = np.zeros(self.Nimpropers,dtype=int)
        for i in range(self.Nimpropers):
            self.read_data_line(datafile,
                                self.improper_labels,self.impropers)

    def read_coeffs(self,datafile,coeffs,N):
        setattr(self,coeffs,{})
        coeffdict = getattr(self,coeffs)
        for i in range(N):
            line = self.read(datafile)
            label = int(line[0])
            coeffdict[label] = {}
            for k in range(1, len(line)):
                coeffdict[label][k] = float(line[k])

    def validate(self):
       for attribute in self.attributes:
            try:
                a = getattr(self,attribute)
            except AttributeError:
                print 'WARNING: undefined ' + attribute
            
