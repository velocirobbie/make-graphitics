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
        self.attributes = {'masses', 'coords', 'molecules', 
                      'charges', 'types',
                      'xlo', 'xhi', 'ylo', 'yhi', 'zlo', 'zhi',
                      'bonds', 'bond_types','Nbond_types',
                      'angles', 'angle_types','Nangle_types',
                      'dihedrals', 'dihedral_types',
                      'Ndihedral_types',
                      'impropers', 'improper_types',
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
        self.masses = np.zeros(self.Natom_types)
        for i in range(self.Natom_types):
            line = self.read(datafile)
            index = int(line[0]) - 1
            atom_type = int(index)
            atom_mass = float(line[1])
            self.masses[atom_type] = atom_mass

    def read_atoms(self,datafile):
        self.coords = np.zeros((self.Natoms,3))
        self.charges = np.zeros(self.Natoms)
        self.molecules = np.zeros(self.Natoms)
        self.types = np.zeros(self.Natoms)
        for i in range(self.Natoms):
            line = self.read(datafile)
            index = int(line[0]) - 1
            self.coords[index] = line[4:7]
            self.molecules[index] = line[1]
            self.types[index] = line[2]
            self.charges[index] = line[3]

    def read_velocities(self,datafile):
        print '--- Ignoring Velocities --- noone cares'
        for i in range(self.Natoms):
            line = self.read(datafile)

    def read_data_line(self,datafile,types,atoms):
        line = self.read(datafile)
        index = int(line[0]) - 1 # index line up with numpy array
        types[index] = line[1]
        atoms[index] = line[2:]

    def read_bonds(self,datafile):
        self.bonds = np.zeros((self.Nbonds,2))
        self.bond_types = np.zeros((self.Nbonds))
        for i in range(self.Nbonds):
            self.read_data_line(datafile,
                                self.bond_types,self.bonds)

    def read_angles(self,datafile):
        self.angles = np.zeros((self.Nangles,3))
        self.angle_types = np.zeros(self.Nangles)
        for i in range(self.Nangles):
            self.read_data_line(datafile,
                                self.angle_types,self.angles)

    def read_dihedrals(self,datafile):
        self.dihedrals = np.zeros((self.Ndihedrals,4))
        self.dihedral_types = np.zeros(self.Ndihedrals)
        for i in range(self.Ndihedrals):
            self.read_data_line(datafile,
                                self.dihedral_types,self.dihedrals)

    def read_impropers(self,datafile):
        self.impropers = np.zeros((self.Nimpropers,4))
        self.improper_types = np.zeros(self.Nimpropers)
        for i in range(self.Nimpropers):
            self.read_data_line(datafile,
                                self.improper_types,self.impropers)

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
            
