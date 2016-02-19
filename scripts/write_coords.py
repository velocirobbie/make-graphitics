import os

class Writer(object):
    def __init__(self, coords, molecule_labels=[],
            bonds=[],angles=[],torsions=[],box_dimensions=[],
            system_name='comment line'):
        # Takes a numpy 3xN array of atom coordinates and outputs
        # them in different formats for viewing/modelling
        self.coords = coords
        self.molecule = molecule_labels
        self.bonds = bonds
        self.angles = angles
        self.torsions = torsions
        self.size = box_dimensions
        self.system_name = system_name
        
    def write_xyz(
            self,filename='out.xyz'):
        with open(filename,'w') as outfile:
            outfile.write(str(len(self.coords))+'\n'
                    +self.system_name+'\n')
            for atom in self.coords:
                xyz=str(atom[0])+' '+str(atom[1])+' '+str(atom[2])
                outfile.write('C ' + xyz + '\n')
            print 'Coords written to '+str(filename)

    def write_lammps(
            self,filename='data.lammps'):
        # atom_type full
        with open(filename,'w') as outfile:
            outfile.write(
                    '# '+ self.system_name +'\n' +
                    str(len(self.coords)) +' atoms \n')
            if len(self.bonds): outfile.write(
                    str(len(self.bonds)) +' bonds \n')
            if len(self.angles): outfile.write(
                    str(len(self.angles)) +' angles \n')
            if len(self.torsions): outfile.write(
                    str(len(self.torsions))+' dihedrals \n')
            outfile.write('\n'
                    '1 atom types \n')
            if len(self.bonds): outfile.write('1 bond types \n')
            if len(self.angles): outfile.write('1 angle types \n')
            if len(self.torsions): outfile.write(
                    '1 dihedral types \n')
            outfile.write('\n'
                    '0.0 \t'+str(self.size[0])+'\t xlo xhi \n'
                    '0.0 \t'+str(self.size[1])+'\t ylo yhi \n'
                    '0.0 \t'+str(self.size[2])+'\t zlo zhi \n'
                    '\n'
                    'Masses \n'
                    '\n'
                    '1 12.01 \n'
                    '\n'
                    'Atoms \n'
                    '\n'
                    )
            
            for i in range(len(self.coords)):
                outfile.write(
                        str(i+1)+'\t ' +              # atom ID
                        str(self.molecule[i]) +       # molecule ID
                        ' 1'                          # atom type
                        ' 0 '+                        # atom charge
                        str(self.coords[i][0])+' \t' +# x
                        str(self.coords[i][1])+' \t' +# y
                        str(self.coords[i][2])+' \n'  # z
                        )            
            
            if len(self.bonds):
                outfile.write('\n Bonds \n \n')
                for i in range(len(self.bonds)):
                    outfile.write(
                            str(i+1)+'\t' +           # bond ID
                            ' 1 '  +                    # bond type
                            str(self.bonds[i][0])+' \t'+# atom 1
                            str(self.bonds[i][1])+' \n'# atom 2
                            )
             
            if len(self.angles):
                outfile.write('\n Angles \n \n')
                for i in range(len(self.angles)):
                    outfile.write(
                            str(i+1)+'\t' +           # angle ID
                            ' 1 ' +                    # angle type
                            str(self.angles[i][0])+' \t'+
                            str(self.angles[i][1])+' \t'+
                            str(self.angles[i][2])+' \n'
                            )
            
            if len(self.torsions):
                outfile.write('\n Dihedrals \n \n')
                for i in range(len(self.torsions)):
                    outfile.write(
                            str(i+1)+'\t'+
                            ' 1 '+
                            str(self.torsions[i][0])+' \t'+
                            str(self.torsions[i][1])+' \t'+
                            str(self.torsions[i][2])+' \t'+
                            str(self.torsions[i][3])+' \n'
                            )
                        

                

            print 'Coords written to '+filename


