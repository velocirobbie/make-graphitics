import os

class Writer(object):
    def __init__(self, coords, molecule_labels):
        # Takes a numpy 3xN array of atom coordinates and outputs
        # them in different formats for viewing
        self.coords = coords
        self.molecule = molecule_labels

    def write_xyz(
            self,filename='out.xyz',system_name='comment line'):
        with open(filename,'w') as outfile:
            outfile.write(str(len(self.coords))+'\n'+system_name+'\n')
            for atom in self.coords:
                xyz=str(atom[0])+' '+str(atom[1])+' '+str(atom[2])
                outfile.write('C ' + xyz + '\n')
            print 'Coords written to '+str(filename)

    def write_lammps(
            self,system_size,filename='data.lammps',system_name='comment line'):
        # atom_type full
        with open(filename,'w') as outfile:
            outfile.write(
                    '# '+ system_name +'\n' +
                    str(len(self.coords)) +' atoms \n' 
                    '\n'
                    '1 atom types \n' 
                    '\n'
                    '0.0 \t'+str(system_size[0])+'\t xlo xhi \n'
                    '0.0 \t'+str(system_size[1])+'\t ylo yhi \n'
                    '0.0 \t'+str(system_size[2])+'\t zlo zhi \n'
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
                        str(i+1)+'\t ' +                  # atom ID
                        str(self.molecule[i]) +         # molecule ID
                        ' 1'                            # atom type
                        ' 0 '+                           # atom charge
                        str(self.coords[i][0])+' \t' +  # x
                        str(self.coords[i][1])+' \t' +  # y
                        str(self.coords[i][2])+' \n'    # z
                        )            
            print 'Coords written to '+filename


