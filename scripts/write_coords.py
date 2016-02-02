import os

class Writer(object):
    def __init__(self, coords):
        # Takes a numpy 3xN array of atom coordinates and outputs
        # them in different formats for viewing
        self.coords = coords

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
        with open(filename,'w') as outfile:
            outfile.write(
                    '# '+ system_name +
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
                    '/n'
                    'Pair Coeffs \n' 
                    '\n'
                    '1  0.1     3.56 \n' 
                    '\n'
                    'Atoms \n'
                    '\n'
                    )
            
            i=1
            for atom in self.coords:
                outfile.write(str(i)+'\t 1 1 0 '+str(atom[0])+'\t '+str(atom[1])+'\t '+str(atom[2])+' \n')
                i+=1
            
            print 'Coords written to '+str(filename)


