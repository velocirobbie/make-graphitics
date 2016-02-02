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
            outfile.write('# '+ system_name +'\n \n')
            outfile.write(str(len(self.coords))+' atoms \n \n')
            outfile.write('1 atom types \n \n')
            outfile.write('0.0'+str(system_size[0])+' xlo xhi \n')
            outfile.write('0.0'+str(system_size[1])+' ylo yhi \n')
            outfile.write('0.0'+str(system_size[2])+' zlo zhi \n \n')
            outfile.write('Pair Coeffs \n \n')
            outfile.write('1   0.1 3.56    0.1     3.56 \n \n')
            outfile.write('Atoms \n \n')
            i=1
            for atom in self.coords:
                outfile.write(str(i)+' 1 1 '+str(atom[0])+str(atom[1])+str(atom[2])+' \n')
                i+=1
            print 'Coords written to '+str(filename)
