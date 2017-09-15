import os
import numpy as np

class Writer(object):
    def __init__(self, sim,
            system_name='comment line'):
        # Takes a numpy 3xN array of atom coordinates and outputs
        # them in different formats for viewing/modelling
        self.coords = sim.coords
        self.atom_labels = sim.atom_labels
        self.natom_types = len(np.unique(self.atom_labels))
        self.molecule = sim.molecule_labels
        self.charges = sim.atom_charges
        self.bonds = sim.bonds
        self.bond_labels = sim.bond_labels
        self.nbond_types = len(np.unique(self.bond_labels))
        self.angles = sim.angles
        self.angle_labels = sim.angle_labels
        self.nangle_types = len(np.unique(self.angle_labels))
        self.torsions = sim.torsions
        self.torsion_labels = sim.torsion_labels
        self.ntorsion_types = len(np.unique(self.torsion_labels))
        self.impropers = sim.impropers
        self.improper_labels = sim.improper_labels
        self.nimproper_types = len(np.unique(self.improper_labels))
        self.size = sim.box_dimensions
        self.system_name = system_name
        
        self.atom_masses = []
        for atom in np.unique(self.atom_labels):
            if atom == 1: self.atom_masses.append(12.01)
            if atom == 2: self.atom_masses.append(1.00)
            if atom == 3: self.atom_masses.append(12.01)
            if atom == 4: self.atom_masses.append(16)
            if atom == 5: self.atom_masses.append(1)
            if atom == 6: self.atom_masses.append(12)
            if atom == 7: self.atom_masses.append(16)


    def write_xyz(self,filename='out.xyz'):
        with open(filename,'w') as outfile:
            outfile.write(str(len(self.coords))+'\n'
                    +self.system_name+'\n')
            for i in range(len(self.coords)):
                xyz=(str(self.coords[i][0])+' '+
                     str(self.coords[i][1])+' '+
                     str(self.coords[i][2]))
                if self.atom_labels[i] == 1: atom_label = 'C '
                elif self.atom_labels[i] == 2: atom_label = 'H '
                elif self.atom_labels[i] == 3: atom_label = 'C '
                elif self.atom_labels[i] == 5: atom_label = 'H '
                elif self.atom_labels[i] == 4: atom_label = 'O '
                elif self.atom_labels[i] == 6: atom_label = 'C '
                elif self.atom_labels[i] == 7: atom_label = 'O '
                else: atom_label = str(self.atom_labels[i])+' '
                outfile.write(atom_label + xyz + '\n')
            print 'Coords written to '+str(filename)

    def write_lammps(
            self,filename='data.lammps'):
        # atom_type full
        with open(filename,'w') as outfile:
            outfile.write(
                    '# '+ self.system_name +'\n' +
                    str(len(self.coords)) +' atoms \n'+
                    str(len(self.bonds)) +' bonds \n'+
                    str(len(self.angles)) +' angles \n'+
                    str(len(self.torsions))+' dihedrals \n'+
                    str(len(self.impropers))+' impropers \n'
                    '\n'+
                    str(self.natom_types)+' atom types \n'+
                    str(self.nbond_types)+' bond types \n'+
                    str(self.nangle_types)+' angle types \n'+
                    str(self.ntorsion_types)+' dihedral types \n'+
                    str(self.nimproper_types)+' improper types \n'+
                    '\n'
                    '0.0 \t'+str(self.size[0])+'\t xlo xhi \n'
                    '0.0 \t'+str(self.size[1])+'\t ylo yhi \n'
                    '0.0 \t'+str(self.size[2])+'\t zlo zhi \n'
                    '\n'
                    'Masses \n'
                    '\n')
            for i in range(self.natom_types):
                outfile.write(
                      str(i+1)+'\t '+
                      str(self.atom_masses[i])+'\n')
            outfile.write('\n Atoms \n \n')
            
            for i in range(len(self.coords)):
                outfile.write(
                        str(i+1)+'\t ' +             # atom ID
                        str(self.molecule[i])+'\t '+ # molecule ID
                        str(self.atom_labels[i])+'\t '+#atom type
                        str(self.charges[i])+'\t '+#atomcharg
                        str(self.coords[i][0])+'\t ' +# x
                        str(self.coords[i][1])+'\t ' +# y
                        str(self.coords[i][2])+'\n '  # z
                        )            
            
            if len(self.bonds):
                outfile.write('\n Bonds \n \n')
                for i in range(len(self.bonds)):
                    outfile.write(
                            str(i+1)+'\t ' +           # bond ID
                            str(self.bond_labels[i])+'\t '+
                            str(self.bonds[i][0])+'\t '+# atom 1
                            str(self.bonds[i][1])+'\n'# atom 2
                            )
             
            if len(self.angles):
                outfile.write('\n Angles \n \n')
                for i in range(len(self.angles)):
                    outfile.write(
                            str(i+1)+'\t ' +         # angle ID
                            str(self.angle_labels[i])+'\t '+
                            str(self.angles[i][0])+'\t '+
                            str(self.angles[i][1])+'\t '+
                            str(self.angles[i][2])+'\n'
                            )
            
            if len(self.torsions):
                outfile.write('\n Dihedrals \n \n')
                for i in range(len(self.torsions)):
                    outfile.write(
                            str(i+1)+'\t'+
                            str(self.torsion_labels[i])+' \t'+
                            str(self.torsions[i][0])+' \t'+
                            str(self.torsions[i][1])+' \t'+
                            str(self.torsions[i][2])+' \t'+
                            str(self.torsions[i][3])+' \n'
                            )
            
            if len(self.impropers):
                outfile.write('\n Impropers \n \n')
                for i in range(len(self.impropers)):
                    outfile.write(
                            str(i+1)+'\t'+
                            str(self.improper_labels[i])+'\t'+
                            str(self.impropers[i][0])+' \t'+
                            str(self.impropers[i][1])+' \t'+
                            str(self.impropers[i][2])+' \t'+
                            str(self.impropers[i][3])+' \n'
                            )

            print 'Coords written to '+filename


