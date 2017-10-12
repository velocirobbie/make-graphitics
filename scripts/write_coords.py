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

        if hasattr(sim,'masses'): 
            self.masses = sim.masses
        if hasattr(sim,'bond_coeffs'): 
            self.bond_coeffs = sim.bond_coeffs
        if hasattr(sim,'angle_coeffs'): 
            self.angle_coeffs = sim.angle_coeffs
        if hasattr(sim,'dihedral_coeffs'): 
            self.dihedral_coeffs = sim.dihedral_coeffs
        if hasattr(sim,'improper_coeffs'): 
            self.improper_coeffs = sim.improper_coeffs
        if hasattr(sim,'pair_coeffs'): 
            self.pair_coeffs = sim.pair_coeffs


    def write_xyz(self,filename='out.xyz'):
        with open(filename,'w') as outfile:
            outfile.write(str(len(self.coords))+'\n'
                    +self.system_name+'\n')
            for i in range(len(self.coords)):
                xyz=(str(self.coords[i][0])+' '+
                     str(self.coords[i][1])+' '+
                     str(self.coords[i][2]))
                if   self.atom_labels[i] in [1,8,11]: 
                    atom_label = 'C '
                elif self.atom_labels[i] in [2,5]: 
                    atom_label = 'H '
                elif self.atom_labels[i] in [3]: 
                    atom_label = 'N '
                elif self.atom_labels[i] in [4,6,7,9,10]: 
                    atom_label = 'O '
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
                    '\n')
            if hasattr(self,'masses'):
                outfile.write('\n Masses \n \n')
                for i in self.masses:
                    outfile.write(
                            str(i)+'\t'+
                            str(self.masses[i])+'\n'
                               )
            
            if hasattr(self,'pair_coeffs'):    
                outfile.write('\n Pair Coeffs \n \n')
                for i in self.pair_coeffs:
                    outfile.write(
                            str(i)+'\t'+
                            str(self.pair_coeffs[i][1])+'\t'+
                            str(self.pair_coeffs[i][2])+'\n'
                            )
 
            if hasattr(self,'bond_coeffs'):
                outfile.write('\n Bond Coeffs \n \n')
                for i in self.bond_coeffs:
                    outfile.write(
                            str(i)+'\t'+
                            str(self.bond_coeffs[i][1])+'\t'+
                            str(self.bond_coeffs[i][2])+'\n'
                            )
            if hasattr(self,'angle_coeffs'):
                outfile.write('\n Angle Coeffs \n \n')
                for i in self.angle_coeffs:
                    outfile.write(
                            str(i)+'\t'+
                            str(self.angle_coeffs[i][1])+'\t'+
                            str(self.angle_coeffs[i][2])+'\n'
                            )
            if hasattr(self,'dihedral_coeffs'):    
                outfile.write('\n Dihedral Coeffs \n \n')
                for i in self.dihedral_coeffs:
                    outfile.write(
                            str(i)+'\t'+
                            str(self.dihedral_coeffs[i][1])+'\t'+
                            str(self.dihedral_coeffs[i][2])+'\t'+
                            str(self.dihedral_coeffs[i][3])+'\t'+
                            str(self.dihedral_coeffs[i][4])+'\n'
                            )

            if hasattr(self,'improper_coeffs'):    
                outfile.write('\n Improper Coeffs \n \n')
                for i in self.improper_coeffs:
                    outfile.write(
                            str(i)+'\t'+
                            str(self.improper_coeffs[i][1])+'\t'+
                            str(self.improper_coeffs[i][2])+'\n'
                            )



                
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


