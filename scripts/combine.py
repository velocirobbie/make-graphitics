import numpy as np

class Combine(object):
    def __init__(self,sim1,sim2):
        # combine two Crystal objects
        # keep cell size from sim1
        natoms1 = len(sim1.coords)
        nmols1 = np.amax(sim1.molecule_labels)
        
        self.box_dimensions = sim1.box_dimensions
        self.angle_labels = self.join(sim1.angle_labels,
                                      sim2.angle_labels)
        self.angles = self.stack(sim1.angles,sim2.angles+natoms1)
        self.atom_charges = self.join(sim1.atom_charges,
                                      sim2.atom_charges)
        self.atom_labels = self.join(sim1.atom_labels,
                                     sim2.atom_labels)
        self.bond_labels = self.join(sim1.bond_labels,
                                     sim2.bond_labels)
        self.bonds = self.stack(sim1.bonds,sim2.bonds+natoms1)
        self.coords = self.stack(sim1.coords,sim2.coords)
        self.molecule_labels=self.join(
                sim1.molecule_labels,
                list(np.array(sim2.molecule_labels)+nmols1))
        self.torsion_labels = self.join(sim1.torsion_labels,
                                        sim2.torsion_labels)
        self.torsions = self.stack(sim1.torsions,
                                   sim2.torsions+natoms1)
        self.improper_labels = self.join(sim1.improper_labels,
                                         sim2.improper_labels)
        self.impropers = self.stack(sim1.impropers,
                                    sim2.impropers+natoms1)

    def stack(self,array1,array2):
        return np.vstack((array1,array2))

    def join(self,list1,list2):
        return list1 + list2
