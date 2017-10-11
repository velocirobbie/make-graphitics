import numpy as np

class Combine(object):
    def __init__(self,sim1,sim2):
        # combine two Crystal objects
        # keep cell size from sim1
        natoms1 = len(sim1.coords)
        nmols1 = np.amax(sim1.molecule_labels)
        
        self.box_dimensions = sim1.box_dimensions
        self.atom_charges = self.join(sim1.atom_charges,
                                      sim2.atom_charges)
        
        for attr in ['coords','bonds','angles','torsions','impropers']:
            attr1 = getattr(sim1,attr)
            attr2 = getattr(sim2,attr)
            setattr(self,attr,self.stack(attr1,attr2))
        
        for attr in ['atom_labels','bond_labels','angle_labels',
                     'torsion_labels','improper_labels']:
            attr1 = getattr(sim1,attr)
            attr2 = getattr(sim2,attr)
            setattr(self,attr,self.join(attr1,attr2))

        self.molecule_labels=self.join(
                sim1.molecule_labels,
                list(np.array(sim2.molecule_labels)+nmols1))
        
        for coeff in ['pair_coeffs','bond_coeffs','angle_coeffs',
                      'torsion_coeffs','improper_coeffs','masses']:
            if hasattr(sim1,coeff):
                setattr(self,coeff,getattr(sim1,coeff))
            print coeff
            print getattr(self,coeff)

    def stack(self,array1,array2):
        return np.vstack((array1,array2))

    def join(self,list1,list2):
        return list1 + list2

    def reduce(self, sim1, sim2, labels, coeffs):
        pass
