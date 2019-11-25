from sim import Sim
import numpy as np
import copy

class Combine(Sim):
    def __init__(self,sim1,sim2):
        # combine two Crystal objects
        # keep cell size from sim1
        natoms1 = len(sim1.coords)
        nmols1 = np.amax(sim1.molecule_labels)
        
        self.vdw_defs = copy.deepcopy(sim1.vdw_defs) 
        #self.pair_coeffs = sim1.pair_coeffs
        #self.masses = sim1.masses
        for i in sim2.vdw_defs:
            exists_in_sim1 = 0
            for j in sim1.vdw_defs:
                if sim1.vdw_defs[j] == sim2.vdw_defs[i]:
                    exists_in_sim1 += 1
                    sim2.atom_labels = self.replace_labels(sim2.atom_labels,i,j)
            if not exists_in_sim1:
                new_label = max(self.vdw_defs.keys()) + 1
        #        self.pair_coeffs[new_label] = sim2.pair_coeffs[i]
        #        self.masses[new_label] = sim2.masses[i]
                sim2.atom_labels = self.replace_labels(sim2.atom_labels,i,new_label)
                self.vdw_defs[new_label] = sim2.vdw_defs[i]

            elif exists_in_sim1 > 1:
                raise Exception(exists_in_sim1)
        
        self.box_dimensions = sim1.box_dimensions
 
        #for thing in ['bond','angle','improper']:
            #coeff = thing+'_coeffs'
            #types = thing+'_types'
            #labels= thing+'_labels'
#            N     = 'N'+thing+'_types'
        #    self.combine_coeff(sim1, sim2, coeff,types,labels)
        #self.combine_coeff(sim1, sim2,'dihedral_coeffs','dihedral_types','dihedral_labels')

        setattr(self,'coords',self.stack(getattr(sim1,'coords'),
                                         getattr(sim2,'coords')))
        for attr in ['bonds','angles','dihedrals','impropers']:
            attr1 = getattr(sim1,attr)
            attr2 = getattr(sim2,attr) + len(getattr(sim1,'coords'))
            #print attr,len(attr1),len(attr2)
            setattr(self,attr,self.stack(attr1,attr2))
            #print len(getattr(self,attr))
        
        for attr in ['atom_charges','atom_labels','bond_labels','angle_labels',
                     'dihedral_labels','improper_labels']:
            attr1 = getattr(sim1,attr)
            attr2 = getattr(sim2,attr)
            setattr(self,attr,self.join(attr1,attr2))

        self.molecule_labels=self.join(
                sim1.molecule_labels,
                list(np.array(sim2.molecule_labels)+nmols1))
 
        
    def combine_coeff(self, sim1, sim2, coeff, types, labels):
        new_types  = copy.deepcopy(getattr(sim1, types))
        new_coeffs = copy.deepcopy(getattr(sim1, coeff))
        types1 = getattr(sim1,types)
        types2 = getattr(sim2,types)
        for i in range(len(types2)):
            def2  = [sim2.vdw_defs[a] for a in types2[i]]
            exists_in_sim1 = 0
            for j in range(len(types1)):
                def1 = [sim1.vdw_defs[a] for a in types1[j]]
                if (def1 == def2) or (def1 == list(reversed(def2))):
                    print 'matched',coeff,i+1,j+1,def1,def2
                    exists_in_sim1 += 1
                    new_labels = self.replace_labels(getattr(sim2,labels),i+1,j+1)
                    setattr(sim2,labels,new_labels)
            if not exists_in_sim1:
                new_label = max(new_coeffs.keys()) + 1
                new_coeff = getattr(sim2,coeff)[i+1]
                new_coeffs[new_label] = new_coeff
                new_labels = self.replace_labels(getattr(sim2,labels),i+1,new_label)
                setattr(sim2,labels,new_labels)
                new_types += [ getattr(sim2,types)[i] ]
            elif exists_in_sim1 > 1:
                raise IndexError(exists_in_sim1)
        setattr(self,coeff,new_coeffs)
        setattr(self,types,new_types)

    def replace_labels(self,labels,a,b):
        count = 0
        for label in range(len(labels)):
            if labels[label] == a:
                labels[label] = b
                count +=1
        return labels

    def stack(self,array1,array2):
        return np.vstack((array1,array2))

    def join(self,list1,list2):
        return list1 + list2

    def reduce(self, sim1, sim2, labels, coeffs):
        pass
