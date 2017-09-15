class Molecule(object):
    """A molecule or motif to be projected onto lattice points
       Structure and bonding is defined within a derivative of this class"""

    def cell_shape(self):
        raise NotImplementedError('cell_shape not defined for this motif')
    def cell_coords(self):
        raise NotImplementedError('cell_coords not defined for this motif')
    def assign_molecules(self):
        raise NotImplementedError('assing_molecules not defined for this motif')
    def assign_atom_labels(self):
        raise NotImplementedError('assign_atom_labels not defined for this motif')
    def assign_atom_charges(self):
        raise NotImplementedError('assing_atom_charges not defined for this motif')
    def assign_bonds(self):
        raise NotImplementedError('assign_bonds not defined for this motif')
    def connection_types(self):
        raise NotImplementedError('connection_types not defined for this motif')
