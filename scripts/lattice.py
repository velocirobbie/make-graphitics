import numpy as np

class Lattice(object):
    def __init__(self,cell_dimensions):
        self.cell_dimensions = cell_dimensions

    def lattice_size_vdw(self,vdw_cutoff):
        # Will result in a cell at least as big as 2 * vdw_cutoff
        # in all directions
        Nlattice_points = []
        for axis in self.cell_dimensions:
            Npoints_on_axis = int(np.ceil((2*vdw_cutoff)/axis))
            Nlattice_points.append(Npoints_on_axis)
        return Nlattice_points

    def lattice_size_layers(self,vdw_cutoff,N_layers):
        # Define a cell as big as 2 * vdw_cutoff in x y direction
        # and a specified number of cells in z
        Nlattice_points= []
        Nlattice_points.append(int(np.ceil(2*vdw_cutoff/self.cell_dimensions[0])))
        Nlattice_points.append(int(np.ceil(2*vdw_cutoff/self.cell_dimensions[1])))
        Nlattice_points.append(N_layers)
        return Nlattice_points

    def create_lattice_points(self,Nlattice_points):
        a,b,c = self.cell_dimensions
        lattice_points = [] 
        for x in range(Nlattice_points[0]):
            for y in range(Nlattice_points[1]):
                for z in range(Nlattice_points[2]):
                     point = [x*a,y*b,z*c]
                     lattice_points.append(point)
        lattice_points = np.array(lattice_points)
        return lattice_points

    def cell_onto_lattice(self,cell_coords,lattice_points):
        N_atoms = len(cell_coords) * len(lattice_points)
        atoms = np.empty([N_atoms,3],float)
        #counter = 0
        for i,lattice_point in enumerate(lattice_points):
            start = i * N_atoms
            end = (i+1 * N_atoms)
            atoms[start:end] = np.array(cell_coords) + np.array(lattice_point)
            #for atom in cell_coords:
            #    new_atom = np.array(atom) + np.array(lattice_point)
            #    atoms[counter] = new_atom
            #    counter += 1
        return atoms

    def system_size(self,Nlattice_points):
        box_vectors = np.array(Nlattice_points)*np.array(self.cell_dimensions)
        box_dimensions = np.vstack((np.zeros(3),box_vectors)).transpose()
        return box_dimensions
