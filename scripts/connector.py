import numpy as np

class Connector(object):
    def __init__(self):
        pass

    def index_cell(self,cell_position,lattice_dimensions,
                   atoms_per_cell):
        N = atoms_per_cell
        total = cell_position[2] * N
        total += cell_position[1] * lattice_dimensions[2] * N
        total += (cell_position[0] * lattice_dimensions[2] 
                * lattice_dimensions[1] * N)
        return total

    def graphite_bonds(self,lattice_dimensions):
        internal_bonds = np.array([[1,2],[2,3],[3,4],
                                   [5,6],[6,7],[8,5]])
        bonds = np.empty((0,2),dtype=int)
        # loop through all cells
        for x in range(lattice_dimensions[0]):
         for y in range(lattice_dimensions[1]):
          for z in range(lattice_dimensions[2]):
            cell_position = [x,y,z]
            xcell_position,ycell_position,xycell_position = (
             find_adjacent_cells(cell_position,lattice_dimensions))
            
            # Add bonds within cell
            i = self.index_cell([x,y,z],lattice_dimensions,8)
            bonds = np.vstack((bonds,internal_bonds+i))
            # Add bonds that cross cell boundries
            bonds = self.add_cross_bond(lattice_dimensions,
                    xcell_position,[3,2],bonds,i)
            bonds = self.add_cross_bond(lattice_dimensions,
                    xycell_position,[4,1],bonds,i)
            bonds = self.add_cross_bond(lattice_dimensions,
                    ycell_position,[4,1],bonds,i)
            
            bonds = self.add_cross_bond(lattice_dimensions,
                    xcell_position,[7,6],bonds,i)
            bonds = self.add_cross_bond(lattice_dimensions,
                    xcell_position,[8,5],bonds,i)
            bonds = self.add_cross_bond(lattice_dimensions,
                    ycell_position,[7,8],bonds,i)

        return bonds

    def angles(self,bonds):
        angles = np.empty((0,3),dtype=int)
        N = int(np.amax(bonds)) # Number of atoms
        
        for centre in range(1,N+1):
            neighbours = find_neighbours(bonds,centre)
            for i in range(len(neighbours)):
                for j in range(i+1,len(neighbours)):
                    angle = [neighbours[i],centre,neighbours[j]]
                    angles = np.vstack((angles,angle))
        return angles

    def torsions(self,bonds):
        torsions = np.empty((0,4),dtype=int)
        for bond in bonds:
            neighbours1 = find_neighbours(bonds,bond[0])
            neighbours1.remove(bond[1])
            neighbours2 = find_neighbours(bonds,bond[1])
            neighbours2.remove(bond[0])
            if len(neighbours1) and len(neighbours2):
                for neighbour1 in neighbours1:
                    for neighbour2 in neighbours2:
                        torsion = [neighbour1,bond[0],
                                   bond[1],neighbour2]
                        torsions = np.vstack((torsions,torsion))
        return torsions       

    def add_cross_bond(self,lattice_dimensions,
                   cell_position,atoms,bonds,i):
        # cell_position : the adjoining cell
        # atoms [i,j]: ith atom in cell, jth atom in adjoining
        atoms[0] += i
        atoms[1] += self.index_cell(cell_position,
                                   lattice_dimensions,8)
        return np.vstack((bonds,atoms))



def find_connections(bonds,centre):
    connections = np.where(bonds==centre)
    connections =np.vstack((connections[0],connections[1]))
    return connections.transpose()
            
def find_neighbours(bonds,centre):
    connections = find_connections(bonds,centre)
    neighbours = []
    for connection in connections:
        # Find atom connected to centre
        neighbour = bonds[connection[0]][connection[1]-1]
        neighbours.append(neighbour)
    return neighbours

def find_adjacent_cells(cell_position,lattice_dimensions):
    x,y,z = cell_position
    # Account for periodic system 
    if x == lattice_dimensions[0]-1: 
        xcell_position=[0,y,z]
    else: xcell_position = [x+1,y,z]
    
    if y == lattice_dimensions[1]-1: 
        ycell_position=[x,0,z]
    else: ycell_position = [x,y+1,z]          
                                                    
    xycell_position = [xcell_position[0],   
                       ycell_position[1],z]

    return xcell_position, ycell_position, xycell_position
       






