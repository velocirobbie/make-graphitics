import numpy as np
from math import pi

class Island(object):
    def __init__(self):
        self.atoms = []
        self.coords = []
        self.area = None
        self.boundary = None

    def populate_coords(self,all_coords):
        N = len(self.atoms)
        self.coords = np.empty((N,3))
        for i,atom in enumerate(self.atoms):
            self.coords[i] = all_coords[atom-1]

    def calc_boundary(self,alpha):
        self.boundary = outer_polygon(self.coords[:,0:2], alpha)

    def calc_polygon_area(self):
        self.area = self.boundary.area

    def com(self):
        return self.coords[:,0:2].mean(0)

    def natoms(self):
        return len(self.atoms)

    def diameter(self):
        return 2 * np.sqrt(self.area/pi)

def build_bond_network(bonds, atom_types):
    N = len(atom_types)
    bond_network = {i+1:{'type':type_,'bonded_to':[]} for i,type_ in enumerate(atom_types)}
    for bond in bonds:
        bond_network[bond[0]]['bonded_to'] += [bond[1]]
        bond_network[bond[1]]['bonded_to'] += [bond[0]]
    return bond_network

def flood_island(index, bond_network, island_labels, atom_types, island_index):
    import sys
    sys.setrecursionlimit(10000)
    island = Island()

    def add_neighbours(atom):
        island.atoms += [atom]
        island_labels[atom-1] = island_index
        neighbours = bond_network[atom]['bonded_to']
        for neighbour in neighbours:
            atom_type = atom_types[neighbour-1]
            already_included = island_labels[neighbour-1]
            if (atom_type==1) and (not already_included):
                add_neighbours(neighbour)
    add_neighbours(index+1)

    return island, island_labels

def find_islands_by_flood(sim):
    bonds = sim.bonds #2xN array
    atom_types = sim.atom_labels

    N = len(atom_types)
    bond_network = build_bond_network(bonds, atom_types)
    islands = []

    # array recording if atoms are in an island
    # 0=not an island (graphene); 1,2,3.. = in island N
    island_labels = np.zeros(N)

    for i in range(N):
        # if not already counted in an island, and a aromatic carbon
        if (island_labels[i] == 0) and (atom_types[i] == 1) :
            # found new island
            island_index = len(islands)+1
            island, island_labels = flood_island(i, bond_network, island_labels,atom_types, island_index)
            islands += [island]

    map(lambda island: island.populate_coords(sim.coords), islands)

    return islands

def outer_polygon(coords,alpha):
    if len(coords) < 3: return 0
    from scipy.spatial import Delaunay
    from shapely.ops import cascaded_union, polygonize
    import shapely.geometry as geometry

    def alpha_shape(coords, alpha):
        def add_edge(edges, edge_points, coords, i, j):
            """
            Add a line between the i-th and j-th points,
            if not in the list already
            """
            if (i, j) in edges or (j, i) in edges:
                # already added
                return
            edges.add( (i, j) )
            edge_points.append(coords[ [i, j] ])
        tri = Delaunay(coords)
        edges = set()
        edge_points = []
        # loop over triangles:
        # ia, ib, ic = indices of corner points of the
        # triangle
        for ia, ib, ic in tri.vertices:
            pa = coords[ia]
            pb = coords[ib]
            pc = coords[ic]
            # Lengths of sides of triangle
            a = np.sqrt((pa[0]-pb[0])**2 + (pa[1]-pb[1])**2)
            b = np.sqrt((pb[0]-pc[0])**2 + (pb[1]-pc[1])**2)
            c = np.sqrt((pc[0]-pa[0])**2 + (pc[1]-pa[1])**2)
            if a > alpha or b > alpha or c > alpha:
                pass
            else:
                add_edge(edges, edge_points, coords, ia, ib)
                add_edge(edges, edge_points, coords, ib, ic)
                add_edge(edges, edge_points, coords, ic, ia)
        m = geometry.MultiLineString(edge_points)
        triangles = list(polygonize(m))
        return cascaded_union(triangles), edge_points

    if len(coords) < 4:
        # When you have a triangle, there is no sense
        # in computing an alpha shape.
        return geometry.MultiPoint(list(coords)).convex_hull
    else:
        concave_hull, edge_points = alpha_shape(coords, alpha=alpha)
        return concave_hull

def simple_area(coords):
    bond_length = 1.42 # Angstroms, in graphene
    atom_area = pi * (1.42/2)**2
    return len(coords) * atom_area

def strip_small_islands(islands,min_atoms):
    new_islands = []
    for island in islands:
        if island.natoms() >= min_atoms:
            new_islands += [island]
    return new_islands

def write_islands_xyz(islands):
    N = sum([island.natoms() for island in islands])

    with open('islands.xyz','w') as f:
        f.write(str(N)+'\n')
        f.write('islands\n')
        for i,island in enumerate(islands):
            for coord in island.coords:
                f.write(str((i+1)%10)+'\t')
                for axis in coord:
                    f.write(str(axis)+'\t')
                f.write('\n')

def write_gnuplot(islands):
    # write object file for polygons
    with open('island_objects.sh','w') as f:
        for i,island in enumerate(islands):
            f.write("set object "+str(i+1)+" polygon from \\\n")
            coords = list(island.boundary.exterior.coords)
            for point in coords:
                f.write("\t"+str(point[0])+","+str(point[1])+" to \\\n")
            f.write("\t"+str(coords[-1][0])+","+str(coords[-1][1])+" \n")
            f.write("set object "+str(i+1)+" fc rgb '#000000' fillstyle solid lw 0\n\n")
    # write coords of all island atoms
    with open('island_coords.dat','w') as f:
        for island in islands:
            for coord in island.coords:
                f.write(str(coord[0])+"\t "+str(coord[1])+"\n")
    # draw circles of an island's approximate area
    with open('island_area_circle.dat','w') as f:
        for island in islands:
            f.write(str(island.com()[0])+"\t "+str(island.com()[1])+"\t "+str(island.diameter()/2)+"\n")
    # plot with:
    # gnuplot> load 'island_objects.sh'
    # gnuplot> pl 'island_coords.dat' ,'island_area_circle.dat' w circ

def calc_island_sizes(sim):
    islands = find_islands_by_flood(sim)
    islands = strip_small_islands(islands,6)
    write_islands_xyz(islands)

    map(lambda island: island.calc_boundary(3), islands)
    map(lambda island: island.calc_polygon_area(), islands)
    write_gnuplot(islands)

    sizes = [island.diameter() for island in islands]
    return sizes

