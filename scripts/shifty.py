import numpy as np
from write_coords import Writer

class Shifter(object):
    def __init__(self, target, output_style,data):
        self.data = data
        self.coords = data[0]
        self.mol = data[3]
        self.target = target
        self.output_style = output_style
        
        if target == 'top':
            self.target = np.amax(data[3])

    def z_shift(self,start,end,step):
        shift_range = np.arange(start,end,step)
        for shift in shift_range:
            shifted_coords = self.move_molecule(0,0,shift)
            self.write_shifted_coords(shifted_coords,shift)
    
    def in_plane_shift(self,direction,start,end,step):
        # direction form [x,y] - should be a unit vector
        shift_range = np.arange(start,end,step)
        for shift in shift_range:
            shifted_coords = self.move_molecule(
                    shift*direction[0],shift*direction[1],0)
            self.write_shifted_coords(shifted_coords,shift)

    def move_molecule(self,x_shift,y_shift,z_shift):
        new_coords = np.empty(np.shape(self.coords))
        for i in range(len(self.coords)):
            if self.mol[i] == self.target:
                new_coords[i][0] = self.coords[i][0] + x_shift
                new_coords[i][1] = self.coords[i][1] + y_shift
                new_coords[i][2] = self.coords[i][2] + z_shift
            else:
                new_coords[i] = self.coords[i]
        return new_coords

    def write_shifted_coords(self,shifted_coords,shift):
        temp_data = self.data
        temp_data[0] = shifted_coords
        writer = Writer(*temp_data)
        if self.output_style == 'xyz':
            writer.write_xyz('shift_'+str(shift)+'.xyz')
        elif self.output_style == 'lammps':
            writer.write_lammps('data.shift_'+str(shift))

