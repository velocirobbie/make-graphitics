import numpy as np
from write_coords import Writer

class Shifter(object):
    def __init__(
            self,coords,molecule_labels,target,
            output_style,cell_dimensions):
        
        self.coords = coords
        self.mol = molecule_labels
        self.target = target
        self.output_style = output_style
        self.cell = cell_dimensions

    def z_shift(self,start,end,step):
        shift_range = np.arange(start,end,step)
        for shift in shift_range:
            shifted_coords = self.move_molecule(0,0,shift)
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
        writer = Writer(shifted_coords,self.mol)
        if self.output_style == 'xyz':
            writer.write_xyz('shift_'+str(shift)+'.xyz',str(shift))
        elif self.output_style == 'lammps':
            writer.write_lammps(self.cell,'data.shift_'+str(shift))

