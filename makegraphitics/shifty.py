import numpy as np
from math import cos, sin, pi
from write_coords import Writer


class Shifter(object):
    def __init__(self, sim, output_style="xyz", target=0):
        self.sim = sim
        self.coords = sim.coords
        self.mol = sim.molecule_labels
        self.target = target
        self.output_style = output_style

        if not target:
            self.target = np.amax(self.mol)

    def rotate(self, end, steps):
        rotate_range = np.arange(0, end, steps)
        for rotation in rotate_range:
            rotated_coords = self.rotate_molecule(rotation * pi / 180.0)
            self.write_shifted_coords(rotated_coords, rotation)

    def rotate_molecule(self, theta):
        sint = sin(theta)
        cost = cos(theta)
        new_coords = np.empty(np.shape(self.coords))
        for i in range(len(self.coords)):
            if self.mol[i] == self.target:
                a = cost * self.coords[i][0] - sint * self.coords[i][1]
                b = sint * self.coords[i][0] + cost * self.coords[i][1]
                c = self.coords[i][2]
                new_coords[i] = [a, b, c]
            else:
                new_coords[i] = self.coords[i]
        return new_coords

    def z_shift(self, start, end, step):
        shift_range = np.arange(start, end, step)
        for shift in shift_range:
            shifted_coords = self.move_molecule(0, 0, shift)
            self.write_shifted_coords(shifted_coords, shift)

    def in_plane_shift(self, direction, start, end, step):
        # direction form [x,y] - should be a unit vector
        shift_range = np.arange(start, end, step)
        for shift in shift_range:
            shifted_coords = self.move_molecule(
                shift * direction[0], shift * direction[1], 0
            )
            self.write_shifted_coords(shifted_coords, shift)

    def move_molecule(self, x_shift, y_shift, z_shift):
        new_coords = np.empty(np.shape(self.coords))
        for i in range(len(self.coords)):
            if self.mol[i] == self.target:
                new_coords[i][0] = self.coords[i][0] + x_shift
                new_coords[i][1] = self.coords[i][1] + y_shift
                new_coords[i][2] = self.coords[i][2] + z_shift
            else:
                new_coords[i] = self.coords[i]
        return new_coords

    def write_shifted_coords(self, shifted_coords, shift):
        temp_sim = self.sim
        temp_sim.coords = shifted_coords
        writer = Writer(temp_sim)
        if self.output_style == "xyz":
            writer.write_xyz("shift_" + str(shift) + ".xyz")
        elif self.output_style == "lammps":
            writer.write_lammps("data.shift_" + str(shift))
