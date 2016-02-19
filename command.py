from process import process
from scripts.write_coords import Writer
from scripts.shifty import Shifter

system_name=        'Comment line'
crystal_parameters= 'OPLS' # listed in config file
layer_or_bulk=      'bulk'

data = process(crystal_parameters, layer_or_bulk) + [system_name]

#output = Writer(*data)
#output.write_xyz('bulk.xyz')
#output.write_lammps('data.bulk')
    
data[5][2] = 100  # set z spacing 
output = Shifter('top','xyz',*data)
#output.in_plane_shift([1,0],0,3.0,0.02)
output.z_shift(-0,10,1)

