from process_coronene import process
from scripts.write_coords import Writer
from scripts.shifty import Shifter

system_name=        'Driedling - coronene zshift'
crystal_parameters= 'Driedling' # listed in config file

data = process(crystal_parameters) + [system_name]

#output = Writer(*data)
#output.write_xyz('coronene.xyz')
#output.write_lammps('data.bulk')
    
data[-2][2] = 100  # set z spacing 
output = Shifter('top','lammps',data)
#output.in_plane_shift([1,0],-5,5,0.1)
output.z_shift(-0.5,3,0.1)
output.z_shift(3,10,0.5)


