from scripts import *
from analisys import *
import numpy as np
import sys

go_data_file = sys.argv[1]

go = ReadLammpsData(go_data_file)

island_sizes = calc_island_sizes(go)
print island_sizes
