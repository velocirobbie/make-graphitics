from scripts import *
from analisys import calc_island_size
import numpy as np
import sys

go_data_file = sys.argv[1]

go = ReadLammpsData(go_data_file)


