import sys


class OPLS_Reader(object):
    def __init__(self, datafile):
        self.setup_dicts()
        self.main = {
            "atom    ": self.add_mass,
            "vdw     ": self.add_pair,
            "bond    ": self.add_bond,
            "angle   ": self.add_angle,
            "torsion ": self.add_dihedral,
            "charge  ": self.add_charge,
            "imptors ": self.add_improper,
        }

        with open(datafile) as f:
            for line in f:
                if len(line) > 7:
                    self.readline(line)
                else:
                    pass
                    # blank line
                    # print 'BLANK',line

    def readline(self, line):
        func = self.main.get(line[0:8], False)
        if not func:
            pass
            # data we can't parse
            # print 'WRONG',line
        else:
            func(line)

    def add_mass(self, line):
        self.vdw_type["vdw"] += [int(line[12:15])]
        self.vdw_type["type"] += [int(line[15:20])]
        self.mass["a"] += [int(line[12:15])]
        self.mass["m"] += [float(line[65:72])]

    def add_pair(self, line):
        line = line.split()
        self.pair["a"] += [int(line[1])]
        self.pair["s"] += [float(line[2])]
        self.pair["e"] += [float(line[3])]

    def add_bond(self, line):
        line = line.split()
        self.bond["a1"] += [int(line[1])]
        self.bond["a2"] += [int(line[2])]
        self.bond["k"] += [float(line[3])]
        self.bond["r"] += [float(line[4])]

    def add_angle(self, line):
        line = line.split()
        self.angle["a1"] += [int(line[1])]
        self.angle["a2"] += [int(line[2])]
        self.angle["a3"] += [int(line[3])]
        self.angle["k"] += [float(line[4])]
        self.angle["r"] += [float(line[5])]

    def add_dihedral(self, line):
        line = line.split()
        self.dihedral["a1"] += [int(line[1])]
        self.dihedral["a2"] += [int(line[2])]
        self.dihedral["a3"] += [int(line[3])]
        self.dihedral["a4"] += [int(line[4])]
        self.dihedral["k1"] += [float(line[5])]
        self.dihedral["k2"] += [float(line[8])]
        self.dihedral["k3"] += [float(line[11])]
        if len(line) == 17:
            self.dihedral["k4"] += [float(line[14])]
        else:
            self.dihedral["k4"] += [0.0]

    def add_charge(self, line):
        line = line.split()
        self.charge["a"] += [int(line[1])]
        self.charge["q"] += [float(line[2])]

    def add_improper(self, line):
        line = line.split()
        self.improper["a1"] += [int(line[1])]
        self.improper["a2"] += [int(line[2])]
        self.improper["centre"] += [int(line[3])]
        self.improper["a3"] += [int(line[4])]
        self.improper["k"] += [float(line[5])]
        self.improper["r"] += [float(line[6])]

    def setup_dicts(self):
        self.vdw_type = {"vdw": [], "type": []}
        self.pair = {"a": [], "s": [], "e": []}
        self.bond = {"a1": [], "a2": [], "k": [], "r": []}
        self.angle = {"a1": [], "a2": [], "a3": [], "k": [], "r": []}
        self.dihedral = {
            "a1": [],
            "a2": [],
            "a3": [],
            "a4": [],
            "k1": [],
            "k2": [],
            "k3": [],
            "k4": [],
        }
        self.improper = {"a1": [], "a2": [], "centre": [], "a3": [], "k": [], "r": []}
        self.charge = {"a": [], "q": []}
        self.mass = {"a": [], "m": []}


# a = OPLS_Reader(sys.argv[1])
