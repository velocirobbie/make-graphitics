import os
from opls_reader import OPLS_Reader


class Parameterise(object):
    def __init__(self, crystal, vdw_defs=None, forcefield="OPLS"):
        # use vdw_defs of atoms to paramterise the bonded and nonbonded coefficients
        if not vdw_defs:
            try:
                self.vdw_defs = crystal.vdw_defs
            except AttributeError:
                raise Exception(
                    "The simulation you want to parameterise requires vdw_defs ,"
                    + "provide as an attribute to the simulation or as a keyword arg"
                )
        else:
            self.vdw_defs = vdw_defs
            crystal.vdw_defs = vdw_defs

        self.retrieve_ff_data(forcefield)

        self.type_defs = {}
        for label in self.vdw_defs:
            vdw = self.vdw_defs[label]
            self.type_defs[label] = self.vdw_type["type"][
                self.vdw_type["vdw"].index(vdw)
            ]
        print "Atom label -> OPLS vdw definitions: \t", self.vdw_defs
        print "Atom label -> OPLS type definitions: \t", self.type_defs

        try:
            crystal.bond_types, crystal.angle_types, crystal.dihedral_types, crystal.improper_types
        except AttributeError:
            crystal.generate_connections()
        crystal.bond_coeffs = self.match_bonds(crystal.bond_types)
        crystal.angle_coeffs = self.match_angles(crystal.angle_types)
        crystal.dihedral_coeffs = self.match_dihedrals(crystal.dihedral_types)
        crystal.improper_coeffs = self.match_impropers(crystal.improper_types)
        crystal.pair_coeffs = self.match_pairs()
        crystal.masses = self.match_masses()

    def match_charges(self):
        charge_data = self.charge_data
        charge_coeffs = {"a": [], "q": []}
        for label in self.vdw_defs:
            found = 0
            for j in range(len(charge_data["a"])):
                atom_data = charge_data["a"][j]
                if self.vdw_defs[label] == atom_data:
                    found += 1
                    charge_coeffs["a"] += [label]
                    charge_coeffs["q"] += [charge_data["q"][j]]
            if found != 1:
                raise ValueError("WRONG", label, "\t found ", found, " entries")
        return charge_coeffs

    def match_masses(self):
        mass_data = self.mass_data
        mass_coeffs = {}
        for label in self.vdw_defs:
            found = 0
            for j in range(len(mass_data["a"])):
                atom_data = mass_data["a"][j]
                if self.vdw_defs[label] == atom_data:
                    found += 1
                    mass_coeffs[label] = mass_data["m"][j]
            if found != 1:
                raise ValueError("WRONG", label, "\t found ", found, " entries")
        return mass_coeffs

    def match_pairs(self):
        pair_data = self.pair_data
        pair_coeffs = {}
        for label in self.vdw_defs:
            found = 0
            for j in range(len(pair_data["a"])):
                atom_data = pair_data["a"][j]
                if self.vdw_defs[label] == atom_data:
                    found += 1
                    pair_coeffs[label] = {}
                    pair_coeffs[label][1] = pair_data["e"][j]
                    pair_coeffs[label][2] = pair_data["s"][j]
            if found != 1:
                raise ValueError("WRONG", label, "\t found ", found, " entries")
        return pair_coeffs

    def match_bonds(self, bond_types):
        bond_data = self.bond_data
        bond_coeffs = {}
        for i in range(len(bond_types)):
            a1 = self.type_defs[bond_types[i][0]]
            a2 = self.type_defs[bond_types[i][1]]
            atoms = [a1, a2]
            found = 0
            for j in range(len(bond_data["k"])):
                atom_data = [bond_data["a1"][j], bond_data["a2"][j]]
                flag1 = atoms == atom_data
                flag2 = atoms == list(reversed(atom_data))
                if flag1 or flag2:
                    found += 1
                    bond_coeffs[i + 1] = {}
                    bond_coeffs[i + 1][1] = bond_data["k"][j]
                    bond_coeffs[i + 1][2] = bond_data["r"][j]
            if found != 1:
                raise ValueError("WRONG", atoms, "\t found ", found, " entries")
        return bond_coeffs

    def match_angles(self, angle_types):
        def search_angles(angle_data, atoms, angle_coeffs):
            found = 0
            for j in range(len(angle_data["k"])):
                atom_data = [
                    angle_data["a1"][j],
                    angle_data["a2"][j],
                    angle_data["a3"][j],
                ]
                flag1 = atoms == atom_data
                flag2 = atoms == list(reversed(atom_data))
                if flag1 or flag2:
                    found += 1
                    angle_coeffs[i + 1] = {}
                    angle_coeffs[i + 1][1] = angle_data["k"][j]
                    angle_coeffs[i + 1][2] = angle_data["r"][j]
            return angle_coeffs, found

        angle_data = self.angle_data
        angle_coeffs = {}
        questionable_substitutions = 0
        for i in range(len(angle_types)):
            a1 = self.type_defs[angle_types[i][0]]
            a2 = self.type_defs[angle_types[i][1]]
            a3 = self.type_defs[angle_types[i][2]]
            atoms = [a1, a2, a3]
            # print atoms, angle_types[i]
            angle_coeffs, found = search_angles(angle_data, atoms, angle_coeffs)

            if found == 0:
                for j in range(3):
                    if atoms[j] == 48:
                        atoms[j] = 47  # Try alkene carbon instead of aromatic
                    if atoms[j] == 49:
                        atoms[j] = 46  # Try alkene H
                angle_coeffs, found = search_angles(angle_data, atoms, angle_coeffs)
                if found == 1:
                    print found, [a1, a2, a3], atoms, "sub aromatic -> alkene"
            if found == 0:
                if atoms[1] == 47:
                    if atoms[0] == 13:
                        atoms[0] = 47
                    elif atoms[2] == 13:
                        atoms[0] = 47
                    angle_coeffs, found = search_angles(angle_data, atoms, angle_coeffs)
                    if found == 0:
                        print found, [a1, a2, a3], atoms, "made 13 -> 47 sub"
            if found != 1:
                raise ValueError("WRONG", atoms, "\t found ", found, " entries")
                # print 'WRONG',atoms,'\t found ',found,' entries'
        if questionable_substitutions != 0:
            print "made ", questionable_substitutions, " questionable angle subs"
        return angle_coeffs

    def match_dihedrals(self, dihedral_types):
        def add_coeff(dihedral_data, dihedral_coeffs, j):
            N = len(dihedral_coeffs)
            dihedral_coeffs[N + 1] = {}
            dihedral_coeffs[N + 1][1] = dihedral_data["k1"][j]
            dihedral_coeffs[N + 1][2] = dihedral_data["k2"][j]
            dihedral_coeffs[N + 1][3] = dihedral_data["k3"][j]
            dihedral_coeffs[N + 1][4] = dihedral_data["k4"][j]
            return dihedral_coeffs

        def check_wildcards(atoms, dihedral_data, dihedral_coeffs):
            found = 0
            t = len(dihedral_data["k1"])
            for i in range(t):
                a1 = dihedral_data["a1"][i]
                a2 = dihedral_data["a2"][i]
                a3 = dihedral_data["a3"][i]
                a4 = dihedral_data["a4"][i]
                if a1 == 0:
                    flag1 = atoms[1:4] == [a2, a3, a4]
                    flag2 = atoms[0:3] == [a4, a3, a2]
                    if flag1 or flag2:
                        found += 1
                        dihedral_coeffs = add_coeff(dihedral_data, dihedral_coeffs, i)
                        break
                if a4 == 0:
                    flag1 = atoms[0:3] == [a1, a2, a3]
                    flag2 = atoms[1:4] == [a4, a2, a1]
                    if flag1 or flag2:
                        foudn += 1
                        dihedral_coeffs = add_coeff(dihedral_data, dihedral_coeffs, i)
                        break
                if a1 == 0 and a4 == 0:
                    flag1 = [atoms[1], atoms[2]] == [a2, a3]
                    flag2 = [atoms[2], atoms[1]] == [a3, a2]
                    if flag1 or flag2:
                        found += 1
                        dihedral_coeffs = add_coeff(dihedral_data, dihedral_coeffs, i)
                        break
            return dihedral_coeffs, found

        def search_dihedrals(dihedral_data, atoms, dihedral_coeffs):
            found = 0
            for j in range(len(dihedral_data["k1"])):
                atom_data = [
                    dihedral_data["a1"][j],
                    dihedral_data["a2"][j],
                    dihedral_data["a3"][j],
                    dihedral_data["a4"][j],
                ]
                flag1 = atoms == atom_data
                flag2 = atoms == list(reversed(atom_data))
                if flag1 or flag2:
                    found += 1
                    dihedral_coeffs = add_coeff(dihedral_data, dihedral_coeffs, j)
            if found == 0:
                dihedral_coeffs, found = check_wildcards(
                    atoms, dihedral_data, dihedral_coeffs
                )
            return dihedral_coeffs, found

        dihedral_data = self.dihedral_data
        dihedral_coeffs = {}
        questionable = 0
        for i in range(len(dihedral_types)):
            a1 = self.type_defs[dihedral_types[i][0]]
            a2 = self.type_defs[dihedral_types[i][1]]
            a3 = self.type_defs[dihedral_types[i][2]]
            a4 = self.type_defs[dihedral_types[i][3]]
            atoms = [a1, a2, a3, a4]
            # print atoms

            dihedral_coeffs, found = search_dihedrals(
                dihedral_data, atoms, dihedral_coeffs
            )

            if found == 0:
                for j in range(4):
                    if atoms[j] == 48:
                        atoms[j] = 47  # Try alkene carbon instead of aromatic
                    if atoms[j] == 49:
                        atoms[j] = 46  # Try alkene H
                dihedral_coeffs, found = search_dihedrals(
                    dihedral_data, atoms, dihedral_coeffs
                )
            if found == 0:
                change = False
                if atoms[0] == 5:
                    atoms[0] = 46
                    change = True
                if atoms[3] == 5:
                    atoms[3] = 46
                dihedral_coeffs, found = search_dihedrals(
                    dihedral_data, atoms, dihedral_coeffs
                )
                if found != 0:
                    print found, [a1, a2, a3, a4], atoms, "made 5 -> 46 sub"

            if found == 0:
                if atoms[1:4] == [13, 47, 3] or atoms[0:3] == [3, 47, 13]:
                    atoms = [0, 13, 47, 47]
                    dihedral_coeffs, found = search_dihedrals(
                        dihedral_data, atoms, dihedral_coeffs
                    )
                if found != 0:
                    print found, [a1, a2, a3, a4], atoms, "made 13,47,47 sub"

            if found == 0:
                if atoms[1] == 13 and atoms[2] == 13:
                    if atoms[0] == 47:
                        atoms[0] = 3
                    elif atoms[3] == 47:
                        atoms[0] = 3
                    dihedral_coeffs, found = search_dihedrals(
                        dihedral_data, atoms, dihedral_coeffs
                    )
                    if found != 0:
                        print found, [a1, a2, a3, a4], atoms, "sub 47 -> 3"

            if found == 0:
                if atoms == [13, 47, 5, 7] or atoms == [7, 5, 47, 13]:
                    atoms = [7, 5, 47, 47]
                    dihedral_coeffs, found = search_dihedrals(
                        dihedral_data, atoms, dihedral_coeffs
                    )
                    if found != 0:
                        print found, [a1, a2, a3, a4], atoms, "made phenol sub"

            if found != 1:
                raise ValueError(
                    "Torsion", [a1, a2, a3, a4], "found ", found, " entries"
                )
                # print 'WRONG',[a1,a2,a3,a4],'\t found ',found,' entries'
        if questionable != 0:
            print "made ", questionable, " questionable dihedral substitutions"
        return dihedral_coeffs

    def match_impropers(self, improper_types):
        def search_impropers(improper_data, centre, neighbours, improper_coeffs):
            found = 0
            for j in range(len(improper_data["k"])):
                data_neighbours = [
                    improper_data["a1"][j],
                    improper_data["a2"][j],
                    improper_data["a3"][j],
                ]
                data_centre = improper_data["centre"][j]
                flag1 = data_centre == centre
                # Test if other improper atom in connected to centre
                flag2 = set(neighbours) >= {data_neighbours[2]}
                flag3 = data_neighbours == [0, 0, 0]
                if flag1 and (flag2 or flag3):
                    found += 1
                    improper_coeffs[i + 1] = {}
                    improper_coeffs[i + 1][1] = improper_data["k"][j]
                    improper_coeffs[i + 1][2] = improper_data["r"][j]
            return improper_coeffs, found

        improper_data = self.improper_data
        improper_coeffs = {}
        questionable_substitutions = 0
        for i in range(len(improper_types)):
            centre = self.type_defs[improper_types[i][0]]
            a1 = self.type_defs[improper_types[i][1]]
            a2 = self.type_defs[improper_types[i][2]]
            a3 = self.type_defs[improper_types[i][3]]
            neighbours = [a1, a2, a3]

            improper_coeffs, found = search_impropers(
                improper_data, centre, neighbours, improper_coeffs
            )
            if found != 1:
                raise ValueError(
                    "Improper", centre, neighbours, "found ", found, " entries"
                )
                # print 'WRONG',centre,neighbours,'\t found ',found,' entries'
        return improper_coeffs

    def retrieve_ff_data(self, forcefield):
        if forcefield == "OPLS":
            paramfile = os.path.dirname(__file__) + "/params/oplsaa.prm"
            data = OPLS_Reader(paramfile)
        elif forcefield == "ReaxFF":
            # paramters can be found https://github.com/lammps/lammps/tree/master/potentials
            raise Exception("No need to generate paramters for this forcefield")
        else:
            raise Exception(forcefield, "forcefield not implemented")
        self.vdw_type = data.vdw_type
        self.bond_data = data.bond
        self.angle_data = data.angle
        self.dihedral_data = data.dihedral
        self.improper_data = data.improper
        self.pair_data = data.pair
        self.mass_data = data.mass
        self.charge_data = data.charge
