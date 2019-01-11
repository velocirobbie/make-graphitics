# Graphite Graphene and Graphene Oxide Builder

These scripts can create many graphitic structures for use in atomistic modelling.

Available structures:
- Hexagonal graphene flake
- Rectangular perodic graphene sheet
- Periodic graphite 
- Graphene and graphite oxide

Output:
- .xyz 
- lammps data file

Automatically parameterise by forcefields:
- OPLS
- GraFF

## Examples

1) Make a rectangular graphene sheet that extends through periodic boundaries. Parameterised with OPLS and outputs to .xyz for easy veiwing with VMD and a LAMMPS data file.
``python2.7 make_graphene_sheet.py``
Size of the sheet can be specified in `make_graphene_sheet.py`.

2) Make a hexagonal flake of graphene oxide. Parameterised with OPLS and outputs to .xyz for easy veiwing with VMD and a LAMMPS data file.
``python2.7 make_GO_flake.py``
There are several tunable parameters in `make_GO_flake.py` including
- flake radius
- C/O target ratio, `ratio`
- Rate at which new nodes are added, `new_island_freq`
- output snapshots of the oxidation process every N steps with `video=N`. Viewed in VMD with `topo readvarxyz out.xyz`


