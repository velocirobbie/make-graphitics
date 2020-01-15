# Graphite, Graphene and Graphene Oxide Builder

MakeGraphitics is a library to ceate various atomistic graphitic structures for molecular dynamics.

Available structures:
- Hexagonal graphene flake
- Rectangular graphene flake
- Rectangular perodic graphene sheet (no edges)
- Periodic graphite 
- Graphene and graphite oxide

Output:
- .xyz 
- lammps data file

Automatically parameterise by forcefields:
- OPLS
- GraFF
- ReaxFF

## Install

Clone this repository. Install using Python2.7. Run the tests to check the installation has worked.
```
git clone https://github.com/velocirobbie/make-graphitics
cd make-graphitics
python setup.py install
pytest
```

A conda environment is provided if you do not have the right packages. If you have conda set up, execute these commands to create a working python environment before the install setp.
```
conda env create --file graphene-env.yml
conda activate graphene
```

## Examples

See the scripts in the `examples/` directory for a number of sample structures.

1) Make a rectangular graphene sheet that extends through periodic boundaries. Parameterised with OPLS and outputs to .xyz for easy veiwing with VMD and a LAMMPS data file.
```
python2.7 graphene_sheet.py
```
Size of the sheet can be specified in `graphene_sheet.py`.

2) Make a hexagonal flake of graphene oxide. Parameterised with OPLS and outputs to .xyz for easy veiwing with VMD and a LAMMPS data file.
```
python2.7 GO_flake.py
```
There are several tunable parameters in `GO_flake.py` that you may be interested in. Including:
- flake radius
- C/O target ratio, `ratio`
- Rate at which new nodes are added, `new_island_freq`
- output snapshots of the oxidation process every N steps with `video_xyz=N`. Viewed in VMD with `topo readvarxyz out.xyz`

## Notes on the Oxidiser 

The Oxidiser takes a graphitic structure and attempts to oxidise it by the process described in (Yang, Angewandte Chemie, 2014; Sinclair, 2019). The algorithm proceeds as follows:

1) If hydrogens exist (i.e. edge of a flake), 1/4 are changed to alcohol groups and 1/4 to carboxyl groups (Lerf and Klinowski model). These values can be changed by passing the Oxidiser object the optional arguments: ` edge_OHratio = 0.25, edge_carboxyl_ratio = 0.5`.

2) The reactivity of every possible site is calculated. This is done by using a ranodom forest approach to extend the data set of GO reactivites given by Yang et al. We do not take into account the reactivity of the edges.

3) A site is oxidised at random weighted by each site's reactivity. The chance of an oxidisation producing an alcohol or epoxy group on the surface is by default 50:50, but can be specified by passing Oxidiser the optional argument: `surface_OHratio = 0.5`

4) The time elapsed between oxidations is estimated from the reactivity of the site that has been oxidised.

5) New nodes are added proportionally to `time_elapsed * new_island_freq`. Note this can be 0. The reasons for doing this are outlined in (Sinclair 2019).

6) Steps 2-5 are repeated until the target C/O ratio is reached or no new sites are available, usually C/O ~ 1.7 . We recommend setting the target, `ratio`, to over 2 as this is what is seen experimentally.

## Notes on Parameterisation

Not all the bonded interactions that can occur in graphene oxide are included in the OPLS parameterisation. We make some neccesary like for like atom-type substitutions to get around this problem. It is not ideal but common practice in molecular dynamics. The substitutions used are outputed after a parameterisation step. Each substitution line outputs the origional atom types, the atom types used to parameterise them, and a summary string that you can use to find in the script `scripts/params.py`. Substitutions keep atom types as close to the origional as possible e.g. replaces an aromatic C with an alkene C, whcih are both sp2 carbon atoms. 

## More structure examples

More examples of building structures with this script are in the `examples` directory. Files must be copied to the source directory and executed there.

Note that differenct structures can be combined into one simulation object with `Combine`. Also coordinates can be manipulated before writing to a lammps file. An examploe of this is shown in `peel_sim.py`.

# Citing

The work contained here has been published in some of my own papers e.g.

 - Graphene–Graphene Interactions: Friction, Superlubricity, and Exfoliation https://doi.org/10.1002/adma.201705791

 - Modeling Nanostructure in Graphene Oxide: Inhomogeneity and the Percolation Threshold https://doi.org/10.1021/acs.jcim.9b00114

 - The Role of Graphene in Enhancing the Material Properties of Thermosetting Polymers https://doi.org/10.1002/adts.201800168

I would appreciate a citation if you any of the code in any published work :) You could cite the graphene oxide structure paper, this github page (if the journal allows), or the latest release on the zenodo repository

```
@article{sinclair2019modelling,
  title={Modelling nanostructure in graphene oxide: inhomogeneity and the percolation threshold},
  author={Sinclair, Robert Callum and Coveney, Peter Vivian},
  journal = {Journal of Chemical Information and Modeling},
  volume = {59},
  number = {6},
  pages = {2741-2745},
  year = {2019},
  doi = {10.1021/acs.jcim.9b00114},
}
@misc{make-graphitics-github,
    url = {https://github.com/velocirobbie/make-graphitics},
    howpublished = {\url{https://github.com/velocirobbie/make-graphitics}},
    note = {Accessed: \today},
    author = {Sinclair, Robert C.},
    year = {2019}
}
@misc{make-graphitics_zenodo,
  author    = {Sinclair, Robert. C. },
  title     = {make-graphitics},
  version   = {0.1.0},
  publisher = {Zenodo},
  year      = {2019},
  doi       = {10.5281/zenodo.2548538}
}

