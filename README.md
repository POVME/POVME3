# POVME

POVME (Pocket VOlume MEasurer) is a Python package for extracting actionable information from ensembles of protein structures for use in drug design.

Medicinal chemists might use this to find representative confirmations of a binding pocket that encompass the diversity found in a simulation or ensemble of structures.

POVME is a way to compress a large amount of protein confirmations into clusters or families

POVME is a tool to map and navigate a confirmational space of protein binding pockets.

Is a free, open source binding pocket analysis tool intended for use in trajectory analysis.

**[This is a very, very early version of POVME3.0. Target date for actual release is September 2016. Until then there is no expectation of support or correctness of results, and I recommend that only developers try to use this.]**

## Download

```bash
git clone
cd POVME
```

You will also need [VMD](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD) and [pip](link)

## Example

```bash
cd examples
cd basic_example
python ../../POVME2.py sample_input.ini
```

<!-- TODO explain what is happening in the example -->

Once this runs, you will have an output directory named `basic_example_output`.

Within this directory, unzip `basic_example_volume_trajectory.pdb.gz` and open with [VMD](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD).

This file is a volumetric trajectory, which is a series of `.pdb` frames representing the shape of the binding pocket in each snapshot. If you load this along with the original snapshots, you will be able to visualize the binding pocket trajectory.

## Using VMD

Open VMD

File, new molecule

Browse, basic_example, output, volumetric_trajectory.pdb

load

load files for new molecule

browse, basic_example, `4NSS.pdb`

(`4NSS.pdb` is protein input file, but pocket volumes dont make sense unless you see them overlayed on the protein.)

load

Press play on bottom right to see it go.

### Additional

You probably want to turn the speed down

Show the protein as a "NewCartoon"

- graphics > representations > Drawing Method > NewCartoon
- Apply

Show binding pocket as spheres

- graphics > representations
- select `basic_example_volume_trajectory.pdb` from selected molecule
- drawing method > vdw
- apply  



## Dependencies

`Peel`, heavily based on [Binana](jacob durant), enables us to color the binding site shape with relevant chemical features, such as hydrogen bond donors/acceptors, potential pi-stacking interactions, and hydrophobic pockets. This is a standalone library that will eventually be broken out into its own package. `Peel` contains the `Algebra` class, which enables comparisons of and mathematical operations on binding pocket shapes.

`Clustering` is a package that ensembles of binding site shapes and perform clustering and principal component analysis on them. Clustering can be used to find metastable binding site shapes; principal component analysis can be used to find correlated motions in binding pockets.

[`Pocket ID`](jacob durant) is a library that runs successive, course iterations of POVME analysis on a whole protein to find binding pockets.

[`PyMolecule`](jacob durant) (soon to be replaced with something more reasonable, probably MDTraj or Prody) is a class to read `.pdb` files, which is the standard format for protein structures.

## Usage

## Examples

## Contribute

## Todo

- [ ]
