# POVME

POVME (Pocket VOlume MEasurer) is a Python package for extracting actionable information from ensembles of protein structures for use in drug design.

Medicinal chemists might use this to find representative confirmations of a binding pocket that encompass the diversity found in a simulation or ensemble of structures.

POVME is a way to compress a large amount of protein confirmations into clusters or families

POVME is a tool to map and navigate a confirmational space of protein binding pockets.

Is a free, open source binding pocket analysis tool intended for use in trajectory analysis.

**[This is a very, very early version of POVME3.0. Target date for actual release is September 2016. Until then there is no expectation of support or correctness of results, and I recommend that only developers try to use this.]**

## Installation

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
