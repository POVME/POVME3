# POVME

We present a substantial update to the open-source POVME binding pocket analysis software. New capabilities of POVME 3.0 include a flexible chemical coloring scheme for feature identification, post-analysis tools for comparing large ensembles of pockets (e.g., from molecular dynamics simulations), and the introduction of scripts and methods that facilitate binding pocket comparison and analysis. We envision the use of this software for visualization of binding pocket dynamics, selection of representative structures for ensemble docking, and incorporation of molecular dynamics results into ligand design efforts.

## Install

If POVME2 is installed on your system, we recommend making a separate python environment. POVME3 has been tested on builds using miniconda. For our testing, we use:

```wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p miniconda2
source miniconda2/bin/activate
pip install povme
```

Note that this method will create a separate python build. You will need to run "source miniconda2/bin/activate" each time you want to use POVME3.



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

## Dependencies

`Peel`, heavily based on [Binana](jacob durant), enables us to color the binding site shape with relevant chemical features, such as hydrogen bond donors/acceptors, potential pi-stacking interactions, and hydrophobic pockets. This is a standalone library that will eventually be broken out into its own package. `Peel` contains the `Algebra` class, which enables comparisons of and mathematical operations on binding pocket shapes.

`Clustering` is a package that ensembles of binding site shapes and perform clustering and principal component analysis on them. Clustering can be used to find metastable binding site shapes; principal component analysis can be used to find correlated subpockets in binding sites.

[`Pocket ID`](jacob durant) is a library that runs successive, coarse iterations of POVME analysis on a whole protein to find binding pockets.

[`PyMolecule`](jacob durant) (soon to be replaced with something more reasonable, probably MDTraj or Prody) is a class to read `.pdb` files, which is the standard format for protein structures.

## Usage

## Examples

## Contribute

## Todo

- [ ]
