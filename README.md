# POVME

We present a substantial update to the open-source POVME binding pocket analysis software. New capabilities of POVME 3.0 include a flexible chemical coloring scheme for feature identification, post-analysis tools for comparing large ensembles of pockets (e.g., from molecular dynamics simulations), and the introduction of scripts and methods that facilitate binding pocket comparison and analysis. We envision the use of this software for visualization of binding pocket dynamics, selection of representative structures for ensemble docking, and incorporation of molecular dynamics results into ligand design efforts.

## Install

If POVME2 is installed on your system, we recommend making a separate python environment. POVME3 has been tested on builds using miniconda. For our testing, we use:

```wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p miniconda2
source miniconda2/bin/activate
pip install povme
```

Note that this method will create a separate python build. You will need to run ```source miniconda2/bin/activate``` each time you want to use POVME3.



## Example

The POVME Git repository comes with examples and test cases which are not included in the pip install.

```git clone https://github.com/POVME/POVME.git
cd POVME/POVME/examples/ligand_example/
POVME3.py sample_POVME_input.ini
```


Once this runs, you will have an output directory named `POVME_test_run`.

We recommend that you visualize the results using [VMD](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD). Open the POVME output using VMD with the following command: ```vmd -m POVME_test_run/POVME_volume_trajectory.pdb 1BYQ_every250.pdb```. Under the Graphics-->Representations menu in VMD, show the ```0: POVME_volume_trajectory.pdb``` molecule using the Drawing Method "VDW". Now press the play button in the bottom right corner of the VMD Main window to watch the pocket trajectory.

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
