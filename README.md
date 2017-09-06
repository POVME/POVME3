# POVME

[![CircleCI](https://circleci.com/gh/POVME/POVME.svg?style=svg)](https://circleci.com/gh/POVME/POVME)

We present a substantial update to the open-source POVME binding pocket analysis software. New capabilities of POVME 3.0 include a flexible chemical coloring scheme for feature identification, post-analysis tools for comparing large ensembles of pockets (e.g., from molecular dynamics simulations), and the introduction of scripts and methods that facilitate binding pocket comparison and analysis. We envision the use of this software for visualization of binding pocket dynamics, selection of representative structures for ensemble docking, and incorporation of molecular dynamics results into ligand design efforts.

This document is specific to POVME 3.0. Users interested in POVME 2.0 will find resources [here](http://rocce-vm0.ucsd.edu/data/sw/hosted/POVME/)

## Install

If POVME2 is installed on your system, we recommend making a separate python environment. We normally install POVME3.0 using [miniconda](https://conda.io/docs/install/quick.html):

```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
bash Miniconda2-latest-Linux-x86_64.sh -b -p miniconda2
source miniconda2/bin/activate
pip install povme
```

Note that this method will create a separate python build. You will need to run ```source miniconda2/bin/activate``` in each terminal where you want to use POVME3.



# Examples Directory

The POVME Git repository comes with examples and test cases which are not included in the pip install.

```
git clone https://github.com/POVME/POVME.git
cd POVME/POVME/examples/
```

## Basic example
```
cd basic_example
POVME3.py sample_input.ini 
```

This example shows the "classic" operation of POVME, using a geometrically-defined inclusion sphere. If you open the "sample_input.ini" text file, you will find the operating parameters. The minimum input required for POVME to run is the input trajectory name and inclusion region.

Once this runs, you will have an output directory named `POVME_test_run`.

We recommend that you visualize the results using [VMD](http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD). Open the POVME output using VMD with the following command: 

```vmd -m POVME_test_run/POVME_volume_trajectory.pdb 1BYQ_every250.pdb```

Under the Graphics-->Representations menu in VMD, show the ```0: POVME_volume_trajectory.pdb``` molecule using the Drawing Method "VDW". Now press the play button in the bottom right corner of the VMD Main window to watch the pocket trajectory.


## Ligand-defined inclusion region example
```
cd ligand_example/
POVME3.py sample_POVME_input.ini
```

POVME 3.0 now allows users to define the inclusion region of a pocket using a ligand residue name. The pocket will then be defined in all grid points within 3 Angstroms of the ligand atoms in the loaded PDB trajectory. Note that this residue name must match the one given in the input PDB trajectory.


## Clustering and PCA example
```
cd analysis_workflow_example/
source runWorkflow.sh
```
The bulk of the new capabilities of POVME 3.0 are in separate scripts. Three of these are showcased in the analysis workflow example. 

This example runs POVME on 5 trajectories taken from the POVME 3.0 paper's HSP90 simulations. Each of these trajectory PDB files has 5 frames, and has had the ligand removed. After running POVME on these trajectories, three post-processing scripts are run:

* binding_site_overlap.py calculates the similarity of all of the analyzed frames 
* cluster.py processes the binding_site_overlap matrix and performs hierarchical clustering  
   * This example is programmed to yield five clusters, as specified in the "-n" argument to cluster.py  
   * A heatmap showing which frames belong to which cluster is displayed when cluster.py finishes running  
   * Combined and individual-simulation transition maps are displayed  
   * The most representative frames from each cluster are output in the ```3-post_analysis/ALL/cluster#``` subdirectories   
   * The average pocket shape of each cluster can be visualized in vmd by running ```vmd -e visualizeAll.vmd``` in the ```3-post_analysis/ALL``` subdirectory, and showing the second representation in each loaded object  
   * Text files of the cluster members and representatives are written, with each line corresponding to one cluster  
* pocketPointsPca.py runs principal component analysis of the analyzed frames  
   * Scatterplots of each simulations position in PC space are shown  
   * A plot of the explained variance for each PC is shown, as well as a line measuring the cumulative total  
   * The first 10 principal components can be visualized by running ```vmd -e loadAllPcs.vmd```  



## Miscellaneous

`Peel`, heavily based on Binana (by Jacob Durrant), enables the coloring of binding sites with relevant chemical features, such as hydrogen bond donors/acceptors, potential pi-stacking interactions, and hydrophobic pockets. This is a standalone library that will eventually be broken out into its own package. `Peel` contains the `Algebra` class, which enables comparisons of and mathematical operations on binding pocket shapes. The coloring scheme is intended to be for visualization only, and has not been validated for any quantitative purpose.

`Clustering` is a package that ensembles of binding site shapes and perform clustering and principal component analysis on them. Clustering can be used to find metastable binding site shapes; principal component analysis can be used to find correlated subpockets in binding sites.

`PyMolecule`(by Jacob Durrant) is a lightweight class to read PDB files.

