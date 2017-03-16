import POVME.packages.pymolecule.pymolecule as pymolecule
import numpy
import random
import os
import shutil

# This program demonstrates all the features of pymolecule 2.0, a python framework
# for loading, saving, and manipulating 3D molecular data.

########## First, set up the example directory ##########
if os.path.exists('./example_output'): shutil.rmtree('./example_output')
os.mkdir('./example_output')

########## Saving and Loading Molecular Data ##########

# First, create a molecule object.
molecule = pymolecule.Molecule()
#molecule.fileio.load_pym_into('./example_output/pym_example.pym') # remove this later

# Now, load a PDB file into that object. Calculate the bonds by distance, and
# reindex the serial and resseq fields. Definitions used::
#   pymolecule.Molecule.fileio.load_pdb_into
#   pymolecule.Molecule.fileio.load_pdb_into_using_file_object
#   pymolecule.Molecule.fileio.serial_reindex
#   pymolecule.Molecule.fileio.resseq_reindex
#   pymolecule.Molecule.atoms_and_bonds.create_bonds_by_distance
#   pymolecule.Molecule.selections.select_all_atoms_bound_to_selection
molecule.fileio.load_pdb_into('./example.pdb', True, True, True)

# Save the molecule object as a PDB file. The serial and resseq fields have
# already been reindexed, so we'll skip that. Definitions used:
#   pymolecule.Molecule.fileio.save_pdb
molecule.fileio.save_pdb('./example_output/first_save.pdb', False, False, False)
# You may wish to verify that ./example_output/first_save.pdb has CONECT data and
# sequential serial/resseq fields

# You can also "save" to a string instead of to a file. Definitions used:
#   pymolecule.Molecule.fileio.save_pdb
print molecule.fileio.save_pdb('FILENAME DOES NOT MATTER', False, False, True)[:997]

# Pymolecule can also save the data in the pym format, which is generally
# faster than the PDB format. Definitions used::
#   pymolecule.Molecule.fileio.save_pym
molecule.fileio.save_pym('./example_output/pym_example.pym')

# You can also load from the pym format. Let's create a new molecule
# object and load the pym file into it. Definitions used::
#   pymolecule.Molecule.fileio.load_pym_into
mol_from_pym = pymolecule.Molecule()
mol_from_pym.fileio.load_pym_into('./example_output/pym_example.pym')

########## Atom Selections ##########

# Pymolecule has some excelent selection functions. This example shows
# how to select all the nitrogen and oxygen atoms in PRO residues.
# Definitions used::
#   pymolecule.Molecule.selections.select_atoms
PRO_C_N_selection = molecule.selections.select_atoms({'element_stripped':['C','N'], 'resname_stripped':'PRO'})

# You can also invert that selection. Definitions used:
#   pymolecule.Molecule.selections.invert_selection
not_PRO_C_N_selection = molecule.selections.invert_selection(PRO_C_N_selection)

# You can also get the selection of all the atoms near an existing
# selection (using a user-defined distance cutoff). For example,
# let's get all the atoms within 3 A of the PRO_C_N_selection above.
# Definitions used:
#   pymolecule.Molecule.selections.select_atoms_near_other_selection
near_PRO_C_N_selection = molecule.selections.select_atoms_near_other_selection(PRO_C_N_selection, 3.0)

# Note that near_PRO_C_N_selection doesn't include the atoms specified by
# PRO_C_N_selection. If you want them all together, you can use numpy's
# set-theory routines.
near_PRO_C_N_selection_including_PRO_C_N_selection = numpy.union1d(near_PRO_C_N_selection, PRO_C_N_selection)

# You can also get a selection that includes all the atoms that are in
# the same residues as the atoms of another selection. Definitions used:
#   pymolecule.Molecule.selections.select_atoms_in_same_residue
same_whole_residue_as_near_PRO_C_N_selection_including_PRO_C_N_selection = molecule.selections.select_atoms_in_same_residue(near_PRO_C_N_selection_including_PRO_C_N_selection)

# It's also possible to select all atoms within a bounding box. Definitions
# used:
#   pymolecule.Molecule.selections.select_atoms_in_bounding_box
bounding_box = numpy.array([[-5.0, 60.0, 85.0], [15.0, 80.0, 105.0]])
bounding_box_selection = molecule.selections.select_atoms_in_bounding_box(bounding_box)

# If your pymolecule.Molecule object actually contains multiple molecules,
# you can select all the atoms that belong to the same molecule as a
# user-specified atom. For example, let's select all the atoms that belong
# to the same molecule as the first atom in our pymolecule.Molecule object.
# (This should just be chain A). Note that this only works if you've
# specified the bonds. Definitions used:
#   pymolecule.Molecule.selections.select_atoms_from_same_molecule
same_molecule_as_first_atom_selection = molecule.selections.select_atoms_from_same_molecule([0])

# You can do the same with entire selections, instead of individual indecies.
# For example, let's select all the molecules containing either the first
# or the last atom. This should effectively select both chain A and B.
# Definitions used:
#   pymolecule.Molecule.selections.select_atoms_from_same_molecule
same_molecule_as_first_or_last_atom_selection = molecule.selections.select_atoms_from_same_molecule(numpy.array([0,len(molecule.information.get_coordinates())-1]))

# You can also get the selections of all the unique molecules. This can be
# used to divide a pymolecule.Molecule object containing multiple molecules
# (e.g., chains) into separate molecules. In this case, it will divide
# molecule into chain A and B, which are not connected by any bond and so
# are separate molecules. Note that, again, bonds must have been defined.
# Definitions used:
#   pymolecule.Molecule.selections.selections_of_constituent_molecules
unique_molecule_selections = molecule.selections.selections_of_constituent_molecules()

# You can get a seleciton of all atoms that are bonded to a given atom
# selection. For example, let's get a selection for all the atoms bound to
# the first and 26th atom (not including those atoms themselves). Definitions
# used:
#   pymolecule.Molecule.select_all_atoms_bound_to_selection
bonded_atoms_selection = molecule.selections.select_all_atoms_bound_to_selection(numpy.array([0,25]))

# You can select the atoms of a branch. You need only specify the root atom
# of the branch and the first atom of the branch attached to that root,
# to indicate the directionality. For example, the 3186th atom of molecule
# is a TRP alpha carbon, and the 3189th atom is the TRP beta carbon.
# so using the 3186th and 3189th atom as the root and directionality atom,
# respectively, you can select the alpha carbon + TRP sidechain.
# Definitions used:
#   pymolecule.Molecule.select_branch
TRP_side_chain_and_alpha_carbon = molecule.selections.select_branch(3186, 3189)

# Another common task is to select all atoms. Definitions used:
#   pymolecule.Molecule.selections.select_all
all_atoms = molecule.selections.select_all()

# Any selection can be converted into a molecule object. Definitions used:
#   pymolecule.Molecule.selections.get_molecule_from_selection
PRO_C_N_molecule = molecule.selections.get_molecule_from_selection(PRO_C_N_selection)
not_PRO_C_N_molecule = molecule.selections.get_molecule_from_selection(not_PRO_C_N_selection)
all_atoms_molecule = molecule.selections.get_molecule_from_selection(all_atoms)
near_PRO_C_N_molecule = molecule.selections.get_molecule_from_selection(near_PRO_C_N_selection)
near_PRO_C_N_selection_including_PRO_C_N_molecule = molecule.selections.get_molecule_from_selection(near_PRO_C_N_selection_including_PRO_C_N_selection)
same_whole_residue_as_near_PRO_C_N_selection_including_PRO_C_N_molecule = molecule.selections.get_molecule_from_selection(same_whole_residue_as_near_PRO_C_N_selection_including_PRO_C_N_selection)
same_molecule_as_first_atom_molecule = molecule.selections.get_molecule_from_selection(same_molecule_as_first_atom_selection)
same_molecule_as_first_or_last_atom_molecule = molecule.selections.get_molecule_from_selection(same_molecule_as_first_or_last_atom_selection)
bonded_atoms_molecule = molecule.selections.get_molecule_from_selection(bonded_atoms_selection)
TRP_side_chain_and_alpha_carbon_molecule = molecule.selections.get_molecule_from_selection(TRP_side_chain_and_alpha_carbon)
bounding_box_molecule = molecule.selections.get_molecule_from_selection(bounding_box_selection)
unique_molecules = [molecule.selections.get_molecule_from_selection(sel) for sel in unique_molecule_selections]

# Like any molecule object, these can also be saved to files. Note that here we
# are not reindexing the serial or resseq fields.
PRO_C_N_molecule.fileio.save_pdb('./example_output/selection_PRO_C_N.pdb', False, False)
not_PRO_C_N_molecule.fileio.save_pdb('./example_output/selection_not_PRO_C_N.pdb', False, False)
all_atoms_molecule.fileio.save_pdb('./example_output/selection_all_atoms.pdb', False, False)
near_PRO_C_N_molecule.fileio.save_pdb('./example_output/selection_near_PRO_C_N.pdb', False, False)
near_PRO_C_N_selection_including_PRO_C_N_molecule.fileio.save_pdb('./example_output/selection_near_PRO_C_N_selection_including_PRO_C_N.pdb', False, False)
same_whole_residue_as_near_PRO_C_N_selection_including_PRO_C_N_molecule.fileio.save_pdb('./example_output/selection_same_whole_residue_as_near_PRO_C_N_selection_including_PRO_C_N.pdb', False, False)
same_molecule_as_first_atom_molecule.fileio.save_pdb('./example_output/selection_same_molecule_as_first_atom.pdb', False, False)
same_molecule_as_first_or_last_atom_molecule.fileio.save_pdb('./example_output/selection_same_molecule_as_first_or_last_atom.pdb', False, False)
bonded_atoms_molecule.fileio.save_pdb('./example_output/selection_bonded_atoms.pdb', False, False)
TRP_side_chain_and_alpha_carbon_molecule.fileio.save_pdb('./example_output/selection_TRP_side_chain_and_alpha_carbon.pdb', False, False)
bounding_box_molecule.fileio.save_pdb('./example_output/bounding_box.pdb', False, False)
for index, mol in enumerate(unique_molecules): mol.fileio.save_pdb('./example_output/unique_mols_' + str(index + 1) + '.pdb', False, False)

# There's also a shortcut for making a copy of an existing molecule object,
# so you don't have to create a selection and then a molecule from that selection.
# Definitions used::
#   pymolecule.Molecule.copy
mol_copy = PRO_C_N_molecule.copy()

# Again, like any molecule object, this copy can also be saved.
mol_copy.fileio.save_pdb('./example_output/selection_PRO_C_N_copy.pdb')

# Another common task is to separate a molecule by chain or residue. Pymolecule
# provides special functions to get all these selections. First, I'll demonstrate
# how to get all the chain selections. Definitions used:
#   pymolecule.Molecule.selections.get_chain_divisions
chain_selections = molecule.selections.selections_of_chains()

# chain_selections is a dictionary, where the keys are the chainids and the values
# are the corresponding selections. So, to save each chain to a separate file:
molecule_A = molecule.selections.get_molecule_from_selection(chain_selections['A'])
molecule_A.fileio.save_pdb("./example_output/2HU4_chain_A.pdb")
molecule_B = molecule.selections.get_molecule_from_selection(chain_selections['B'])
molecule_B.fileio.save_pdb("./example_output/2HU4_chain_B.pdb")

# Similarly, you can get the selections of all the residues. Definitions used:
#   pymolecule.Molecule.selections.selections_of_residues
residue_selections = molecule.selections.selections_of_residues()

# residue_selections is a dictionary, where they keys have the form
# {resname}-{resseq}-{chainid} (e.g., THR-461-B) and the values are the corresponding
# selections. So, to save the first five residues to separate files:
for index, residueid in enumerate(residue_selections.keys()):
    if index > 4: break
    residue = molecule.selections.get_molecule_from_selection(residue_selections[residueid])
    residue.fileio.save_pdb('./example_output/2HU4_residue_' + str(index+1) + ".pdb")

# You can also select atoms from two different molecules that come within a user-
# specified distance of each other. Let's identify atoms from chain A and B that come
# within 10 A of each other using select_close_atoms_from_different_molecules_atom_by_atom.
# If the third parameter is set to True, this function performs a simple pair-wise
# distance comparison. It's fairly fast, but will not work well when comparing two
# molecules that are very large. Definitions used:
#   pymolecule.Molecule.selections.select_close_atoms_from_different_molecules
molecule_A_close_atoms_selection, molecule_B_close_atoms_selection = molecule_A.selections.select_close_atoms_from_different_molecules(molecule_B, 10.0, True)

molecule_A_close_atoms = molecule_A.selections.get_molecule_from_selection(molecule_A_close_atoms_selection)
molecule_B_close_atoms = molecule_B.selections.get_molecule_from_selection(molecule_B_close_atoms_selection)

molecule_A_close_atoms.fileio.save_pdb("./example_output/chain_A_close_atoms.pdb")
molecule_B_close_atoms.fileio.save_pdb("./example_output/chain_B_close_atoms.pdb")

# If the third parameter is set to False, a more sophisticated method for identifying
# juxtaposed atoms is used. This method encompasses the whole molecule, each chain,
# and each residue in a bounding sphere and then checks the spheres for overlap to
# minimize the number of expensive pair-wise distance comparisons required. The first
# time it's run, it's slower than when a pairwise comparison is prformed. However, the
# location and radii of the bounding spheres are saved, so that it's much faster the
# second time it's run. This function is also ideal when trying to detect steric
# clashes between two very large molecules. Definitions used:
#   pymolecule.Molecule.selections.select_close_atoms_from_different_molecules
molecule_A_close_atoms_selection_v2, molecule_B_close_atoms_selection_v2 = molecule_A.selections.select_close_atoms_from_different_molecules(molecule_B, 10.0, False)

molecule_A_close_atoms_v2 = molecule_A.selections.get_molecule_from_selection(molecule_A_close_atoms_selection_v2)
molecule_B_close_atoms_v2 = molecule_B.selections.get_molecule_from_selection(molecule_B_close_atoms_selection_v2)

molecule_A_close_atoms_v2.fileio.save_pdb("./example_output/chain_A_close_atoms_v2.pdb")
molecule_B_close_atoms_v2.fileio.save_pdb("./example_output/chain_B_close_atoms_v2.pdb")

########## Manipulating Atomic Coordinates ##########

# PyMolecule also has some useful functions for manipulating molecular data (i.e., 
# rotating and translating molecules, etc.). In order to undo any of these
# manipulations, we must first set the undo point. Definitions used:
#   pymolecule.Molecule.manipulation.set_coordinate_undo_point
molecule.manipulation.set_coordinate_undo_point()

# Now, let's translate the molecule (without rotation) so that the 5th atom is
# at (1.0, 1.0, 1.0) and save it. Definition used:
#   pymolecule.Molecule.manipulation.set_atom_location
molecule.manipulation.set_atom_location(4, numpy.array([1.0, 1.0, 1.0]))
molecule.fileio.save_pdb('./example_output/2HU4_5th_atom_to_pt.pdb')

# We can also move the molecule 2.0 A up the x axis. Definitions used:
#   pymolecule.Molecule.manipulation.translate_molecule
molecule.manipulation.translate_molecule(numpy.array([2.0, 0.0, 0.0]))
molecule.fileio.save_pdb('./example_output/2HU4_translate_2_down_x.pdb')

# There are also functions for rotating functions. For example, here's how
# to rotate a molecule 45 degress around a line connecting the points (1,1,1)
# and (10,0,10). Definitions used:
#   pymolecule.Molecule.manipulation.rotate_molecule_around_a_line_between_points
molecule.manipulation.rotate_molecule_around_a_line_between_points(numpy.array([[1.0, 1.0, 1.0]]), numpy.array([[10.0, 10.0, 10.0]]), numpy.pi*0.25)
molecule.fileio.save_pdb('./example_output/2HU4_rotation_test_1.pdb')

# You can also rotate the molecule 45 degrees around the line formed by
# connecting the 4th and 5th atom. Definitions used:
#   pymolecule.Molecule.manipulation.rotate_molecule_around_a_line_between_atoms
molecule.manipulation.rotate_molecule_around_a_line_between_atoms(3, 4, numpy.pi*0.25)
molecule.fileio.save_pdb('./example_output/2HU4_rotation_test_2_about_4th_5th_atom.pdb')

# The molecule can also be pivoted around a point, in this case (10.0, 10.0, 10.0)
# Definitions used:
#   pymolecule.Molecule.manipulation.rotate_molecule_around_pivot_point
molecule.manipulation.rotate_molecule_around_pivot_point(numpy.array([10.0,10.0,10.0]), numpy.pi*0.25, numpy.pi*0.25, numpy.pi*0.25)
molecule.fileio.save_pdb('./example_output/2HU4_rotation_test_about_pivot_pt.pdb')

# Or the molecule can be pivoted around the coordinates of a given atom, in this case
# the 5th atom. Definitions used:
#   pymolecule.Molecule.manipulation.rotate_molecule_around_pivot_atom
molecule.manipulation.rotate_molecule_around_pivot_atom(4, numpy.pi*0.25, numpy.pi*0.25, numpy.pi*0.25)
molecule.fileio.save_pdb('./example_output/2HU4_rotation_test_about_pivot_5th_atom.pdb')

########## Obtaining Molecular Information ##########

# There are also PyMolecule functions for getting information about the molecule
# like the center of mass, the geometric center, the total mass, etc.
# Definitions used:
#   pymolecule.Molecule.information.get_center_of_mass
#   pymolecule.Molecule.information.get_geometric_center
#   pymolecule.Molecule.information.get_total_mass
#   pymolecule.Molecule.information.get_total_number_of_atoms
#   pymolecule.Molecule.information.get_total_number_of_heavy_atoms
#   pymolecule.Molecule.information.get_bounding_box
print "Center of mass:", molecule.information.get_center_of_mass()
print "Geometric center:", molecule.information.get_geometric_center()
print "Total mass:", molecule.information.get_total_mass()
print "Total number of atoms:", molecule.information.get_total_number_of_atoms()
print "Total number of heavy atoms:", molecule.information.get_total_number_of_heavy_atoms()
print "Bounding box:", molecule.information.get_bounding_box()

# You can also get some information about what atoms are bonded to selected atoms.
# For example, get the number of hydrogen atoms bonded to the 1st atom:
# Definitions used:
#   pymolecule.Molecule.atoms_and_bonds.get_number_of_bond_partners_of_element
print "Number of hydrogen atoms bound to first atom:", molecule.atoms_and_bonds.get_number_of_bond_partners_of_element(0, "H")

# You can also get the index of the first atom with a given element bound to
# a specified atom. Definitions used:
#   pymolecule.Molecule.atoms_and_bonds.get_index_of_first_bond_partner_of_element
print "Index of the first hydrogen atom bound to first atom:", molecule.atoms_and_bonds.get_index_of_first_bond_partner_of_element(0, "H")

# Functions also exist for determining whether a given atom belongs to a protein,
# DNA, or RNA molecule. Consider the 3943rd atom, which is known to belong to a
# protein. Definitions used:
#   pymolecule.Molecule.information.belongs_to_protein
#   pymolecule.Molecule.information.belongs_to_dna
#   pymolecule.Molecule.information.belongs_to_rna
print "3943rd atom belongs to protein?", molecule.information.belongs_to_protein(3942)
print "3943rd atom belongs to DNA?", molecule.information.belongs_to_dna(3942)
print "3943rd atom belongs to RNA?", molecule.information.belongs_to_rna(3942)

########## Constructing Molecules Programatically ##########

# It's also possible to construct molecules programatically, rather than
# loading them from files. For example, here's how to create a new
# molecule object and add five carbon atoms. Definitions used:
#   pymolecule.Molecule.atoms_and_bonds.add_atom
new_mol = pymolecule.Molecule()
new_mol.atoms_and_bonds.add_atom(serial=1, element="C", coordinates = numpy.array([0,0,0]))
new_mol.atoms_and_bonds.add_atom(serial=2, element="C", coordinates = numpy.array([5,5,0]))
new_mol.atoms_and_bonds.add_atom(serial=3, element="C", coordinates = numpy.array([10,0,0]))
new_mol.atoms_and_bonds.add_atom(serial=4, element="C", coordinates = numpy.array([10,10,0]))
new_mol.atoms_and_bonds.add_atom(serial=5, element="C", coordinates = numpy.array([0,10,0]))

# You can also add bonds programatically. Definitions used:
#   pymolecule.Molecule.atoms_and_bonds.add_bond
new_mol.atoms_and_bonds.add_bond(0,2)
new_mol.atoms_and_bonds.add_bond(2,3)
new_mol.atoms_and_bonds.add_bond(3,4)
new_mol.atoms_and_bonds.add_bond(4,0)

# Let's save the file so you can examine it, taking note that the bonds
# are present in the PDB CONECT entries.
new_mol.fileio.save_pdb('./example_output/from_scratch_v1.pdb',False,False)

# There are also functions for deleting atoms and bonds. Definitions used:
#   pymolecule.Molecule.atoms_and_bonds.delete_atom
#   pymolecule.Molecule.atoms_and_bonds.delete_bond
new_mol.atoms_and_bonds.delete_atom(1)
new_mol.atoms_and_bonds.delete_bond(2,1)

# Let's save that PDB as well so you can compare them:
new_mol.fileio.save_pdb('./example_output/from_scratch_v2.pdb',False,False)

########## Comparing Multiple Molecules ##########

# PyMolecule also provides useful functions for comparing multiple molecule objects.
# Let's see if there are any steric clashes between chain A and B, using a large
# 10 A cutoff. As above, if the third parameter is True, a simple pair-wise distance
# comparison is performed. Definitions used:
#   pymolecule.Molecule.other_molecule.steric_clash_with_another_molecule
print "Clash according to steric_clash_with_another_molecule_atom_by_atom:", molecule_A.other_molecule.steric_clash_with_another_molecule(molecule_B, 10.0, True)

# If the third parameter is False (as above), the more sophisticated method for
# identifying juxtaposed atoms is used. Definitions used:
#   pymolecule.Molecule.other_molecule.steric_clash_with_another_molecule
print "Clash according to steric_clash_with_another_molecule:", molecule_A.other_molecule.steric_clash_with_another_molecule(molecule_B, 10.0, False)

# PyMolecule can also be used to identify the minimum distance between two molecules.
# The second parameter is like the third above, determining whether pair-wise or more
# sophisticated methods are used to identify distances. Definitions used:
#   pymolecule.Molecule.other_molecule.get_distance_to_another_molecule
print 'Closest inter-chain distance according to pairwise comparison:', molecule_A.other_molecule.get_distance_to_another_molecule(molecule_B, True)
print 'Closest inter-chain distance according to sophisticated method:', molecule_A.other_molecule.get_distance_to_another_molecule(molecule_B, False)
#### I THINK THE ABOVE, THE BETTER METHOD MIGHT TAKE LONGER. WHY??? CHECK GENERALLY.

# PyMolecule can calculate RMSD values. To demonstrate this, let's look at the first
# 100 atoms of the two chains.
molecule_A_first_hundred = molecule_A.selections.get_molecule_from_selection(numpy.arange(0,100,1,dtype=int))
molecule_B_first_hundred = molecule_B.selections.get_molecule_from_selection(numpy.arange(0,100,1,dtype=int))

# The atoms of these two molecules are in the same order, so we can determine which
# atoms are equivalent based on order alone. Definitions used:
#   pymolecule.Molecule.other_molecule.get_rmsd_order_dependent
print "RMSD between two chains, based on order: ", molecule_A_first_hundred.other_molecule.get_rmsd_order_dependent(molecule_B_first_hundred)

# You can also specify which atoms are equivalent. This is accomplished using a tuple
# of two numpy.array, where each array is comprised of the indecies of self and
# other_mol, respectively, such that equivalent atoms are listed in order. So, for
# example, if (atom1/mol1 = atom6/mol2) and (atom2/mol1 = atom3/mol2) than the tethers
# would be ([1,2], [6,3]). Let's just generate these tethers randomly to demonstrate
# how this functionality is used. Definitions used:
#   pymolecule.Molecule.other_molecule.get_rmsd_equivalent_atoms_specified
mol1_index_in_order = numpy.arange(0,len(molecule_A_first_hundred.information.get_coordinates()),1,dtype=int)
mol2_index_in_order = numpy.arange(0,len(molecule_B_first_hundred.information.get_coordinates()),1,dtype=int)
numpy.random.shuffle(mol1_index_in_order)
numpy.random.shuffle(mol2_index_in_order)

print "RMSD between two chains, equiv atoms specified:", molecule_A_first_hundred.other_molecule.get_rmsd_equivalent_atoms_specified(molecule_B_first_hundred, (mol1_index_in_order, mol2_index_in_order))

##### YOU NEED ORDER INDEPENDENT ONE TOO!!!! WITH FINGERPRINTS???? ######

# PyMolecule also has a function to merge two molecules. Let's recreate the original
# molecule object, with chains A and B in the same object. Note the importance in this
# case of reindexing the serial field when you save, or the bonds may be incorrect.
# Definitions used:
#   pymolecule.Molecule.other_molecule.merge_with_another_molecule
merged = molecule_A.other_molecule.merge_with_another_molecule(molecule_B)
merged.fileio.save_pdb('./example_output/merged.pdb',True,False)

# You can also align portions of molecules to minimize the RMSD distance between them.
# For example, let's align the two chains by the first TRP residue in chain A and the
# last in chain B. Definitions used:
#   pymolecule.Molecule.other_molecule.get_other_molecule_aligned_to_this

molecule_A_first_TRP_selection = molecule_A.selections.select_atoms_in_same_residue(molecule_A.selections.select_atoms({'resname_stripped':'TRP'})[0])
molecule_B_first_TRP_selection = molecule_B.selections.select_atoms_in_same_residue(molecule_B.selections.select_atoms({'resname_stripped':'TRP'})[-1])
tethers = numpy.vstack((molecule_A_first_TRP_selection, molecule_B_first_TRP_selection))
molecule_B_aligned_to_A = molecule_A.other_molecule.get_other_molecule_aligned_to_this(molecule_B, tethers)
molecule_B_aligned_to_A.fileio.save_pdb('./example_output/last_TRP_of_chain_B_aligned_to_first_TRP_of_chain_A.pdb', False, False, False)

########## Random Geometry Functions ##########

# Most geometry functions should be done through numpy. However, there are a
# few functions that have dedicated functions, since they aren't necessarily
# trivial to do in numpy, or because of historical reasons (i.e., they were
# included in older versions of pymolecule)

# For example, you can get the angle between three points. Definitions used:
#   pymolecule.Molecule.geometry.get_angle_between_three_points
print "Angle in radians between 1st, 2nd, and 3rd atoms:", molecule.geometry.get_angle_between_three_points(molecule.information.get_coordinates()[0], molecule.information.get_coordinates()[1], molecule.information.get_coordinates()[2])

# You can also get the dihedral angle between four points. Definitions used:
#   pymolecule.Molecule.geometry.get_dihedral_angle
print "Dihedral angle in radians between 1st, 2nd, 3rd, and 4th atoms:", molecule.geometry.get_dihedral_angle(molecule.information.get_coordinates()[0], molecule.information.get_coordinates()[1], molecule.information.get_coordinates()[2], molecule.information.get_coordinates()[3])

# There's also a function to determine whether or not four points are in
# a plane (within a given tolerance). For example, the 202nd, 203rd, 207th,
# and 208th atoms in molecule belong to the same TRP sidechain and so are
# roughly planar. Definitions used:
#   pymolecule.Molecule.geometry.is_planar
#   pymolecule.Molecule.geometry.get_planarity_deviation
#four_atoms_from_TRP_sidechain_selection = numpy.array([201, 202, 206, 207])
print "Four atoms in the same TRP sidechain are roughly planar?", molecule.geometry.is_planar(molecule.information.get_coordinates()[201], molecule.information.get_coordinates()[202], molecule.information.get_coordinates()[206], molecule.information.get_coordinates()[207], 0.2)
