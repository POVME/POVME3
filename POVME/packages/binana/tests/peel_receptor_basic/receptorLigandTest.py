import POVME.packages.binana.peel as peel
import POVME.packages.pymolecule.pymolecule as pymolecule
#import packages.pymolecule.pymolecule as pymolecule
#import numpy
my_params = peel.defaultParams
#my_protein = pymolecule.Molecule()
ligand = peel.PDB()
ligand.LoadPDB('4FSM-results_50221585.pdb')
my_ligand_peel = peel.peel(ligand, my_params, isLigand=True)
my_ligand_peel.write_vmd_script('visualize_lig.vmd', peel.defaultParams)

protein = pymolecule.Molecule()
protein.fileio.load_pdb_into('4FSM_leaped.pdb', serial_reindex=True, resseq_reindex=False)
#protein = peel.PDB()
#protein.LoadPDB('4FSM_leaped.pdb')
my_protein_peel = peel.peel(protein, my_params, isLigand=False)
my_protein_peel.write_vmd_script('visualize_prot.vmd', peel.defaultParams)

#povmeMap = numpy.ones((13,13,13))

#my_peel.color_povme_map(povmeMap, [0,13,0,13,0,13], 1)
my_ligand_feature_maps = my_ligand_peel.create_feature_maps([-5,30,-20,5,-5,20], .5)
my_protein_feature_maps = my_protein_peel.create_feature_maps([-5,30,-20,5,-5,20], .5)

my_protein_feature_maps['occupancy'].write_dx_file('PROT_OCC.dx')
my_protein_feature_maps['hbondDonor'].write_dx_file('PROT_HBD.dx')
my_protein_feature_maps['hbondAcceptor'].write_dx_file('PROT_HBA.dx')
my_protein_feature_maps['aromatic'].write_dx_file('PROT_ARO.dx')
my_protein_feature_maps['hydrophobicity'].write_dx_file('PROT_HPBTY.dx')

my_ligand_feature_maps['occupancy'].write_dx_file('LIG_OCC.dx')
my_ligand_feature_maps['hbondDonor'].write_dx_file('LIG_HBD.dx')
my_ligand_feature_maps['hbondAcceptor'].write_dx_file('LIG_HBA.dx')
my_ligand_feature_maps['aromatic'].write_dx_file('LIG_ARO.dx')
my_ligand_feature_maps['hydrophobicity'].write_dx_file('LIG_HPBTY.dx')

print("Done!")
