import POVME.packages.binana.peel as peel
import POVME.packages.pymolecule.pymolecule as pymolecule
#import numpy
my_params = peel.defaultParams
my_protein = pymolecule.Molecule()
my_protein.fileio.load_pdb_into('2gfc_reduced.pdb', serial_reindex=True, resseq_reindex=False)

my_peel = peel.peel(my_protein, my_params)
my_peel.write_vmd_script('visualize_2gfc.vmd', peel.defaultParams)

#povmeMap = numpy.ones((13,13,13))

#my_peel.color_povme_map(povmeMap, [0,13,0,13,0,13], 1)
# If given no "features" argument, create_feature_maps will create maps of all features
my_feature_maps = my_peel.create_feature_maps([10,80,-30,45,-30,30], 1)

my_feature_maps['hbondAcceptor'].write_pdb('HBA.pdb')
my_feature_maps['hbondAcceptor'].write_dx_file('HBA.dx')
my_feature_maps['hbondDonor'].write_pdb('HBD.pdb')
my_feature_maps['hbondDonor'].write_dx_file('HBD.dx')
my_feature_maps['aromatic'].write_pdb('ARO.pdb')
my_feature_maps['aromatic'].write_dx_file('ARO.dx')
my_feature_maps['hydrophobic'].write_pdb('HPB.pdb')
my_feature_maps['hydrophobic'].write_dx_file('HPB.dx')
my_feature_maps['hydrophilic'].write_pdb('HPL.pdb')
my_feature_maps['hydrophilic'].write_dx_file('HPL.dx')
#my_feature_maps['hydrophobicity'].write_pdb('HPBTY.pdb')
#my_feature_maps['hydrophobicity'].write_dx_file('HPBTY.dx')


print("Done!")
