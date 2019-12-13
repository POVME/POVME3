import POVME.packages.binana.peel as peel
#import packages.pymolecule.pymolecule as pymolecule
#import numpy
my_params = peel.defaultParams
#my_protein = pymolecule.Molecule()
ligand = peel.PDB()

ligand.LoadPDB('3CZ_reduced.pdb')
#my_protein.fileio.load_pdb_into('3CZ_reduced.pdb', bonds_by_distance=True, serial_reindex=True, resseq_reindex=False)


my_peel = peel.peel(ligand, my_params, isLigand=True)
my_peel.write_vmd_script('visualize_lig.vmd', peel.defaultParams)

#povmeMap = numpy.ones((13,13,13))

#my_peel.color_povme_map(povmeMap, [0,13,0,13,0,13], 1)
my_feature_maps = my_peel.create_feature_maps([50,70,-30,-10,45,65], 1)

my_feature_maps['occupancy'].write_dx_file('OCC.dx')
my_feature_maps['hbondAcceptor'].write_pdb('HBA.pdb')
#There are no hbond donor groups on the ligand
#my_feature_maps['hbondDonor'].write_pdb('HBD.pdb') 
my_feature_maps['aromatic'].write_pdb('ARO.pdb')
my_feature_maps['aromatic'].write_dx_file('ARO.dx')
my_feature_maps['hydrophobic'].write_pdb('HBC.pdb')
my_feature_maps['hydrophilic'].write_pdb('HPL.pdb')
#This is a bad feature
#my_feature_maps['hydrophobicity'].write_pdb('HPBTY.pdb')
#my_feature_maps['hydrophobicity'].write_dx_file('HPBTY.dx')

#HBAMap = my_feature_maps[0].data
#print HBAMap
#HBAPoints = numpy.transpose(numpy.nonzero(HBAMap))
#print HBAPoints

print("Done!")
