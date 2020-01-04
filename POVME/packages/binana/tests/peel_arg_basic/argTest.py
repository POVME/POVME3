#import numpy
import POVME.packages.pymolecule.pymolecule as pymolecule
import POVME.packages.binana.peel as peel
import sys
#print sys.modules




my_params = peel.defaultParams

#DEFAULT
#arg = pymolecule.FileIO.load_pdb_into( 'arg.pdb', bonds_by_distance=True, serial_reindex = True, resseq_reindex=False)
arg = pymolecule.Molecule()
arg.fileio.load_pdb_into( 'arg.pdb', bonds_by_distance=True, serial_reindex = True, resseq_reindex=False)

my_peel = peel.peel(arg, my_params)
my_peel.write_vmd_script('visualize_arg.vmd', peel.defaultParams)


my_feature_maps = my_peel.create_feature_maps([-10,15,-5,15,-10,15], 1)
my_feature_maps['hbondAcceptor'].write_pdb('HBA.pdb')
my_feature_maps['hbondAcceptor'].write_dx_file('HBA.dx')
my_feature_maps['hbondDonor'].write_pdb('HBD.pdb')
my_feature_maps['hbondDonor'].write_dx_file('HBD.dx')
my_feature_maps['occupancy'].write_dx_file('OCC.dx')
#my_feature_maps['aromatic'].write_pdb('ARO.pdb') #don't bother with this
my_feature_maps['hydrophobic'].write_pdb('HBC.pdb')
my_feature_maps['hydrophilic'].write_pdb('HPL.pdb')
#my_feature_maps['hydrophobicity'].write_pdb('HPBTY.pdb')
#my_feature_maps['hydrophobicity'].write_dx_file('HPBTY.dx')
#HBAMap = my_feature_maps[0].data
#print HBAMap
#HBAPoints = numpy.transpose(numpy.nonzero(HBAMap))
#print HBAPoints

print("Done!")
