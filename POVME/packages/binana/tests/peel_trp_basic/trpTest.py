import POVME.packages.binana.peel as peel
import POVME.packages.pymolecule.pymolecule as pymolecule
my_params = peel.defaultParams
trp = pymolecule.Molecule()
trp.fileio.load_pdb_into( 'trp.pdb', bonds_by_distance=True, serial_reindex = True, resseq_reindex=False)
my_peel = peel.peel(trp, my_params)
my_peel.write_vmd_script('visualize_trp.vmd', peel.defaultParams)


my_feature_maps = my_peel.create_feature_maps([-20,20,-20,20,-20,20], 0.5)
my_feature_maps['hbondAcceptor'].write_pdb('HBA.pdb')
my_feature_maps['hbondDonor'].write_pdb('HBD.pdb')
my_feature_maps['aromatic'].write_pdb('ARO.pdb')

print("Done!")
