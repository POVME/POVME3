import POVME.packages.binana.peel as peel
import POVME.packages.pymolecule.pymolecule as pymolecule


# The purpose of this test case is to show the versatility of the string-input scoring functions.
# Users can use functions from the numpy package to construct score components.
# Please note that these functions are meant to show the versatility of the software
# and do not have any actual meaning.

part_1 = pymolecule.Molecule()
part_1.fileio.load_pdb_into('helixPt1.pdb')

part_2 = pymolecule.Molecule()
part_2.fileio.load_pdb_into('helixPt2.pdb')

#peel_1.write_vmd_script_file('peel1.vmd'')

peel_1 = peel.peel(part_1, peel.defaultParams)
fmaps_1 = peel_1.create_feature_maps([25,55,-25,5,-15,15], 1, ['occupancy','hbondDonor','hbondAcceptor','aromatic','hydrophobic','hydrophilic','hydrophobicity'])
peel_2 = peel.peel(part_2, peel.defaultParams)
fmaps_2 = peel_2.create_feature_maps([25,55,-25,5,-15,15], 1, ['occupancy','hbondDonor','hbondAcceptor','aromatic','hydrophobic','hydrophilic','hydrophobicity'])





my_algebra = peel.algebra()
print()
print()
print('The scoring function components are originally:')
print('\n'.join(my_algebra.getScoreFuncs()))
print()
print("However, I'd prefer them to be")
my_score_funcs = ['-2 * hydrophobic_A * hydrophobic_B',
                  'hbondAcceptor_A * numpy.sqrt(hbondDonor_B)',
                  'numpy.sqrt(hbondDonor_A) * hbondAcceptor_B',
                  '-aromatic_A - aromatic_B',
                  '0.0001 * numpy.exp( hydrophobicity_A * hydrophobicity_B )']
print('\n'.join(my_score_funcs))


my_algebra.setScoreFuncs(my_score_funcs)

my_score_maps, my_scores = my_algebra.scoreAll(fmaps_1, fmaps_2)

print('The scores for each term are')
print('\n'.join(['%.5e\t%s'%(score, func) for func, score in zip(my_score_funcs, my_scores)]))

my_score_maps[0].write_dx_file('term1.dx')
my_score_maps[1].write_dx_file('term2.dx')
my_score_maps[2].write_dx_file('term3.dx')
my_score_maps[3].write_dx_file('term4.dx')
my_score_maps[4].write_dx_file('term5.dx')
