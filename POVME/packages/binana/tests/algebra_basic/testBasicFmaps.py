import POVME.packages.binana.peel as peel
import POVME.packages.pymolecule.pymolecule as pymolecule
import numpy
features = ['aromatic','hbondDonor','hbondAcceptor','hydrophobic','hydrophilic','hydrophobicity','occupancy']

mostly_empty = pymolecule.Molecule()
mostly_empty.fileio.load_pdb_into('mostly_empty.pdb')

my_peel = peel.peel(mostly_empty, peel.defaultParams)

one_point_fmaps_1 = my_peel.create_feature_maps([0,0,0,0,0,0], 1, features)
one_point_fmaps_2 = my_peel.create_feature_maps([0,0,0,0,0,0], 1, features)

print()
print('======Beginning scoring ======')
for i in range(len(features)):
    #So the point of featuremap_1 for aromatic will be set to 0, hBondDonor 2, etc...
    one_point_fmaps_1[features[i]].setData(numpy.array([[[2*i]]]))
    #So the point of featuremap_2 for aromatic will be set to 1, hBondDonor 3, etc...
    one_point_fmaps_2[features[i]].setData(numpy.array([[[(2*i)+1]]]))
    print('Setting %s to be %i for featuremap_A and %i for featuremap_B' %(features[i], i*2, (i*2)+1))


my_algebra = peel.algebra()
my_score_funcs = ['hbondAcceptor_A * hbondDonor_B', 'hbondAcceptor_B * hbondDonor_A']
my_algebra.setScoreFuncs(my_score_funcs)
my_score_maps, my_scores = my_algebra.scoreAll(one_point_fmaps_1, one_point_fmaps_2)
print('\n'.join(['%.5e\t%s'%(score, func) for func, score in zip(my_score_funcs, my_scores)]))



print()
print('Now trying out meaner score functions')
print() 

mean_score_funcs = ['-2 * hydrophobic_A * hydrophobic_B',
                    'hbondAcceptor_A * numpy.sqrt(hbondDonor_B)',
                    'numpy.sqrt(hbondDonor_A) * hbondAcceptor_B',
                    '-aromatic_A - aromatic_B',
                    '(hbondAcceptor_A + hbondDonor_B) * (aromatic_A + aromatic_B)',
                    '(hbondAcceptor_B + hbondDonor_A) * (aromatic_A + aromatic_B)',
                    '0.0001 * pow(10, hydrophobicity_A * hydrophobicity_B )']



my_algebra.setScoreFuncs(mean_score_funcs)
my_score_maps, my_scores = my_algebra.scoreAll(one_point_fmaps_1, one_point_fmaps_2)
print('\n'.join(['%.5e\t%s'%(score, func) for func, score in zip(mean_score_funcs, my_scores)]))

print() 
print('Now trying out previously buggy score functions')
print()

buggy_score_funcs = ['hbondAcceptor_A * hbondDonor_B', 'hbondAcceptor_B * hbondDonor_A','-(hbondDonor_A * hbondDonor_B)']

my_algebra.setScoreFuncs(buggy_score_funcs)
my_score_maps, my_scores = my_algebra.scoreAll(one_point_fmaps_1, one_point_fmaps_2)
print('\n'.join(['%.5e\t%s'%(score, func) for func, score in zip(buggy_score_funcs, my_scores)]))
