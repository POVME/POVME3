import POVME.packages.binana.peel as peel
import POVME.packages.pymolecule.pymolecule as pymolecule

part_1 = pymolecule.Molecule()
part_1.fileio.load_pdb_into('helixPt1.pdb')

part_2 = pymolecule.Molecule()
part_2.fileio.load_pdb_into('helixPt2.pdb')

#peel_1.write_vmd_script_file('peel1.vmd'')

peel_1 = peel.peel(part_1, peel.defaultParams)
fmaps_1 = peel_1.create_feature_maps([25,55,-25,5,-15,15], 1, ['hbondDonor','hbondAcceptor'])
peel_2 = peel.peel(part_2, peel.defaultParams)
fmaps_2 = peel_2.create_feature_maps([25,55,-25,5,-15,15], 1, ['hbondDonor','hbondAcceptor'])

#peel_1.write_vmd_script_file('peel_1.vmd')
#peel_2.write_vmd_script_file('peel_2.vmd')

fmaps_1['hbondDonor'].write_dx_file('hbd1.dx')
fmaps_1['hbondAcceptor'].write_dx_file('hba1.dx')
fmaps_2['hbondDonor'].write_dx_file('hbd2.dx')
fmaps_2['hbondAcceptor'].write_dx_file('hba2.dx')

my_algebra = peel.algebra()
my_algebra.setScoreFuncs(['hbondAcceptor_A * hbondDonor_B',
                          'hbondAcceptor_B * hbondDonor_A'])
my_score_maps, my_scores = my_algebra.scoreAll(fmaps_1, fmaps_2)

my_score_maps[0].write_dx_file('term1.dx')
my_score_maps[1].write_dx_file('term2.dx')
