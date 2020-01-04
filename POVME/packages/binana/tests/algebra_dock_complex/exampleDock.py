
import itertools
import POVME.packages.binana.peel as peel
import numpy
import sys
import os
import scipy.ndimage.interpolation as sni
import pickle
import csv


profile = False

if profile:
    from line_profiler import LineProfiler

    def do_profile(follow=[]):
        def inner(func):
            def profiled_func(*args, **kwargs):
                try:
                    profiler = LineProfiler()
                    profiler.add_function(func)
                    for f in follow:
                        profiler.add_function(f)
                    profiler.enable_by_count()
                    return func(*args, **kwargs)
                finally:
                    profiler.print_stats()
            return profiled_func
        return inner






features = ['aromatic','hbondAcceptor','hbondDonor','hydrophilic','hydrophobic','hydrophobicity','occupancy','adjacency']

#receptorName = sys.argv[1]
#ligandName = sys.argv[2]
receptorName = "2QO8"
ligandName = "2QO8-results_11032.pdb"

receptorData = {}
#ligandData = {}
for feature in features:
    #thisReceptorName = 'proteinFM/%s%s_.npy' %(receptorName, feature)
    thisReceptorName = '%s%s_.npy' %(receptorName, feature)
    thisReceptorData = numpy.load(thisReceptorName)
    receptorData[feature] = peel.featureMap.fromPovmeList(thisReceptorData, 1.)


#with open('ligandFM/%s_featureMaps.cPickle' %(ligandName)) as my_file:
#    ligandData = cPickle.load(my_file)

my_params = peel.defaultParams

ligand_molecule = peel.PDB()
ligand_molecule.LoadPDB(ligandName)
ligand_peel = peel.peel(ligand_molecule, my_params, isLigand=True)
my_algebra = peel.algebra()

#scoreMaps, scores = my_algebra.scoreAll(receptorData, ligandData)



number = str(5)
nsteps = str(500)

outputFile = number + "_" + nsteps
#os.system('/lv_scratch/j5wagner/projects/POVME_docking_inputs/points_on_sphere/4d_points_on_sphere.o ' + number + " " + nsteps + " > " + outputFile + '.out')
spherePoints = numpy.genfromtxt(outputFile + '.out',skip_header=3, skip_footer=0, usecols=(4,5,6,7), 
comments='}')
#spherePoints = [[0,1,0,0]]
print(spherePoints)



#n_cubes = 50
#ligandDim = ligandData['aromatic'].getShape()
#receptorDim = receptorData['aromatic'].shape
#receptorSize = receptorData['aromatic'].getData().size
#print ligandDim, receptorDim, receptorSize
#cubeLength = (receptorSize / n_cubes) ** .33333333
#print cubeLength
#x = []
#y = []
#z = []
#for n in numpy.arange(0,receptorDim[0]-ligandDim[0],cubeLength):
#	x.append(n)
#for n in numpy.arange(0,receptorDim[1]-ligandDim[1],cubeLength):
#	y.append(n)
#for n in numpy.arange(0,receptorDim[2]-ligandDim[2],cubeLength):
#	z.append(n)
translation_list = list(itertools.product(list(range(-8,9,2)),list(range(-8,9,2)),list(range(-8,9,2))))
#print translation_list
my_score_functions = ['hydrophobic_A * hydrophobic_B',
                     'hydrophilic_A * hydrophilic_B', 
                     'hydrophilic_A * hydrophobic_B', 
                     'hydrophobic_A * hydrophilic_B', 
                     'hbondDonor_A * hbondAcceptor_B',
                     'hbondAcceptor_A * hbondDonor_B',
                     'hbondAcceptor_A * hbondAcceptor_B',
                     'hbondDonor_A * hbondDonor_B',
                     'aromatic_A * aromatic_B',
                     #'aromatic_A * hbondDonor_A * hbondAcceptor_B',
                     #'aromatic_A * hbondAcceptor_A * hbondDonor_B',
                     #'aromatic_B * hbondDonor_A * hbondAcceptor_B',
                     #'aromatic_B * hbondAcceptor_A * hbondDonor_B',
                     #'hydrophobicity_A * hydrophobicity_B',
                     'occupancy_A * occupancy_B',
                     'adjacency_A * occupancy_B',
                     'occupancy_A * adjacency_B',
                     'adjacency_A * adjacency_B']
my_algebra.setScoreFuncs(my_score_functions)

#a = my_algebra.dockOne(receptorData, ligandData, translation_list, spherePoints)
if profile:
    @do_profile(follow=[peel.algebra.dockPeel, peel.algebra.scoreAll, peel.algebra.scoreOne, peel.algebra.commonizeVolume, peel.featureMap.snap_to_grid, peel.featureMap.point_to_nearest_index, peel.featureMap.index_to_coord])
    def do_docking():
        result = my_algebra.dockPeel(receptorData, ligand_peel, translation_list, spherePoints)
        return result

else:
    def do_docking():
        result = my_algebra.dockPeel(receptorData, ligand_peel, translation_list, spherePoints)
        return result

a=do_docking()

with open('results_%s.cPickle' %(ligandName), 'w') as fo:
    pickle.dump(a, fo)



'''import packages.binana.peel as peel
import numpy
import sys
import cPickle

features = ['aromatic','hbondAcceptor','hbondDonor','hydrophilic','hydrophobic','hydrophobicity','occupancy']

receptorName = sys.argv[1]
ligandName = sys.argv[2]

receptorData = {}
#ligandData = {}
for feature in features:
    thisReceptorName = 'proteinFM/%s%s_.npy' %(receptorName, feature)
    thisReceptorData = numpy.load(thisReceptorName)
    receptorData[feature] = peel.featureMap.fromPovmeList(thisReceptorData, 1.)


with open('ligandFM/%s_featureMaps.cPickle' %(ligandName)) as my_file:
    ligandData = cPickle.load(my_file)

my_algebra = peel.algebra()

scoreMaps, scores = my_algebra.scoreAll(receptorData, ligandData)

print my_algebra.getScoreFuncs()
print scores
'''
