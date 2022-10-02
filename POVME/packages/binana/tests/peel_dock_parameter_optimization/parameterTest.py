# This is a meta-script made with the goal of quantitatively assessing decisions made in the creation of feature fields and the POVME docking algorithm.
# It takes as inputs a POVME/PEEL arun script, a folder of protein PDBs, and a folder of ligand PDBs that are in the right location to fit into their respective protein PDBs.
# It then runs peel to classify the shape and feature maps of the protein pockets, then tries docking the ligand PDBs into the proper protein feature map, with a certain distance between each sampled translation and a set number of rotations.
# Then, it analyzes the "goodness" of the dockings.
# Since we don't know how the different interactions should be weighted given the changing natures of the color fields, it tries to optimize the coefficients attached to each term in such a way as to
# place the rank the correct placement of the ligand first among all sampled dockings (0 translation from the original spot and minimum rotation).
# The script then produces a scatter plot of the results, with the x-axis showing distance from the correct position and the y-axis showing score.
# It also produces a plot of the steric clash scores for each ligand, with the x-axis representing the different docked conformations and the y-axis showing their steric clash score.
# It returns the sum of the ranks of the correct placements. The best combination of parameters would minimize the sum of the correct ranks (for instance, the correct placements should all be in 1st place given a perfect scoring function).

import copy
import pickle
import os
import os.path
import glob
import gzip
import itertools
import numpy
import parameterTestHelpers
import POVME.packages.binana.peel as peel
import pylab
import POVME.packages.pymolecule.pymolecule as pymolecule
import scipy.optimize
import sys

#arunScript = '/home/j5wagner/git/AmaroLabCodeDevel20140815/arun'
proteinPdbDir = '/lv_scratch/j5wagner/projects/POVME_docking_inputs/proteinPdbs'
ligandPdbDir = '/lv_scratch/j5wagner/projects/POVME_docking_inputs/ligandPdbs'


#prefix = 'occScale1.0_occStyleHard_yesNorm_50rots_gridRes1.0'
prefix = sys.argv[1]
nRotations = sys.argv[2]
gridResolution = float(sys.argv[3])
outputDir = '/lv_scratch/j5wagner/projects/POVME_docking_inputs/parameterTests/' + prefix
#proteinFmScriptDir = outputDir + '/proteinFM'
proteinFmOutputDir = outputDir + '/proteinFM'
dockingResultsOutputDir = outputDir + '/dockingResults'
analysisOutputDir = outputDir + '/analysis'

if not os.path.exists(outputDir):
    os.mkdir(outputDir)
    os.mkdir(proteinFmOutputDir)
    os.mkdir(dockingResultsOutputDir)
    os.mkdir(analysisOutputDir)


docking_rotations_list_file = '/lv_scratch/j5wagner/projects/POVME_docking_inputs/rotationsCodebook/%s.out' %nRotations
docking_rotation_list =  numpy.genfromtxt(docking_rotations_list_file, skip_header=3, skip_footer=0, usecols=(4,5,6,7), comments='}')
docking_translation_list = list(itertools.product(numpy.arange(-8./gridResolution,8.001/gridResolution,2./gridResolution),
                                                  numpy.arange(-8./gridResolution,8.001/gridResolution,2./gridResolution),
                                                  numpy.arange(-8./gridResolution,8.001/gridResolution,2./gridResolution)))

print(docking_translation_list[0])
my_score_functions = ['hydrophobic_A * hydrophobic_B',
                     'hydrophilic_A * hydrophilic_B',
                     'hydrophilic_A * hydrophobic_B',
                     'hydrophobic_A * hydrophilic_B',
                     'hbondDonor_A * hbondAcceptor_B',
                     'hbondAcceptor_A * hbondDonor_B',
                     'hbondAcceptor_A * hbondAcceptor_B',
                     'hbondDonor_A * hbondDonor_B',
                     'aromatic_A * aromatic_B',
                     'occupancy_A * occupancy_B',
                     'adjacency_A * occupancy_B',
                     'occupancy_A * adjacency_B',
                     'adjacency_A * adjacency_B']
my_algebra = peel.algebra()
my_algebra.setScoreFuncs(my_score_functions)

my_params = peel.defaultParams

# Figure out which all proteins we're working with

ligandPdbs = glob.glob(ligandPdbDir+'/*pdb')
proteinPdbPrefixes = list(set([i.split('/')[-1].split('-')[0] for i in ligandPdbs]))

# First we generate the protein maps.

allCoords = {}
centerCoords = {}
features = ['aromatic','hbondAcceptor','hbondDonor','hydrophilic','hydrophobic','occupancy', 'adjacency']
for proteinPdbPrefix in proteinPdbPrefixes:
    receptorFmStartedFlagFilename = '%s/started_%s' %(proteinFmOutputDir, proteinPdbPrefix)
    if not os.path.exists(receptorFmStartedFlagFilename):
        os.system('touch %s' %(receptorFmStartedFlagFilename))
        pdbFileName = '%s/%s_leaped.pdb' %(proteinPdbDir, proteinPdbPrefix)
        # We find a coordinate that represents the middle of the binding pocket
        allCoords[proteinPdbPrefix]=[]
        for ligandPdb in ligandPdbs:
            if proteinPdbPrefix in ligandPdb:
                pdbObj = peel.PDB()
                try:
                    pdbObj.LoadPDB(ligandPdb)
                except:
                    os.system('mv %s %s/failed' %(ligandPdb, ligandPdbDir))
                    continue
                for atom in pdbObj.AllAtoms:
                    allCoords[proteinPdbPrefix].append(pdbObj.AllAtoms[atom].coordinates.coords())

        allCoords[proteinPdbPrefix] = numpy.array(allCoords[proteinPdbPrefix])
        centerCoords[proteinPdbPrefix] = [numpy.mean(allCoords[proteinPdbPrefix][:,0]),
                                          numpy.mean(allCoords[proteinPdbPrefix][:,1]),
                                          numpy.mean(allCoords[proteinPdbPrefix][:,2])]
        # Generate the PDB map
        sideLength = 30
        protein_boundaries = [centerCoords[proteinPdbPrefix][0] - sideLength,
                              centerCoords[proteinPdbPrefix][0] + sideLength,
                              centerCoords[proteinPdbPrefix][1] - sideLength,
                              centerCoords[proteinPdbPrefix][1] + sideLength,
                              centerCoords[proteinPdbPrefix][2] - sideLength,
                              centerCoords[proteinPdbPrefix][2] + sideLength]
        protein_molecule = pymolecule.Molecule()
        protein_molecule.fileio.load_pdb_into(pdbFileName, serial_reindex=True, resseq_reindex=False)
        protein_peel = peel.peel(protein_molecule, my_params)
        protein_feature_maps = protein_peel.create_feature_maps(protein_boundaries,
                                                                gridResolution,
                                                                features = features)
        # Now save the featureMaps for this protein
        #with gzip.open('%s/%s_featureMaps.cPickle.gz' %(proteinFmOutputDir, proteinPdbPrefix),'wb') as my_file:
        my_file = gzip.open('%s/%s_featureMaps.cPickle.gz' %(proteinFmOutputDir, proteinPdbPrefix),'wb')
        pickle.dump(protein_feature_maps, my_file, protocol=2)
        my_file.close()
    else:
        #If the job is marked as started, let's see if it's finished. If so, load that receptor map
        if os.path.exists('%s/%s_featureMaps.cPickle.gz' %(proteinFmOutputDir, proteinPdbPrefix)):
            #with gzip.open('%s/%s_featureMaps.cPickle.gz' %(proteinFmOutputDir, proteinPdbPrefix),'rb') as my_file:
            my_file = gzip.open('%s/%s_featureMaps.cPickle.gz' %(proteinFmOutputDir, proteinPdbPrefix),'rb')
            protein_feature_maps = pickle.load(my_file)
            my_file.close()
        #Otherwise, there's probably another process working on making the feature map for this protein.
        #Skip it and go onto the next
        else:
            continue

    # Now that we've got the protein feature maps, let's start docking the ligands!
    for ligandPdb in ligandPdbs:
        ligandPrefix = ligandPdb.split('/')[-1]
        if proteinPdbPrefix in ligandPrefix:
            dockingStartedFlagFilename = '%s/started_%s' %(dockingResultsOutputDir, ligandPrefix)
            if not os.path.exists(dockingStartedFlagFilename):
                os.system('touch %s' %(dockingStartedFlagFilename))
                ligand_molecule = peel.PDB()
                ligand_molecule.LoadPDB(ligandPdb)
                ligand_peel = peel.peel(ligand_molecule, my_params, isLigand=True)

                results = my_algebra.dockPeel(protein_feature_maps, ligand_peel, docking_translation_list, docking_rotation_list)

                #with gzip.open('%s/results_%s.cPickle.gz' %(dockingResultsOutputDir, ligandPrefix), 'wb') as my_file:
                my_file = gzip.open('%s/results_%s.cPickle.gz' %(dockingResultsOutputDir, ligandPrefix), 'wb')
                pickle.dump(results, my_file, protocol=2)
                my_file.close()
            #else:
            #    #with gzip.open('%s/results_%s.cPickle.gz' %(dockingResultsOutputDir, ligandPrefix), 'rb') as my_file:
            #    my_file = gzip.open('%s/results_%s.cPickle.gz' %(dockingResultsOutputDir, ligandPrefix), 'rb')
            #    results = cPickle.load(my_file)
            #    my_file.close()
    #analysisStartedFlagFilename = '%s/started_%s' %(analysisOutputDir, proteinPdbPrefix)

#Analyze scores and optimize terms, then save figures showing how the proper pose rated in terms of combined score and occupancy

def scoreAll(factors):
    this_scoreListDict = copy.deepcopy(template_scoreListDict)
    for ligand in ligands:
        #print len(rotations[ligand]), len(translations[ligand])
        this_scoreListDict[ligand][:,7] = scoreTermsDict[ligand].dot(factors)
        this_scoreListDict[ligand] = this_scoreListDict[ligand][this_scoreListDict[ligand][:,7].argsort()][::-1]
    return this_scoreListDict


def rankFactors(factors):
    #This function takes a list of weighting factors and returns the rank that 0,0,0 gets when using them as multipliers on the score terms
    print(factors)
    this_scoreListDict = scoreAll(factors)
    print('scored!')
    ranks = []
    for ligand in ligands:
        for i, row in enumerate(this_scoreListDict[ligand]):
            #print row
            if row[0] == 0. and row[1] == 0. and row[2] == 0.:
                if row[3] == minRotation[ligand][0] and row[4] == minRotation[ligand][1] and row[5] == minRotation[ligand][2] and row[6] == minRotation[ligand][3]:
                    #print row
                    ranks.append(i)
                    break
    print(sum(ranks))
    return sum(ranks)

ligands = []
allScoreList = []
analysisStartedFlagFilename = '%s/started_total_analysis' %(analysisOutputDir)
if not os.path.exists(analysisStartedFlagFilename):
    os.system('touch %s' %(analysisStartedFlagFilename))
    allScoreListDict = {}
    for filename in glob.glob('%s/results*' %(analysisOutputDir)):
        print(filename)
        filePrefix = filename.replace(analysisOutputDir,'').replace('.cPickle','')
        ligands.append(filePrefix)
        #with gzip.open(filename, 'rb') as my_file:
        my_file = gzip.open(filename,'rb')
        scoreList = pickle.load(my_file)
        my_file.close()

        allScoreListDict[filePrefix] = copy.deepcopy(scoreList)
        allScoreList = allScoreList + scoreList

    minRotation = sorted(list(docking_rotation_list), key=lambda x:abs(x[0]))[0]


    ##Get normalization factors
    numTerms = len(allScoreList[0][2])
    means = [0.0 for i in range(numTerms)]
    for score in allScoreList:
        for i in range(numTerms):
            means[i] += abs(score[2][i]) / len(allScoreList)
    means = numpy.array(means)
    #Avoid dividing by 0
    means[means==0.0] = 1.0


    #This is to make sure that the plot has a reasonable aspect ratio. Otherwise large ligands will stretch the score scale
    ligandMeans = {}
    for ligand in ligands:
        ligandMeans[ligand] = [0.0 for i in range(numTerms)]
        for score in allScoreListDict[ligand]:
            for i in range(numTerms):
                #print abs(score[2][i]) / len(allScoreListDict[ligand])
                ligandMeans[ligand][i] += abs(score[2][i]) / len(allScoreListDict[ligand])
        ligandMeans[ligand] = numpy.array(ligandMeans[ligand])
        #Avoid dividing by 0
        ligandMeans[ligand][ligandMeans[ligand] == 0.0] = 1.


    ###This comes in handy when scoring later
    template_scoreListDict = {}
    scoreTermsDict = {}
    for ligand in ligands:
        template_scoreListDict[ligand] = numpy.zeros((len(docking_rotation_list)*len(docking_translation_list), 8))
        template_scoreListDict[ligand][:,0:3] = numpy.array([[i[0][0], i[0][1], i[0][2]] for i in allScoreListDict[ligand]])
        template_scoreListDict[ligand][:,3:7] = numpy.array([[i[1][0], i[1][1], i[1][2], i[1][3]] for i in allScoreListDict[ligand]])
        scoreTermsDict[ligand] =  numpy.array([i[2] for i in allScoreListDict[ligand]]) / ligandMeans[ligand]


    user_factors = [0., 0., -0., -0., 0.,
                    0., 0., 0., 0., -0.,
                    0., 0., 0.]


    my_optimization = scipy.optimize.minimize(rankFactors, #scoreFactors,
                                              user_factors,
                                              #method = 'COBYLA',
                                              options={'disp':True,
                                                       #'factor':0.1},
                                                       #'m':100,
                                                       'eps': 2},
                                              #'rhobeg':1},
                                              #'pgtol':1},

                                          bounds = ((0,20),(0,20),(-10,10),
                                                    (-10,10),(5,5),(5,5),
                                                    (-10,10),(-10,10),(0,10),
                                                    (-20,-10),(-20,20),(-20,20),
                                                    (-20,20))
                                          )

    print('OPTIMIZATION STATS:', my_optimization)
    optimizedFactors = my_optimization.x
    #with gzip.open('%s/optimization_results.cPickle.gz') as my_file:
    my_file = gzip.open('%s/optimization_results.cPickle.gz', 'wb')
    pickle.dump(my_optimization, my_file, protocol=2)
    my_file.close()



    finalScoreListDict = scoreAll(optimizedFactors)
    finalScoreList = []
    for key in finalScoreListDict:
        finalScoreList.extend(finalScoreListDict[key])
    distScoreList = numpy.array([[numpy.sqrt((i[0]**2)+(i[1]**2)+(i[2]**2)) + i[3], i[7]] for i in finalScoreList])

    pylab.scatter(distScoreList[:,0], distScoreList[:,1], color='black', s=3)

    #for ligand in ligands:
    for row in finalScoreList:
        if row[0] == 0 and row[1] == 0 and row[2] == 0 and tuple(row[3:7]) in list(minRotation.values()):
            pylab.scatter([0],[row[7]], color='r', s=100)

    h, xedges, yedges = numpy.histogram2d(distScoreList[:,0], distScoreList[:,1], bins=50)
    pylab.imshow(h.T, extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                         interpolation = 'nearest',
                         aspect = 'auto',
                         origin = 'lower')
    pylab.savefig('%s/optimization.pdf')



