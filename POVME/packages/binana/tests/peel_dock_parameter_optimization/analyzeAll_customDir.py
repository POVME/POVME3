import pickle
import sys
import numpy
import glob
import gzip
import copy
import sys
import threading
import queue
import time 
ligands = []

allScoreList = []
allScoreListDict = {}
#scoreDict = {}
#resultsDir = '../dockingResults/'
#resultsDir = '../dockingResults_20140816shapes_peelDock_fullRadius_normed/'
resultsDir = sys.argv[1]

def loadPickle(this_q, this_filename):
    filePrefix = this_filename.replace(resultsDir+'/','').replace('.cPickle','')
    fo = gzip.open(this_filename,'rb')
    #this_scoreList = cPickle.load(fo)
    data = fo.read()
    this_scoreList = pickle.loads(data)
    fo.close()
    this_q.put((filePrefix, this_scoreList))
    #this_scoreList2 = copy.deepcopy(this_scoreList)
    #this_q.put((filePrefix, this_scoreList, this_scoreList2))
    
q = queue.Queue()
threadList = []
for filename in glob.glob('%s/results*' %(resultsDir))[:30]:
#for filename in glob.glob('%s/results*' %(resultsDir))[:200]:
    print(filename)
    #filePrefix = filename.replace(resultsDir,'').replace('.cPickle','')
    #ligands.append(filePrefix)
    
    t = threading.Thread(target=loadPickle, args=(q,filename))
    t.daemon = True
    t.start()
    threadList.append(t)
    #fo = gzip.open(filename,'rb')
    #scoreList = cPickle.load(fo)
    #fo.close()

c = 0
while sum([i.is_alive() for i in threadList]) != 0:
    time.sleep(1)
    c += 1
    print('waiting for threads %i sec:' %(c), [i.is_alive() for i in threadList])
while not q.empty():
    prefix, scoreList = q.get()
    ligands.append(prefix)
    #print scoreList
    #scoreList = q.get()
    #allScoreListDict[prefix] = copy.deepcopy(scoreList)
    allScoreListDict[prefix] = scoreList
    allScoreList = allScoreList + scoreList




    #for entry in allScoreList:
    #    scoreDict[(tuple(entry[0]), tuple(entry[1]),filePrefix)] = entry[2]

#allScoreList = numpy.array(allScoreList)
        
rotations = {}
translations = {}
for ligand in ligands:
    rotations[ligand] = set([])
    translations[ligand] = set([])
    for row in allScoreListDict[ligand]:
        translations[ligand].add(tuple(row[0]))
        rotations[ligand].add(tuple(row[1]))
 
#for translation, rotation, ligand in scoreDict.keys():
#    translations[ligand].add(translation)
#    rotations[ligand].add(rotation)
#for key in scoreDict.keys():
#    translations.add(tuple(key[0]))
#    rotations.add(tuple(key[1]))

minRotation = {}
for ligand in ligands:
    #minRotation[ligand] = sorted(list(rotations[ligand]), key=lambda x:abs(x[0]))[0]
    minRotation[ligand] = sorted(list(rotations[ligand]), key=lambda x:abs(x[0]))[-1]
print(minRotation)


##Get normalization factors
numTerms = len(allScoreList[0][2])
means = [0.0 for i in range(numTerms)]
for score in allScoreList:
    for i in range(numTerms):
        means[i] += abs(score[2][i]) / len(allScoreList)
means = numpy.array(means)
#Avoid dividing by 0
means[means==0.0] = 1.

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
print(ligandMeans)
        
###This comes in handy when scoring later
template_scoreListDict = {}
scoreTermsDict = {}
for ligand in ligands:
    template_scoreListDict[ligand] = numpy.zeros((len(rotations[ligand])*len(translations[ligand]), 8))
    template_scoreListDict[ligand][:,0:3] = numpy.array([[i[0][0], i[0][1], i[0][2]] for i in allScoreListDict[ligand]])
    template_scoreListDict[ligand][:,3:7] = numpy.array([[i[1][0], i[1][1], i[1][2], i[1][3]] for i in allScoreListDict[ligand]])
    scoreTermsDict[ligand] =  numpy.array([i[2] for i in allScoreListDict[ligand]]) / ligandMeans[ligand]
   


#def computeScore(scoreList, factors):
#    return sum([scoreList[i]*factors[i]/means[i] for i in range(len(scoreList))])


'''
def getBestTranslations(factors):
    this_bestScoreTransDict= {}
    for translation in translations:
        transScores = []
        bestScore = -99999999
        bestRot = []
        for rotation in rotations:
            key = (translation, rotation)
            if computeScore(scoreDict[key], factors) > bestScore:
                bestRot = key[1]
                bestScore = computeScore(scoreDict[key], factors)
        this_bestScoreTransDict[(tuple(translation), tuple(bestRot))] = [scoreDict[(translation,bestRot)], factors]
        
    this_bestScoreTransList = numpy.array([[key[0][0],
                                       key[0][1],
                                       key[0][2],
                                       computeScore(scoreDict[key], factors)] 
                                      for key in this_bestScoreTransDict.keys()])
    this_bestScoreTransList = numpy.array(sorted(this_bestScoreTransList, key=lambda x: x[3]))
    return this_bestScoreTransList
'''


def scoreAll(factors):
    
#    this_scoreListDict = {}
#    for ligand in ligands:
#        this_scoreListDict[ligand] = numpy.zeros((len(translations[ligand])*len(rotations[ligand]), 8))
#        c=0
#        for translation in translations[ligand]:
#            for rotation in rotations[ligand]:
                
#                this_scoreListDict[ligand][c] = numpy.array([translation[0], translation[1], translation[2], rotation[0], rotation[1], rotation[2], rotation[3], computeScore(scoreDict[(translation, rotation, ligand)], factors)])
#                c += 1
#        this_scoreListDict[ligand] = numpy.array(sorted(this_scoreListDict[ligand], key=lambda x:x[2]))
    #this_scoreListDict = {}
    this_scoreListDict = copy.deepcopy(template_scoreListDict)
    for ligand in ligands:
        #print len(rotations[ligand]), len(translations[ligand])
        this_scoreListDict[ligand][:,7] = scoreTermsDict[ligand].dot(factors)
        #, axis = 0, order = (7))
        #this_scoreListDict[ligand] = numpy.array(sorted(this_scoreListDict[ligand], key=lambda x: x[7], reverse=True))
        #print this_scoreListDict[ligand][:3]
    
    return this_scoreListDict


def rankFactors(factors):
    #This function takes a list of weighting factors and returns the rank that 0,0,0 gets when using them as multipliers on the score terms
    #this_bestScoreTransList = getBestTranslations(factors)
    print(factors)
    this_scoreListDict = scoreAll(factors)
    print('scored!')
    ranks = []
    for ligand in ligands:
        if 1:
            a=this_scoreListDict[ligand]
            #tripleZRows = numpy.logical_and(a[:,0]==0
            tripleZScore = a[(a[:,0]==0) &
                             (a[:,1]==0) &
                             (a[:,2]==0) &
                             (a[:,3]==minRotation[ligand][0]) &
                             (a[:,4]==minRotation[ligand][1]) &
                             (a[:,5]==minRotation[ligand][2]) &
                             (a[:,6]==minRotation[ligand][3]),7]
            rank = numpy.count_nonzero(this_scoreListDict[ligand][:,7]>tripleZScore)            
            ranks.append(rank)
        else:
                
            #this_scoreListDict[ligand] = this_scoreListDict[ligand][this_scoreListDict[ligand][:,7].argsort()][::-1]
            for i, row in enumerate(this_scoreListDict[ligand]):
                #print row
                if row[0] == 0. and row[1] == 0. and row[2] == 0.:
                    if row[3] == minRotation[ligand][0] and row[4] == minRotation[ligand][1] and row[5] == minRotation[ligand][2] and row[6] == minRotation[ligand][3]:
                        #print row
                        #ranks.append(i)
                        tripleZScore = row[7]
                        rank = numpy.count_nonzero(this_scoreListDict[ligand][:,7]>tripleZScore)
                        ranks.append(rank)
                        break
    print(sum(ranks))
    return sum(ranks)
    #for i, row in enumerate(this_bestScoreTransList):
    #    if (row[:3] == numpy.array([0,0,0])).all():
    #        triple0 = i
    #return -triple0


'''
def scoreFactors(factors):
    #This function takes a list of weighting factors and returns a "score" that they get when using them as multipliers on the score terms
    this_bestScoreTransList = getBestTranslations(factors)
    runningScore = 0.0
    for i, row in enumerate(this_bestScoreTransList):
        #runningScore += (numpy.sqrt((row[0]**2) + (row[1]**2) + (row[2]**2))**8) * i
        runningScore += (0.03 + -numpy.exp(-(numpy.sqrt((row[0]**2) + (row[1]**2) + (row[2]**2))**2)/2)) * i
        #runningScore += (0.03 + -numpy.exp(-(numpy.sqrt((row[0]**2) + (row[1]**2) + (row[2]**2))**2)/5)) * row[3]
    return runningScore
'''    

def print_parameters(x):
    print(x)

import scipy.optimize


###Just as a reminder
#my_score_functions =['hydrophobic_A * hydrophobic_B',
#                 2    'hydrophilic_A * hydrophilic_B', 
#                 3    'hydrophilic_A * hydrophobic_B', 
#                 4    'hydrophobic_A * hydrophilic_B', 
#                 5    'hbondDonor_A * hbondAcceptor_B',
#                 6    'hbondAcceptor_A * hbondDonor_B',
#                 7    'hbondAcceptor_A * hbondAcceptor_B',
#                 8    'hbondDonor_A * hbondDonor_B',
#                 9    'aromatic_A * aromatic_B',
#                 10   'occupancy_A * occupancy_B',
#                 11   'adjacency_A * occupancy_B',
#                 12   'occupancy_A * adjacency_B',
#                 13   'adjacency_A * adjacency_B']

## Put in human-inputted factors

#user_factors = [5., 5., -5., -5., 5.,
#                5., 0., 0., 5., -15.,
#                5., 5., 0.]
user_factors = [0., 0., -0., -0., 0.,
                0., 0., 0., 0., -00.,
                0., 0., 0.]


my_optimization = scipy.optimize.minimize(rankFactors, #scoreFactors,
                                          user_factors, 
                                          #method = 'COBYLA',
                                          options={'disp':True,
                                                   #'factor':0.1},
                                                   #'m':100,
                                                   'eps': 2.},
                                          #'rhobeg':1},
                                          #        'pgtol':1},
                                          
                                          #bounds = ((0,20),(0,20),(-10,10),
                                          #          (-10,10),(5,5),(5,5),
                                          #          (-10,10),(-10,10),(0,10),
                                          #          (-100,-10),(-20,20),(-20,20),
                                          #          (-20,20))
                                          bounds = ((0,200),(0,200),(-100,100),
                                                    (-100,100),(0,100),(0,100),
                                                    (-100,100),(-100,100),(0,100),
                                                    (-500,-10),(-200,200),(-200,200),
                                                    (-200,200))
                                          )

print('OPTIMIZATION STATS:', my_optimization)
optimizedFactors = my_optimization.x

'''
my_optimization = scipy.optimize.fmin(rankFactors,user_factors, xtol = 0.01)
optimizedFactors = my_optimization
'''
'''
optimizedFactors = numpy.array([  6.27265804e-01,   4.54168056e+00,   5.46248338e+00,  -1.87883923e+00, 1.44989168e+00,   3.79782140e-01,   5.51136641e-04,  -2.47446270e-03, 1.14398554e+00,   2.64094296e-03,  -2.12394170e-03,  -2.48873856e-03, -9.11841469e-04,   7.13814468e-01,  -1.06753076e+00,   2.37556547e+00, 1.97917698e+00,  -2.62286481e+00])
'''
print('OPTIMIZATION FINISHED!!!!!!')
print(optimizedFactors)

#bestScoreTransList = getBestTranslations(my_optimization.x)

#import matplotlib as mpl
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt

#fig = plt.figure(1)
#fig.clf()
#ax = Axes3D(fig)
#print bestScoreTransList
#p=ax.scatter(bestScoreTransList[:,0],bestScoreTransList[:,1],bestScoreTransList[:,2],c=bestScoreTransList[:,3])
#fig.colorbar(p)
#plt.show()

import pylab

finalScoreListDict = scoreAll(optimizedFactors)
finalScoreList = []
for key in finalScoreListDict:
    finalScoreList.extend(finalScoreListDict[key])

#print finalScoreList
    
#distScoreList = numpy.array([[numpy.sqrt((i[0]**2)+(i[1]**2)+(i[2]**2)),i[7]] for i in finalScoreList])
#distScoreList = numpy.array([[numpy.sqrt((i[0]**2)+(i[1]**2)+(i[2]**2)) + (0.5*numpy.random.random())-0.25,i[7]] for i in finalScoreList])
distScoreList = numpy.array([[numpy.sqrt((i[0]**2)+(i[1]**2)+(i[2]**2)) + (1-abs(i[3])), i[7]] for i in finalScoreList])

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

pylab.savefig('test.pdf')        
pylab.show()
