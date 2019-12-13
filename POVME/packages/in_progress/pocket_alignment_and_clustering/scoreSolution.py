import numpy as np
import sys
import pylab
import random
import copy

if not len(sys.argv) == 3:
    print('usage: python scoreSolution.py attemptFile keyFile')
    sys.exit()






def bootstrap_analysis(list1, list2):
    nBootstraps = 3000
    
    realSimilarity = test_solution(list1, list2)
    list1 = np.array(list1)
    list2 = np.array(list2)
    #list1Copy = copy.deepcopy(list1)
    #list2Copy = copy.deepcopy(list2)
    similarities = []
    for i in range(nBootstraps):
        if i % 1000 == 0:
            print('bootstrap',i)
        newIndices = np.random.randint(0,len(list1),(len(list1)))
        list1Scramble = list1[newIndices]
        #list1Scramble = list(list1)
        #random.shuffle(list1Scramble)
        newIndices = np.random.randint(0,len(list1),(len(list1)))
        list2Scramble = list2[newIndices]
        #print list1Copy
        #random.shuffle(pvl)
        similarity = test_solution(list1Scramble,list2)
        #print similarity
        similarities.append(similarity)
    hist = pylab.hist(similarities, bins=50)
    #realSimilarity = sum(realDifferences)

    betterPop = 0
    for pop,border in zip(hist[0], hist[1][1:]):
        if border < realSimilarity:
            betterPop += pop
    
    pylab.axvline(realSimilarity)
    betterThanPercent = 100*float(betterPop)/nBootstraps
    pylab.text(realSimilarity, 0.5*max(hist[0]),'%.1f' %(betterThanPercent))
    print('Real similarity is better than %r percent of bootstrap permutations' %(betterThanPercent))






def test_solution(groupList1, groupList2, plot = False):
    ASMat = np.zeros((len(groupList1),len(groupList1)),dtype=np.bool)
    for iInd, i in enumerate(groupList1):
        for jInd, j in enumerate(groupList1):
            if i == j:
                ASMat[iInd, jInd] = 1

    #for index, row in enumerate(ASMat):
    #    if sum(row) == 0:
    #        print 'No cluster for entry in row %i' %(index)
    #        1/0

    if plot:
        pylab.subplot(121)
        pylab.imshow(ASMat, interpolation='nearest')

    realMat = np.zeros((len(groupList2),len(groupList2)),dtype=np.bool)
    for iInd, i in enumerate(groupList2):
        for jInd, j in enumerate(groupList2):
            if i == j:
                realMat[iInd, jInd] = 1

    #for index, row in enumerate(realMat):
    #    if sum(row) == 0:
    #        print 'No cluster for entry in row %i' %(index)
    #        1/0

    if plot:
        pylab.subplot(122)
        pylab.imshow(realMat, interpolation='nearest')
        pylab.show()

    attNClusters = len(set(groupList1))
    realNClusters = len(set(groupList2))

    #print '%i clusters proposed, compared to %i clusters actually present' %(attNClusters, realNClusters)

    # Score = correctly predicted linkages (minus the diagonal) / the correct number of linkages
    simScore = float(np.sum(ASMat&realMat) - len(realSolution))# / np.sum(realMat)
    return simScore




    
attemptedSolution = open(sys.argv[1]).readlines()
attemptedSolutionGroups = np.array([int(i.split()[0]) for i in attemptedSolution])
realSolution = open(sys.argv[2]).readlines()
realSolutionGroups = np.array([int(i.split()[0]) for i in realSolution])

if len(attemptedSolution) != len(realSolution):
    print('error: Different number of classifications in the proposed solution and key file')

realScore = test_solution(attemptedSolutionGroups, realSolutionGroups)

print(realScore)

bootstrap_analysis(attemptedSolutionGroups, realSolutionGroups)
pylab.show()
#Doing bootstrap analysis given the number and size of clusters
