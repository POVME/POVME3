import POVME.packages.binana.peel as peel
import numpy as np
import pylab
import sys
import random

randomFactor = 0.001
plot = False
createPdbs = True
createNpys = True
mainPocket = [peel.point([0,0,0]),7]
sidePocket1 = [peel.point([0,0,6]),4]
sidePocket2 = [peel.point([0,6,0]),4]
sidePocket3 = [peel.point([6,0,0]),4]
sidePocket4 = [peel.point([0,4,4]),4]

families = {}
# families[1] = [21,mainPocket,sidePocket1]
# families[2] = [21,mainPocket,sidePocket2]
# families[3] = [21,mainPocket,sidePocket3]
# families[4] = [20,mainPocket,sidePocket4]
# families[5] = [20,mainPocket]
# families[6] = [20,mainPocket,sidePocket1,sidePocket2]
# families[7] = [20,mainPocket,sidePocket1,sidePocket4]
# families[8] = [20,mainPocket,sidePocket1,sidePocket2,sidePocket4]
# families[9] = [20,mainPocket,sidePocket1,sidePocket2,sidePocket3,sidePocket4]
families[1] = [10,mainPocket,sidePocket1]
families[2] = [31,mainPocket,sidePocket2]
#families[3] = [7,mainPocket,sidePocket3]
#families[4] = [46,mainPocket,sidePocket4]
#families[5] = [17,mainPocket]
#families[6] = [25,mainPocket,sidePocket1,sidePocket2]
#families[7] = [2,mainPocket,sidePocket1,sidePocket4]
#families[8] = [18,mainPocket,sidePocket1,sidePocket2,sidePocket4]
#families[9] = [27,mainPocket,sidePocket1,sidePocket2,sidePocket3,sidePocket4]

familyMembership = []
#for index, family in enumerate(families.keys()):
for family in list(families.keys()):
    for replicate in range(families[family][0]):
        outPrefix = 'family%i_rep%02i' %(family, replicate)
        this_fm = peel.featureMap([-20,20,-20,20,-20,20],1)
        for pocket in families[family][1:]:
            randomShift = peel.point(np.random.normal(0,randomFactor,(3)))
            this_fm.add_sphere(pocket[0].point_sum_new(randomShift),
                               pocket[1])
        familyMembership.append((family, outPrefix, this_fm))
        if createPdbs == True: 
            this_fm.write_pdb(outPrefix+'.pdb')
        if createNpys == True:
            np.save(outPrefix+'.npy', this_fm.toPovmeList()[:,:3])


familyMap = np.zeros((len(familyMembership),len(familyMembership)))
for iInd, i in enumerate(familyMembership):
    for jInd, j in enumerate(familyMembership):
        if i[0] == j[0]:
            familyMap[iInd, jInd] = 1

if plot:
    pylab.imshow(familyMap, interpolation='nearest')
    pylab.show()

familyMembershipShuffled = [i for i in familyMembership]
random.shuffle(familyMembershipShuffled)

familyMap = np.zeros((len(familyMembershipShuffled),len(familyMembershipShuffled)))
for iInd, i in enumerate(familyMembershipShuffled):
    if createNpys == True:
        np.save('shuffled%03i.npy' %(iInd),i[2].toPovmeList()[:,:3])
    if createPdbs == True:
        i[2].write_pdb('shuffled%03i.pdb' %(iInd))
    for jInd, j in enumerate(familyMembershipShuffled):
        if i[0] == j[0]:
            familyMap[iInd, jInd] = 1
            
if plot:
    pylab.imshow(familyMap, interpolation='nearest')
    pylab.show()

import time

with open('key_%s' %(time.strftime('%y_%m_%d_%H_%M_%S')),'wb') as of:
    of.write('\n'.join(['%r %r' %(i[0],i[1]) for i in familyMembershipShuffled]))
              
