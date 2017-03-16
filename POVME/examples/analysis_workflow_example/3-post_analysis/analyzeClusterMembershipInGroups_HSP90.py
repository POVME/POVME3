import pylab
import numpy as np
import re

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import Draw
import os
import glob



def fixPdbSdfFormatting(data):
    datasp = data.split('\n')
    readingAtoms = True
    for lineInd, line in enumerate(datasp):
        linesp = line.split()
        if lineInd == 3:
            datasp[lineInd] = line[:16] + ' 0  0  0  0  0' +  line[30:]
        if (lineInd > 3):
            if len(linesp) != 16:
                readingAtoms = False
            if readingAtoms:
                atom = line[32:34].strip()
                datasp[lineInd] = line[:31] + atom.ljust(4) + '  ' + line[35:]
    return '\n'.join(datasp)



data = [i.split() for i in open('cluster_members.csv').readlines()]

print [len(i) for i in data]

nClusters = len(data)
nCols = 3
nRows = (nClusters+(nCols-1))/nCols
print nClusters, nRows, nCols

def getGroupPrefix(runPrefix):
    groupPrefix = runPrefix.split('run')[0].strip('_')
    return groupPrefix


groupPrefix2prefixes2frames = {}
prefix2frames = {}
for cluster in data:
    for member in cluster:
        result = re.findall('([a-zA-Z0-9_-]+)_([0-9]+)', member)
        #print member, result
        prefix = result[0][0].strip('_')
        groupPrefix = getGroupPrefix(prefix)
        if not groupPrefix in groupPrefix2prefixes2frames.keys():
            groupPrefix2prefixes2frames[groupPrefix] = {}
        groupPrefix2prefixes2frames[groupPrefix][prefix] = []
        frame = int(result[0][1])
        prefix2frames[prefix] = prefix2frames.get(prefix,[]) + [frame]
        groupPrefix2prefixes2frames[groupPrefix][prefix] += [frame]
#print prefix2frames

prefixes = prefix2frames.keys()
prefixes.sort()
nPrefixes = len(prefixes)
groupPrefixes = list(set([getGroupPrefix(i) for i in prefixes]))
groupPrefixes.sort()
nGroupPrefixes = len(groupPrefixes)
groupPrefix2row = dict([(groupPrefix, index) for index, groupPrefix in enumerate(groupPrefixes)])
#prefix2row = dict([(prefix, index) for index, prefix in enumerate(prefixes)])
prefix2row = dict([(prefix, groupPrefix2row[getGroupPrefix(prefix)]) for prefix in prefixes])

#frames = {}
prefixAndFrame2Col = {}
groupPrefixColsTaken = {}
for prefix in prefixes:
    groupPrefix = getGroupPrefix(prefix)
    frames = prefix2frames[prefix]
    frames.sort()
    #if prefix in prefixAndFrame2Col.keys():
    if groupPrefix in groupPrefixColsTaken.keys():
        #Note that frame numbering starts at 1 so off-by-1 won't be a problem here
        offset = max(groupPrefixColsTaken[groupPrefix])
    else:
        offset = 0
        
    prefixAndFrame2Col[prefix] = dict([(frame, index+offset) for index, frame in enumerate(frames)])
    groupPrefixColsTaken[groupPrefix] = groupPrefixColsTaken.get(groupPrefix,[])+[i+offset+1 for i, frame in enumerate(frames)]

#memberships = np.zeros((nClusters, nPrefixes, max([len(i) for i in prefix2frames.values()])))
#memberships = np.zeros((nClusters, nGroupPrefixes, max([len(i) for i in prefix2frames.values()])))
memberships = np.zeros((nClusters, nGroupPrefixes, max([max(i) for i in groupPrefixColsTaken.values()])))
print [(groupPrefix, max(i)) for groupPrefix, i in zip(groupPrefixColsTaken.keys(), groupPrefixColsTaken.values())]
#1/0
memberships -= 1
for clusterIndex, cluster in enumerate(data):
    for member in cluster:
        result = re.findall('([a-zA-Z0-9_-]+)_([0-9]+)', member)
        prefix = result[0][0].strip('_')
        frame = int(result[0][1])
        
        row = prefix2row[prefix]
        col = prefixAndFrame2Col[prefix][frame]
        #print member
        if sum(memberships[:,row,col]) > 0:
            print memberships[:,row,col]
            1/0
        memberships[:,row,col] = 0
        memberships[clusterIndex,row,col] = 1
        

    pylab.subplot(nRows, nCols, clusterIndex+1)
    pylab.imshow(memberships[clusterIndex,:,:],
                 interpolation='nearest',
                 vmin=-1, vmax=1,
                 aspect='auto')
    #pylab.yticks(range(len(prefixes)), [[]*4+[prefix] for i, prefix in enumerate(prefixes) if i%5==0])
    #pylab.yticks(range(0,len(prefixes),1), prefixes )
    pylab.yticks(range(0,len(groupPrefixes),1), groupPrefixes )
    pylab.xlabel('Frame')
    pylab.title('Cluster %i' %(clusterIndex))
    


pylab.subplots_adjust(left=0.06,
                      right=0.99,
                      bottom=0.025,
                      top=0.975,
                      hspace=0.250)

pylab.show()
                 
        
#Analyze inter-cluster transitions

#nSimulations = len(prefixes)
nSimulations = memberships.shape[0]
nCols = 3
# First +1 because the first plot is of all the simulations summed up
nRows = (nSimulations+1+(nCols-1))/nCols

allTransitions = np.zeros((nClusters, nClusters))
#transitions = dict(((prefix,np.zeros((nClusters, nClusters))) for prefix in prefixes))
transitions = dict(((groupPrefix,np.zeros((nClusters, nClusters))) for groupPrefix in groupPrefixes))
for simulationIndex, groupPrefix in enumerate(groupPrefixes):
    simulationIndex = groupPrefixes.index(groupPrefix)
    for prefix in groupPrefix2prefixes2frames[groupPrefix].keys():
        currentCluster = None
        sliceStartCol = min(prefixAndFrame2Col[prefix].values())
        sliceEndCol = max(prefixAndFrame2Col[prefix].values())
        #simSlice = memberships[:,simulationIndex,:]
        #for frameIndex in range(simSlice.shape[1]):
        for frameIndex in range(sliceStartCol, sliceEndCol+1):
            nonzeros = np.nonzero(memberships[:,simulationIndex,frameIndex]==1)
            if len(nonzeros[0]) == 0:
                nextCluster = None
            elif len(nonzeros[0]) > 1:
                print nextCluster
                1/0
            else:
                nextCluster = nonzeros[0][0]
            #if not(currentCluster is None) and (nextCluster is None):
            #    print prefix, 'CurrentCluster is %r and nextCluster is %r' %(currentCluster, nextCluster)
            if not(currentCluster is None) and not(nextCluster is None):
                allTransitions[currentCluster, nextCluster] += 1
                transitions[groupPrefix][currentCluster, nextCluster] += 1
                #if currentCluster == 7:
                #    print prefix, 'currentCluster is 7 and nextCluster is %r' %(nextCluster)
            currentCluster = nextCluster
            #print nextCluster
    #pylab.subplot(nRows, nCols, simulationIndex+1)
    #pylab.imshow(np.log10(transitions[prefix]),
    #             interpolation='nearest')
    #pylab.colorbar()
#pylab.show()


# Network diagram of cluster transitions

import networkx as nx
G=nx.Graph()
#nonzeros = np.nonzero(allTransitions)
#print nonzeros
#for row, col in zip(nonzeros[0], nonzeros[1]): 
nCols = 5
nRows = ((nGroupPrefixes) / nCols) 
#pylab.subplot(nRows, nCols, 1)
pylab.subplot(1, 1, 1)

for row in range(nClusters):
    for col in range(nClusters):
        G.add_edge(row,col,weight=allTransitions[row,col]+allTransitions[col,row])

 
elarge=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] >5]
emed=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] <=5 and d['weight'] > 2.5]
esmall=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] <=2.5 and d['weight'] > 0]


#pos = nx.spring_layout(G, k=2./np.sqrt(nClusters))
pos = nx.spring_layout(G, k=3./np.sqrt(nClusters))

#nodes

nFramesPerClusterAllSims = [np.sum(memberships[clusterIndex,:,:]==1) for clusterIndex in range(nClusters)]
nodeSizes = [150. * float(i)/(np.mean(nFramesPerClusterAllSims)) for i in nFramesPerClusterAllSims]
nx.draw_networkx_nodes(G, pos, node_size=nodeSizes)

# edges
nx.draw_networkx_edges(G,pos,edgelist=elarge,
                    width=6)
nx.draw_networkx_edges(G,pos,edgelist=emed,
                       width=6,alpha=1.,edge_color='b',style='dashed')
nx.draw_networkx_edges(G,pos,edgelist=esmall,
                       width=6,alpha=0.3,edge_color='b',style='dashed')

# labels

nx.draw_networkx_labels(G,pos,font_size=20,font_family='sans-serif')
pylab.title('All')
pylab.xticks([],[])
pylab.yticks([],[])
pylab.show()

for simulationIndex, groupPrefix in enumerate(groupPrefixes):
    pylab.subplot(nRows, nCols, simulationIndex+1)
    G=nx.Graph()
    #nonzeros = np.nonzero(transitions[prefix])
    #for row, col in zip(nonzeros[0], nonzeros[1]):
    for row in range(nClusters):
        for col in range(nClusters):
            G.add_edge(row,col,weight=transitions[groupPrefix][row,col]+transitions[groupPrefix][col,row])


    elarge=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] >5]
    emed=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] <=5 and d['weight'] > 2.5]
    esmall=[(u,v) for (u,v,d) in G.edges(data=True) if d['weight'] <=2.5 and d['weight'] > 0]

    #pos = nx.spring_layout(G)

    # nodes
    nFramesPerClusterThisSim = [np.sum(memberships[clusterIndex,simulationIndex,:]==1) for clusterIndex in range(nClusters)]
    print nFramesPerClusterThisSim
    nodeSizes = [150. * float(i)/(np.mean(nFramesPerClusterAllSims)/nGroupPrefixes) for i in nFramesPerClusterThisSim]
    print nodeSizes
    nx.draw_networkx_nodes(G, pos, node_size=nodeSizes)

    
    # edges
    nx.draw_networkx_edges(G,pos,edgelist=elarge,
                           width=6)
    nx.draw_networkx_edges(G,pos,edgelist=emed,
                           width=6,alpha=1.,edge_color='b',style='dashed')
    nx.draw_networkx_edges(G,pos,edgelist=esmall,
                           width=6,alpha=0.3,edge_color='b',style='dashed')

    # labels
    nx.draw_networkx_labels(G,pos,font_size=20,font_family='sans-serif')
    
    pylab.xlim([-0.2, 1.2])
    pylab.ylim([-0.2, 1.2])
    pylab.title(groupPrefix)
    pylab.xticks([],[])
    pylab.yticks([],[])
pylab.show()

#print transitions
## Get simulation similarities by cluster vector overlap
groupClusterVecSimilarities = np.zeros((nGroupPrefixes,nGroupPrefixes))
for index1, groupPrefix1 in enumerate(groupPrefixes):
    
    clusterPops1 = np.diag(transitions[groupPrefix1])
    clusterVect1 = clusterPops1 / np.linalg.norm(clusterPops1)
    for index2, groupPrefix2 in enumerate(groupPrefixes):
        clusterPops2 = np.diag(transitions[groupPrefix2])
        clusterVect2 = clusterPops2 / np.linalg.norm(clusterPops2)
        groupClusterVecSimilarities[index1, index2] = np.dot(clusterVect1,clusterVect2)


## Get simulation similarities by average tanimoto score
index2frameData = open('indexMapToFrames.csv').readlines()

prefix2indices = {}
index2prefix = {}
for line in index2frameData:
    index = line.split(',')[0]
    prefix = line.split('/')[-1].split('_')[0]
    prefix2indices[prefix] = prefix2indices.get(prefix,[]) + [int(index)]
    index2prefix[index] = prefix

sortedPrefixes = prefix2indices.keys()
sortedPrefixes.sort()
groupOverlapSimilarities = np.zeros((nGroupPrefixes, nGroupPrefixes))
tanimoto = np.load('tanimoto_matrix.npy')
for index1, prefix1 in enumerate(sortedPrefixes):
    for index2, prefix2 in enumerate(sortedPrefixes):
        if index2 < index1:
            continue
        indices1 = np.array(prefix2indices[prefix1])
        #print indices1
        indices2 = np.array(prefix2indices[prefix2])
        avSimSimilarity = np.mean(tanimoto[indices1, indices2])
        groupOverlapSimilarities[index1, index2] = avSimSimilarity
        groupOverlapSimilarities[index2, index1] = avSimSimilarity


pylab.subplot(1,3,2)

pylab.imshow(groupClusterVecSimilarities,
             interpolation='nearest',
             vmin=0., vmax=1.,
             origin='lower',
             extent=(-0.5,
                      0.5+groupClusterVecSimilarities.shape[0],
                      -0.5,
                      0.5+groupClusterVecSimilarities.shape[1])
             )

pylab.title('Normalized Cluster Vector Overlap')
maxTicLoc = len(groupPrefixes)+1.
ticLocs = np.arange(0, maxTicLoc,maxTicLoc/(len(groupPrefixes))) 
pylab.xticks(ticLocs, groupPrefixes, rotation=90)
pylab.yticks(ticLocs, groupPrefixes)
#pylab.xticks(np.arange(len(groupPrefixes))+0.5, groupPrefixes)
#pylab.yticks(np.arange(len(groupPrefixes))+0.5, groupPrefixes)
#pylab.colorbar()
#pylab.show()

pylab.subplot(1,3,1)
pylab.imshow(groupOverlapSimilarities,
             interpolation='nearest',
             vmin=0., vmax=1.,
             origin='lower',
             extent=(-0.5,
                      0.5+groupOverlapSimilarities.shape[0],
                      -0.5,
                      0.5+groupOverlapSimilarities.shape[1])
             )

pylab.title('Average Simulation Tanimoto')
maxTicLoc = len(groupPrefixes)+1.
ticLocs = np.arange(0, maxTicLoc,maxTicLoc/(len(groupPrefixes))) 
pylab.xticks(ticLocs, groupPrefixes, rotation=90)
pylab.yticks(ticLocs, groupPrefixes)
#pylab.xticks(range(len(groupPrefixes)), groupPrefixes)
#pylab.yticks(range(len(groupPrefixes)), groupPrefixes)
#pylab.colorbar()
#pylab.show()

pdb2Lig = {#'1BYQ':'ADP',
           '1UYF':'PU1',
           '1UYI':'PUZ',
           '1UYL':'NOLIG',
           #'2VCI':'2GJ',
           #'2WI7':'2KL',
           #'2XHR':'C0P',
           #'2YEJ':'ZZ3',
           #'2YEJ':'XQK',
           #'3B26':'B2L',
           '3D0B':'SNX',
           #'3HEK':'BD0',
           #'3K98':'1RC',
           #'3K99':'PFT',
           #'3RKZ':'06T', 
           #'3RLR':'3RR',
           #'4CWN':'6LV',
           #'4FCR':'0TM',
           #'4LWE':'FJ2',
           '4R3M':'JR9',
           #'4W7T':'3JC'
           } 

pdb2LigMol = {}
import urllib2
sortedPdbList = pdb2Lig.keys()
sortedPdbList.sort()
for pdb in sortedPdbList:
    print pdb, pdb2Lig[pdb]
    if pdb2Lig[pdb] == 'NOLIG':
        pdb2LigMol[pdb] = Chem.Mol()
        continue
    url = 'http://www.rcsb.org/pdb/download/downloadLigandFiles.do?ligandIdList=%s&structIdList=%s&instanceType=all&excludeUnobserved=false&includeHydrogens=false' %(pdb2Lig[pdb],pdb)
    #print url
    response = urllib2.urlopen(url)
    #print "Response received"
    import time
    time.sleep(0.5)
    data = response.read()
    #with open('temp_orig.sdf','wb') as of:
    #    of.write(data)
    data = fixPdbSdfFormatting(data)
    #with open('temp.sdf','wb') as of:
    #    of.write(data)
    #os.system('babel -isdf temp.sdf -x3 -osdf  temp2.sdf >& xx')
    #stio = StringIO.StringIO(data)
    #print data
        
    #suppl = Chem.SDMolSupplier('temp.sdf')

    suppl = Chem.SDMolSupplier()
    suppl.SetData(data)
    for mol in suppl:
        #smiles = Chem.MolToSmiles(mol)
        #print smiles
        m_mol = Chem.MolFromSmiles(mol.GetProp('SMILES'))
        AllChem.Compute2DCoords(m_mol)
        pdb2LigMol[pdb] = Chem.Mol(m_mol)
        break #Only read first mol


legends = []
molList = []

sortedPdbList2 = pdb2LigMol.keys()
sortedPdbList2.sort()
for pdb in sortedPdbList2:
    legends.append(pdb)
    molList.append(pdb2LigMol[pdb])


    
#>>> m1 = Chem.MolFromSmiles('Cc1ccccc1')
#>>> fp1 = AllChem.GetMorganFingerprint(m1,2)
#>>> fp1
#<rdkit.DataStructs.cDataStructs.UIntSparseIntVect object at 0x...>
#>>> m2 = Chem.MolFromSmiles('Cc1ncccc1')
#>>> fp2 = AllChem.GetMorganFingerprint(m2,2)
#>>> DataStructs.DiceSimilarity(fp1,fp2)

from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
molSimilarities = np.zeros((len(sortedPdbList2), len(sortedPdbList2)))
morganFingerprints = False
for idx1, pdb1 in enumerate(sortedPdbList2):
    if morganFingerprints:
        molFp1 = AllChem.GetMorganFingerprint(molList[idx1],3, useFeatures=True)
    else:
        molFp1 = FingerprintMols.FingerprintMol(molList[idx1])
        
    for idx2, pdb2 in enumerate(sortedPdbList2):
        if morganFingerprints:
            molFp2 = AllChem.GetMorganFingerprint(molList[idx2],3, useFeatures=True)
            sim = DataStructs.DiceSimilarity(molFp1, molFp2)
        else:
            molFp2 = FingerprintMols.FingerprintMol(molList[idx2])
            sim = DataStructs.FingerprintSimilarity(molFp1, molFp2)
        molSimilarities[idx1,idx2] = sim

pylab.subplot(1,3,3)
pylab.imshow(molSimilarities,
             interpolation='nearest',
             vmin=0., vmax=1.,
             origin='lower',
             extent=(-0.5,
                      0.5+molSimilarities.shape[0],
                      -0.5,
                      0.5+molSimilarities.shape[1])
             )

pylab.title('Ligand Similarity')
maxTicLoc = len(groupPrefixes)+1.
ticLocs = np.arange(0, maxTicLoc,maxTicLoc/(len(groupPrefixes))) 
pylab.xticks(ticLocs, groupPrefixes, rotation=90)
pylab.yticks(ticLocs, groupPrefixes)
#pylab.xticks(range(len(legends)), legends)
#pylab.yticks(range(len(legends)), legends)
#pylab.colorbar()
pylab.show()

flatGroupOverlapSimilarities = groupOverlapSimilarities[np.triu_indices(nGroupPrefixes,1)]
flatMolSimilarities = molSimilarities[np.triu_indices(nGroupPrefixes,1)]

import random
import scipy.stats
print "Calculating Rank Correlations"
actualTau = scipy.stats.kendalltau(flatGroupOverlapSimilarities,
                                  flatMolSimilarities)
print "Actual Kendall Tau:", actualTau
actualRho = scipy.stats.spearmanr(flatGroupOverlapSimilarities,
                                  flatMolSimilarities)
print "Actual Spearman's Rho:", actualRho
nResamplings = 10000
bootstrapRhos = []
bootstrapTaus = []
#flatGroupOverlapSimilarities = np.random.random(flatGroupOverlapSimilarities.shape)
#flatMolSimilarities = np.random.random(flatMolSimilarities.shape)
for i in range(nResamplings):
    newOrder = range(len(flatMolSimilarities))
    random.shuffle(newOrder)
    scrambledMolSimilarities = [flatMolSimilarities[i] for i in newOrder]
    scrambledRho = scipy.stats.spearmanr(flatGroupOverlapSimilarities,
                                         scrambledMolSimilarities)
    bootstrapRhos.append(scrambledRho[0])
    scrambledTau = scipy.stats.kendalltau(flatGroupOverlapSimilarities,
                                          scrambledMolSimilarities)
    bootstrapTaus.append(scrambledTau[0])
    
nBins = 50
pylab.subplot(1,2,1)
hist = pylab.hist(bootstrapRhos,
                  bins=nBins,
                  cumulative=True,
                  histtype='step')
pylab.axvline(actualRho[0], 0, 1,
              color='r',
              linewidth=4)
pylab.title("Spearman's Rho Bootstrap")
pylab.xlabel("Spearman's Rho")
pylab.ylabel('Population')

pylab.subplot(1,2,2)
hist = pylab.hist(bootstrapTaus,
                  bins=nBins,
                  cumulative=True,
                  histtype='step')
pylab.axvline(actualTau[0], 0, 1,
              color='r',
              linewidth=4)
pylab.title("Kendall Tau Bootstrap")
pylab.xlabel("Kendall Tau")
pylab.ylabel('Population')
pylab.show()                     

'''
#Show cluster representatives
showFirstNClusters = 6
nAngles = 1

nCols = 4
nRows = (showFirstNClusters / 2) + 1
import matplotlib.image as mpimg
import os
#for clusterNum in range(nClusters):
for clusterNum in range(showFirstNClusters):
    thisRow = (clusterNum / 2) + 1
    thisCol = ((nAngles+1) * (clusterNum % 2))+1
    pylab.subplot(nRows,nCols, ((thisRow-1)*nCols)+thisCol)
    pylab.text(0,0,'Cluster %i'%(clusterNum), fontsize=30)
    #for angle in range(1,5):
    for angle in range(1,nAngles+1):
        if not(os.path.exists('clusterPics/cluster%i_render%i.png')):
            os.system('convert clusterPics/cluster%i_render%i.tga clusterPics/cluster%i_render%i.png' %(clusterNum, angle, clusterNum, angle))
        
        img = mpimg.imread('clusterPics/cluster%i_render%i.png' %(clusterNum, angle))
        thisCol = ((nAngles+1) * (clusterNum % 2)) + angle + 1
        pylab.subplot(nRows, nCols,((thisRow-1)*nCols)+thisCol)
        pylab.imshow(img)
        pylab.xticks([],[])
        pylab.yticks([],[])
pylab.show()
'''
