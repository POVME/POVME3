import pylab
import numpy as np
import re
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
        
for clusterIndex, cluster in enumerate(data):
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
