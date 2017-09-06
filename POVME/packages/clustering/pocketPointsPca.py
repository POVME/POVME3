#!python





import numpy as np
import cPickle
import glob
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA as sklearnPCA
import pylab
#import hsp90Helpers
import POVME.packages.binana.peel as peel
import os
import re
pdb2Lig = {#'1BYQ':'ADP',
           '1UYF_every50ns':'PU1',
           '1UYI_every50ns':'PUZ',
           '1UYL_every50ns':'NOLIG',
           #'2VCI':'2GJ',
           #'2WI7':'2KL',
           #'2XHR':'C0P',
           #'2YEJ':'ZZ3',
           #'2YEJ':'XQK',
           #'3B26':'B2L',
           '3D0B_every50ns':'SNX',
           #'3HEK':'BD0',
           #'3K98':'1RC',
           #'3K99':'PFT',
           #'3RKZ':'06T', 
           #'3RLR':'3RR',
           #'4CWN':'6LV',
           #'4FCR':'0TM',
           #'4LWE':'FJ2',
           '4R3M_every50ns':'JR9',
           #'4W7T':'3JC'
           } 
sortedPdbs = pdb2Lig.keys()
sortedPdbs.sort()



every = 1

analysisPath = '../../2-POVME_analysis/'
#groupPrefixes = glob.glob('%s/????' %(analysisPath))
#groupPrefixes = [i.replace(analysisPath+'/','') for i in groupPrefixes]
#print groupPrefixes

# data = open('indexMapToFrames.csv').readlines()
# datasp = [i.split(',') for i in data]
# indices = [int(i[0]) for i in datasp]
# prefixes = [i[1].split('/')[-1].split('_')[0] for i in datasp]

index2frameData = open('indexMapToFrames.csv').readlines()

prefix2indices = {}
prefix2frames = {}
index2prefix = {}
prefixList = []
frameList = []
# Parse indexMapToFrames line, for example
# 0,../../2-POVME_analysis/1UYF_every50ns/1UYF_every50ns_frameInfo/1UYF_every50ns_frame_1.npy
#|index|                                                          |    prefix    |   |frame|     
for line in index2frameData:
    index = int(line.split(',')[0])
    prefix = line.split('/')[-1].split('_frame')[0]
    print prefix
    frame = int(line.split('/')[-1].split('_')[-1].replace('.npy',''))
    frameList.append(frame)
    prefix2indices[prefix] = prefix2indices.get(prefix,[]) + [int(index)]
    prefix2frames[prefix] = prefix2frames.get(prefix,[]) + [int(frame)]
    index2prefix[index] = prefix

#sortedPrefixes = prefix2indices.keys()
#sortedPrefixes.sort()


## Load all relevant frames
print "Loading frames"
pointsSetList = []
frameStringList = []
allPointsSet = set()
prefixList = []
for prefix in sortedPdbs:
    maxFrame = max(prefix2frames[prefix])
    framesToLoad = range(1,maxFrame,every)
    for frame in framesToLoad:
        frameString = '%s_frame_%i' %(prefix, frame)
        frameFilename = '%s/%s/%s_frameInfo/%s.npy' %(analysisPath,
                                                      prefix,
                                                      prefix,
                                                      frameString)
        prefixList.append(prefix)
        frameStringList.append(frameString)
        pointsSet = set([tuple(i) for i in np.load(frameFilename)])
        pointsSetList.append(pointsSet)
        allPointsSet = allPointsSet.union(pointsSet)

        
print "Vectorizing frames"
coord2vectPos = dict([(coord, i) for i, coord in enumerate(allPointsSet)])
vectPos2Coord = dict([(i, coord) for i, coord in enumerate(allPointsSet)])
nVectPos = len(allPointsSet)
nFrames = len(frameStringList)
vecPosMatrix = np.zeros((nFrames, nVectPos), dtype=np.bool)
for f1 in range(nFrames):
    for point in pointsSetList[f1]:
        vecPosMatrix[f1,coord2vectPos[point]] = 1
#print vecPosMatrix
print 'vecPosMatrix.shape is ', vecPosMatrix.shape

print "Making X_std"

normalizeFeatures = 'Mean'
if normalizeFeatures == 'MeanAndStd':
    X_std = StandardScaler().fit_transform(vecPosMatrix)
elif normalizeFeatures == 'Mean':
    X_std = vecPosMatrix - np.mean(vecPosMatrix,0)

print "Making sklearn_pca"
nComponents = 10
sklearn_pca = sklearnPCA(n_components=nComponents)
print "Making Y_sklearn"

# An important note -- we get a hash value for the raw data here
# as a unique identifier of this dataset. Then we save the initial
# PCA object to disk to reduce time on subsequent runs.

sklearn_pca_hash = str(hash(tuple(X_std.flatten()[::100])))[-10:]
#pca_result_hash = 'y_sklearn_hash_%s.cpickle' %(sklearn_pca_hash)
pca_result_hash = 'pca_obj_hash_%s.cpickle' %(sklearn_pca_hash)
if not(os.path.exists(pca_result_hash)):
    #Y_sklearn = sklearn_pca.fit_transform(X_std)
    sklearn_pca.fit(X_std)
    with open(pca_result_hash,'wb') as of:
        cPickle.dump(sklearn_pca, of)
    Y_sklearn = sklearn_pca.transform(X_std)
    #with open(pca_result_hash,'wb') as of:
    #    cPickle.dump(Y_sklearn, of)
else:
    #Y_sklearn = cPickle.load(open(pca_result_hash))
    sklearn_pca = cPickle.load(open(pca_result_hash))
    Y_sklearn = sklearn_pca.transform(X_std)

print 'Y_sklearn.shape', Y_sklearn.shape
for col in range(Y_sklearn.shape[1]):
    print col, np.mean(abs(Y_sklearn[:,col]))

twoD = True
threeD = False
plotExplainedVariance = True

savePcDxFiles = True
nCircImgPlaces = 20
#showMolAtAvgPos = True
#contourNotScatter = False


showMolAtAvgPos = False
contourNotScatter = False

scatters = []
if contourNotScatter:
    colors = ['b', 'c', 'g', 'm', 'r', 'k', '0.5'] 
    markers = ['solid','dashed','dashdot']

if not(contourNotScatter):
    colors = ['b', 'c', 'g', 'm', 'r', 'k'] 
    markers = ['x', 'o', 's', 'v']


colorsMarkers = [[(j,i) for i in markers] for j in colors] # This will be a list of lists
colorsMarkers = [i for sl in colorsMarkers for i in sl] # so we flatten it


# The ligand figure circular arrangement will have a radius of the PC value range times this factor
radiusRangeFactor = 0.65
from matplotlib.offsetbox import AnnotationBbox, OffsetImage
import scipy.spatial.distance as ssd
import itertools
if twoD:
    combos = list(itertools.combinations(range(4),2))
    for index, (pc_x, pc_y) in enumerate(combos):
        xMin = min(Y_sklearn[:,pc_x])
        xMax = max(Y_sklearn[:,pc_x])
        yMin = min(Y_sklearn[:,pc_y])
        yMax = max(Y_sklearn[:,pc_y])
        #print 'xMin, xMax, yMin, yMax', xMin, xMax, yMin, yMax
        # These ranges define how the circular positions of ligand figures are placed. The algorithm needs a circle for assignment, but will convert it to an ellipse for drawing
        
        yRange = yMax - yMin
        xRange = yRange
        #recoverXRangeFactor = (float(xMax - xMin) / yRange) ** 0.5
        recoverXRangeFactor = 1.
        
        #yRange = max(Y_sklearn[:,pc_y]) - min(Y_sklearn[:,pc_y])
        #xOffset = -np.mean(Y_sklearn[:,pc_x])
        #yOffset = -np.mean(Y_sklearn[:,pc_y])
        xOffset = xMin + (xRange/2)
        yOffset = yMin + (yRange/2)
        xyCircPlaces = [[(xRange*np.sin(i))+xOffset,
                         (yRange*np.cos(i))+yOffset] 
                         for i in np.arange(0,2*np.pi, 2*np.pi/nCircImgPlaces)]
        xyCircPlaces = np.array(xyCircPlaces)
        xyCircPlaces *= radiusRangeFactor
        circPlacesTaken = np.zeros((nCircImgPlaces))
        
        # Prioritize which molecule reps to draw first
        toPlot = []
        for pdb, (color, marker) in zip(sortedPdbs, colorsMarkers):
            indices = [i for i,thisPrefix in enumerate(prefixList) if thisPrefix==pdb]
            indices = np.array(indices)
            xs = Y_sklearn[indices,pc_x]
            ys = Y_sklearn[indices,pc_y]
            #angle = np.arctan2(np.mean(xs),np.mean(ys))
            #priority = angle
            xAvg = np.mean(xs) - xOffset
            yAvg = np.mean(ys) - yOffset
            priority = ((xAvg**2) + (yAvg**2)) ** 0.5
            priority = 1.0/priority
            #priority = np.random.random()
            #distanceToSpots = ssd.cdist([[xAvg,yAvg]],
            #                            xyCircPlaces)
            #priority = np.min(distanceToSpots)
            toPlot.append((pdb,color,marker,priority))
        #print 'Priority values:', [i[3] for i in toPlot]
        toPlot.sort(key=lambda x:x[3])
        
        
        scatters = []
        #print colorsMarkers
        #ax = pylab.subplot(1,1,index+1)
        ax = pylab.subplot(1,1,1)
        #for pdb, (color, marker) in zip(sortedPdbs, colorsMarkers):
        for pdb, color, marker, angle in toPlot:
            indices = [i for i,thisPrefix in enumerate(prefixList) if thisPrefix==pdb]
            indices = np.array(indices)
            xs = Y_sklearn[indices,pc_x]
            ys = Y_sklearn[indices,pc_y]

            if not(contourNotScatter):
                scatter = pylab.scatter(xs, ys,
                                        color=color,
                                        marker=marker)
            if contourNotScatter:
                scatter = pylab.plot([], [],
                                     color=color,
                                     linestyle=marker
                                     )

                ### Contour code ###
                nBins = 100
                contourLevel = 0.03
                xBins = np.arange(xMin, xMax+0.0001, (xMax-xMin)/nBins)
                yBins = np.arange(yMin, yMax+0.0001, (yMax-yMin)/nBins)
                hist2d, xEdges, yEdges = np.histogram2d(xs, ys, bins=[xBins, yBins])
                import scipy.ndimage.filters
                hist2d = scipy.ndimage.filters.gaussian_filter(hist2d,3.5)
                print xBins[:-1]
                print
                print yBins[:-1]
                print
                print hist2d
                #xBins += [xBins[-1]+(xRange/25)]
                #yBins += [yBins[-1]+(yRange/25)]
                contour = pylab.contour(xBins[:-1],
                                        yBins[:-1],
                                        hist2d.T,
                                        [contourLevel],
                                        #[contourLevel,100],
                                        label=pdb,
                                        colors=color,
                                        linestyles=marker,
                                        #hatches = ['/','\\','|','-','+','x','o','O','.','*']*2
                                        )
                pylab.clabel(contour,
                             fmt={contourLevel:pdb},
                             #fontsize='large'
                             )
                ### End countour code ###
            pylab.xlabel('PC %s' %(pc_x+1))
            pylab.ylabel('PC %s' %(pc_y+1))
            if showMolAtAvgPos:
                #print ax
                #size = (150,150)
                #img = hsp90Helpers.getLigandAsImg(pdb, size=size)
                img = hsp90Helpers.getLigandAsImg(pdb)
                #pylab.figimage(np.array(img),np.mean(xs), np.mean(ys))
                imageBox = OffsetImage(np.array(img), zoom=1.0)
                xAvg = np.mean(xs)
                yAvg = np.mean(ys)

                distanceToSpots = ssd.cdist([[xAvg,yAvg]],
                                            xyCircPlaces)
                closestSpots = np.argsort(distanceToSpots)[0]
                spotFound = False
                for spotIndex in closestSpots:
                    if circPlacesTaken[spotIndex]==0:
                        thisSpot = xyCircPlaces[spotIndex,:]
                        circPlacesTaken[spotIndex] = 1
                        spotFound = True
                        break
                if spotFound == False:
                    1/0
                
                ab = AnnotationBbox(imageBox, (xAvg, yAvg),
                                    #xybox=size,
                                    #xybox=((1.5+np.random.random((1)))*xAvg,
                                    #       (1.5+np.random.random((1)))*yAvg),
                                    xybox=(((thisSpot[0]-xOffset)*recoverXRangeFactor)+xOffset,
                                           thisSpot[1]),
                                    #xybox=(thisSpot),
                                    #xycoords='data',
                                    #xycoords=(1.5*xAvg, 1.5*yAvg),
                                    #boxcoords="offset points",
                                    #boxcoords=(1.5*xAvg, 1.5*yAvg),
                                    #pad=0.5,
                                    pad=0.1,
                                    #arrowprops=dict(arrowstyle="->",
                                    arrowprops=dict(arrowstyle="simple",
                                                    facecolor='k'
                                                    #connectionstyle="angle,angleA=0,angleB=90,rad=3"
                                                    )
                                    )
                ax.add_artist(ab)

            scatters.append((scatter,pdb))
        scatters.sort(key=lambda x:x[1])
        if not(contourNotScatter):
            pylab.legend([i[0] for i in scatters], sortedPdbs)
        pylab.show()

if threeD:
    
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    for pdb, (color, marker) in zip(sortedPdbs, colorsMarkers):
        indices = [i for i,thisPrefix in enumerate(prefixList) if thisPrefix==pdb]
        indices = np.array(indices)
        xs = Y_sklearn[indices,0]
        ys = Y_sklearn[indices,1]
        zs = Y_sklearn[indices,2]
        scatter = pylab.scatter(xs, ys, zs=zs,
                                color=color,
                                marker=marker)
        scatters.append(scatter)
    pylab.legend(scatters, sortedPdbs)
    pylab.show()

if plotExplainedVariance:
    EVR = sklearn_pca.explained_variance_ratio_
    pylab.bar(np.arange(0.5, 0.5+len(EVR)), EVR)
    pylab.plot([sum(EVR[:i]) for i in range(len(EVR))])
    pylab.ylabel('Explained variance ratio (individual/cumulative)')
    pylab.xlabel('Principal Component')
    pylab.show()





templateVisualizePc = '''mol new {%s} waitfor all
#mol addrep !MOLID!
mol modmaterial 0 top HardPlastic
mol modstyle 0 top Isosurface 0.02500000 0 0 0 1 1
#mol modstyle 0 !MOLID! Isosurface 0.02500000 0 2 1 1 1
#green
mol modcolor 0 top ColorID 12

mol addfile {%s} waitfor all
mol addrep top
#mol addrep !MOLID!
mol modmaterial 1 top HardPlastic
mol modstyle 1 top Isosurface -0.02500000 1 0 0 1 1
#mol modstyle 1 !MOLID! Isosurface -0.02500000 1 2 1 1 1
#red
mol modcolor 1 top ColorID 1


'''

    
if savePcDxFiles:
    loadAllScript = 'display rendermode GLSL\ndisplay projection Orthographic\n'
    print 'Saving Principal Component DX files'
    for pc_ind in range(nComponents):
        print "Saving DX file for PC %i" %(pc_ind+1)
        pcPoints = np.zeros((vecPosMatrix.shape[1],4))
        for vectPos in vectPos2Coord:
            coord = vectPos2Coord[vectPos]
            pcPoints[vectPos,:3] = coord
            pcPoints[vectPos,3] = sklearn_pca.components_[pc_ind,vectPos]
        m_fm = peel.featureMap.fromOffGridPts(pcPoints, 1., justCoords=False)
        pcFilename = 'PC%i.dx' %(pc_ind+1)
        m_fm.write_dx_file(pcFilename)
        loadAllScript += templateVisualizePc.replace('!MOLID!',str(pc_ind)) %(pcFilename, pcFilename)
    with open('loadAllPcs.vmd', 'wb') as of:
        of.write(loadAllScript)
    
print "Done!"

