#!python
# Implementation of Clustering Algorithms in POVME
# By Celia Wong
# Advised by Jeff Wagner
# Amaro Lab, UCSD

import scipy.cluster.vq, scipy.cluster.hierarchy
import argparse
import numpy
import sys
import os
import csv
import copy
import itertools
import collections
#import fnmatch
import pylab
import POVME.packages.binana.peel as peel
#import matplotlib.pyplot

class InputReader():

    def __init__(self):
        #self.coordinates = []
        #self.frames = 0
        self.overlapMatrix = []
        self.prefixToTrajectory = {}
        self.indexToNpyFile = {}
        self.indexToFrame = {}
        self.indexToPrefix = {}

    # Save each frame in the trajectory as a set
    ''' def read_traj(self,traj_file):

        trajectory = open(traj_file,'r')
        frame_coordinates = []
#        temp_coordinates = numoy.array([])

        for line in trajectory:
            if line[0] == 'E':
                self.coordinates.append(frame_coordinates)

#                numpy.append(self.coordinates,frame_coordinates)
                frame_coordinates = []

                self.frames += 1

            elif line[0] == 'A':


            #output_file_name = 'hierarchical_' + command_input['output_name'] + '.csv'

                if line[17] != 'X':
                    frame_coordinates.append((float(line[29:37].strip()),float(line[38:45].strip()),float(line[46:54].strip())))
                    #numpy.append(frame_coordinates,(float(line[29:37].strip()),float(line[38:45].strip()),float(line[46:54].strip())))


        self.coordinates = numpy.array(self.coordinates)

        trajectory.close()
    '''
        
    def read_indexFile(self,indexToFrameFile):
        if indexToFrameFile == None:
            return
        #self.indexToNpyFile = {}
        #self.indexToFrame = {}
        with open(indexToFrameFile) as csvfile:
            fieldnames = ['index','frameFile']
            reader = csv.DictReader(csvfile, fieldnames=fieldnames)
            for row in reader:
                self.indexToNpyFile[int(row['index'])] = row['frameFile']
                if self.indexToFrame != None:
                    try:
                        frameNumber = row['frameFile'].split('frame_')[-1].replace('.npy','')
                        self.indexToFrame[int(row['index'])] = int(frameNumber)
                        framePrefix = row['frameFile'].split('/')[-1].replace('frame_%s.npy'%(frameNumber),'')
                        self.indexToPrefix[int(row['index'])] = framePrefix
                    except:
                        raise Exception("Unable to strip frame number or prefix from input filename %s. Disabling frame number output." %(row['frameFile']))
                        #self.indexToFrame = None
                        #self.indexToPrefix = None
                #self.coordinates.append(numpy.load(row['frame']))
        

    def read_overlap(self,overlap_file):
        overlap_suffix = overlap_file[-4:]
        if overlap_suffix == '.npy':
            self.overlapMatrix = numpy.load(open(overlap_file))
        elif overlap_suffix[-4:] == '.csv':
            overlap = open(overlap_file,'r')
            overlap_values = csv.reader(overlap, delimiter = ',')
    
            for line in overlap_values:
                #self.overlapMatrix.append([float(x) for x in line])
                row = []
                #self.frames += 1
                for value in line:
                    row.append(float(value))
                self.overlapMatrix.append(row)
    
            self.overlapMatrix = numpy.array(self.overlapMatrix)
            #self.oneMinusOverlapMatrix = 1. - numpy.array(self.overlapMatrix)
    
            overlap.close()
        else:
            raise Exception('Unrecognized overlap matrix input file type:', overlap_suffix)

    
    def parse_traj_inputs(self, argst, argsT):
        if argsT != None:
            for argT in argsT:
                data = open(argT).read()
                datasp = data.strip().split()
                for line in datasp:
                    linesp = line.split(':')
                    prefix = linesp[0].strip()
                    trajFile = linesp[1].strip()
                    self.prefixToTrajectory[prefix] = trajFile
        for argt in argst:
            argtsp = argt.split(':')
            prefix = argtsp[0].strip()
            trajFile = argtsp[1].strip()
            self.prefixToTrajectory[prefix] = trajFile
        ## Check to ensure these files exist
        for prefix in self.prefixToTrajectory:
            trajFileName = self.prefixToTrajectory[prefix]
            if not(os.path.exists(trajFileName)):
                raise Exception('ERROR - trajectory file %s doesn\'t exist!' %(trajFileName))

                

class Cluster():
    #def __init__(self,coordinates,frames,overlap_values,frameToFileName):
    def __init__(self,input_reader):
        #self.coordinates = coordinates
        #self.frames = frames
        self.overlap_values = input_reader.overlapMatrix
        self.one_minus_overlap_values = 1. - self.overlap_values
        self.frames = len(self.overlap_values)
        self.indexToFrame = input_reader.indexToFrame
        self.indexToPrefix = input_reader.indexToPrefix
        self.indexToNpyFile = input_reader.indexToNpyFile
        self.prefixToTrajectory = input_reader.prefixToTrajectory
        self.avgLinkage = None
        self.whited_overlap_values = None

        self.do_sanity_checks()

    def do_sanity_checks(self):
        self.check_commandline_inputs()
        self.ensure_file_prefixes_map_to_trajectories()
        #self.ensure_trajectories_exist() #Check performed during -T argument parsing instead
        
    def check_commandline_inputs(self):
        if (self.indexToNpyFile == {}) and (self.prefixToTrajectory != {}):
            raise Exception("ERROR! Given pdb trajectory (-t/T) but not given index file (-i). Output will return matrix indices instead of frame numbers or cluster representative structures.")
        elif self.indexToNpyFile == {}:
            print "Not given index file (-i). Clustering will return matrix indices, but not trajectory frame numbers or members/representatives."
        elif (self.indexToNpyFile != {}) and(self.prefixToTrajectory == {}):
            print "Given index file (-i) but not prefix-to-trajectory mapping (-t or -T). Clustering will return prefix and frame numbers of cluster members, but will not extract representatives."
        elif (self.indexToNpyFile != {}) and (self.prefixToTrajectory != {}):
            print "Given index file (-i) and prefix-to-trajectory mapping (-t or -T). Clustering will return prefix and frame numbers of cluster members, and will extract representatives."

    def ensure_file_prefixes_map_to_trajectories(self):
        if (self.prefixToTrajectory == {}) or (self.indexToNpyFile == {}):
            print "No -i and/or -t/T arguments given - Skipping file-prefix-to-trajectory mapping completeness test"
            return
        else:        
            allPrefixesSet = set(self.indexToPrefix.values())
            for prefix in allPrefixesSet:
                if not(prefix in self.prefixToTrajectory.keys()):
                    raise Exception('File prefix %s not found in -t arguments (which are %r)' %(prefix, self.prefixToTrajectory.keys()))
        return
        
    def ensure_trajectories_exist(self):
        if self.prefixToTrajectory == {}:
            print "No -t/T arguments given. Skipping trajectory-file-existence check"
        else:
            for trajectoryFile in self.prefixToTrajectory.values():
                if not(os.path.exists(trajectoryFile)):
                    raise Exception("Trajectory file %s not found" %(trajectoryFile))
            
            
    def kmeans_cluster(self,number_clusters):

        if self.whited_overlap_values == None:
            self.whited_overlap_values = scipy.cluster.vq.whiten(self.overlap_values)
        frames,result = scipy.cluster.vq.kmeans(self.whited_overlap_values, number_clusters)
        #frames,result = scipy.cluster.vq.kmeans2(self.whited_overlap_values, number_clusters)
        code, dist = scipy.cluster.vq.vq(self.whited_overlap_values,frames)

        print "The clusters are {0}".format(code)
        print code.shape
        list_of_clusters = self.separate_clusters(code)
        #return code
        return list_of_clusters

    def hierarchical_cluster(self, number_clusters):
        if self.avgLinkage == None:
            try:
                overlapHash = str(numpy.sum(self.one_minus_overlap_values.flatten()[::self.one_minus_overlap_values.size/100]))[-7:]
            except:
                overlapHash = str(numpy.sum(self.one_minus_overlap_values.flatten()))[-7:]
            linkageFile = 'avg_linkage_hash_%s.npy' %(overlapHash)
            if os.path.exists(linkageFile):
                self.avgLinkage = numpy.load(linkageFile)
            else:
                self.avgLinkage = scipy.cluster.hierarchy.average(self.one_minus_overlap_values)
                numpy.save(linkageFile, self.avgLinkage)
        result = scipy.cluster.hierarchy.fcluster(self.avgLinkage,
                                                  number_clusters,
                                                  criterion='maxclust')
        
        #result = scipy.cluster.hierarchy.linkage(self.overlap_values)
        # Hierarchical clustering seems to have a bad habit of returning clusters nubered starting from 1 instead of 0. This is undesired behavior.
        result = result - 1
        clusters = self.separate_clusters(result)
        return clusters
    
#        scipy.cluster.hierarchy.dendrogram(result)


    ''' separate_cluster_traj will seperate the original trajectory into the
        number of clusters specified. Each new cluster traj will only contain the
        frames that belong in that cluster.

        cluster_result = list that is number of frames long and contain the cluster
                         each frame is grouped with.
        number_clusters = the number of clusters specified.
        file_name = the file name that was passed into main() as a command line arg.
        traj_file = the original trajectory file containing all frames
    '''
    def separate_cluster_traj(self,cluster_result,number_clusters,file_name,traj_file,output_file):

        list_writes = [None]*number_clusters

        '''Opening n number of clusters (n = number of clusters previously indicated)'''
        for i in range(number_clusters):
            list_writes[i] = open('cluster_'+ str(i)+'.pdb','wb')

        initial_pdb = open(traj_file,'rb')

        current_frame = 0
        current_cluster = cluster_result[current_frame]

        for line in initial_pdb:

            if line[0] == 'E':
                list_writes[current_cluster].write(line)
                if current_frame < len(cluster_result)-1:
                    current_frame += 1
                    current_cluster = cluster_result[current_frame]
            else:
                list_writes[current_cluster].write(line)

        initial_pdb.close()
        for i in range(number_clusters):
            list_writes[i].close()

    ''' Separates the frames into the set clusters
    '''

    def separate_clusters(self,cluster_results):

#        print "cluster results: {0}".format(cluster_results)
        total_num_clusters = len(set(cluster_results))
        list_clusters = [list([]) for i in xrange(total_num_clusters)]

        for i in range(len(cluster_results)):
#            print "Cluster_res for {0} is {1}".format(i, cluster_results[i])
            list_clusters[cluster_results[i]].append(i)
        list_clusters.sort(key=len, reverse=True)
        return list_clusters

    ''' csv file containing differences in binding site is first argument already read and stored into memory by previous command '''

    #def find_centroids(self,binding_volume_matrix,cluster_results,number_clusters,number_frames, indexToFrame):
    def find_centroids(self, list_of_clusters, outputPrefix):
        #number_clusters = len(set(cluster_results))
        #list_of_clusters = self.separate_clusters(cluster_results)
        #print list_of_clusters

        ''' set to some arbitrary large number? '''
        #shortest_average_distance = [1.e20] * number_clusters

        #centroid_list = [[0] for i in xrange(number_clusters)]
        centroid_list = []
        for cluster in list_of_clusters:
            sum_distances = []
            if len(cluster) == 1:
                sum_distances.append(0)
            else:
                cluster = numpy.array(cluster)
                for entry in cluster:
                    allButEntry = cluster[cluster != entry]
                    #print cluster, entry, allButEntry
                    #print self.one_minus_overlap_values[entry,:]
                    totalDist = numpy.sum(self.one_minus_overlap_values[entry,allButEntry])
                    sum_distances.append(totalDist)
                
            #print cluster, sum_distances, numpy.argsort(sum_distances)[0]
            centroid_cluster_index = numpy.argsort(sum_distances)[0]
            centroid_global_index = cluster[centroid_cluster_index]
            centroid_list.append(centroid_global_index)
                
            


        
        if (self.indexToFrame == {}) and (self.indexToNpyFile != {}):
            repsFileName = '%scluster_reps.csv' %(outputPrefix)
            membersFileName = '%scluster_members.csv' %(outputPrefix)
            print "Unable to extract frame numbers from file names. Writing out npy file names to %s and %s" %(repsFileName, membersFileName)
            with open(repsFileName,'wb') as of:
                cluster_rep_file_names = [str(self.indexToNpyFile[i]) for i in centroid_list]
                of.write('\n'.join(cluster_rep_file_names))
            with open(membersFileName,'wb') as of:
                for cluster in list_of_clusters:
                    cluster_member_file_names =[str(self.indexToNpyFile[i]) for i in cluster] 
                    of.write(' '.join(cluster_member_file_names))
                    of.write('\n')
            
        elif (self.indexToFrame == {}) and (self.indexToNpyFile == {}):
            print "No matrix-index-to-trajectory-frame mapping given. Writing out matrix indices"
            with open('%scluster_reps.csv' %(outputPrefix),'wb') as of:
                of.write('\n'.join([str(i) for i in centroid_list]))
            with open('%scluster_members.csv' %(outputPrefix),'wb') as of:
                for cluster in list_of_clusters:
                    of.write(' '.join([str(i) for i in cluster]))
                    of.write('\n')
                    
        elif (self.indexToFrame != {}):
            repsFileName = '%scluster_reps.csv' %(outputPrefix)
            membersFileName = '%scluster_members.csv' %(outputPrefix)
            print "Matrix-index-to-trajectory-frame mapping given. Writing out trajectory frames to %s and %s." %(repsFileName, membersFileName)
            
            with open(repsFileName,'wb') as of:
                cluster_rep_frame_nums = [str(self.indexToFrame[i]) for i in centroid_list]
                cluster_rep_prefixes = [str(self.indexToPrefix[i]) for i in centroid_list]
                cluster_rep_strings = ['_'.join(i) for i in zip(cluster_rep_prefixes, cluster_rep_frame_nums)]
                of.write('\n'.join(cluster_rep_strings))
            with open(membersFileName,'wb') as of:
                for cluster in list_of_clusters:
                    cluster_member_frame_nums =[str(self.indexToFrame[i]) for i in cluster]
                    cluster_member_prefixes = [str(self.indexToPrefix[i]) for i in cluster]
                    cluster_member_strings = ['_'.join(i) for i in zip(cluster_member_prefixes, cluster_member_frame_nums)]
                    of.write(' '.join(cluster_member_strings))
                    of.write('\n')

        if (self.indexToFrame != {}) and (self.prefixToTrajectory != {}):
            print "Extracting trajectory frames"
            matrixIndex2Cluster = {}
            for index, centroid in enumerate(centroid_list):
                matrixIndex2Cluster[centroid] = index
            clusterInd2CentFileName = self.extractFrames(matrixIndex2Cluster, outputPrefix, reps=True)
        else:
            clusterInd2CentFileName = {}
        return clusterInd2CentFileName
        
    def outputAllFrames(self, list_of_clusters, outputPrefix):
        ## check to make sure we'll be able to map all matrix indices to files
        for clusterInd, cluster in enumerate(list_of_clusters):
            #print cluster
            #print indexToFrame.keys()
            for matrixInd in cluster:
                if not(matrixInd in self.indexToFrame.keys()):
                    raise Exception('User requested all frame pdbs to be output to cluster directories, but the program is unable to map all overlap matrix indices to trajectory/frame combinations. Make sure that -t/-T and -i arguments cover all frames and prefixes. Index: %i Cluster: %i' %(matrixInd, clusterInd))

        ## If all mappings exist, extract all relevant frames
        matrixInd2Cluster = {}
        for clusterInd, cluster in enumerate(list_of_clusters):
            for matrixInd in cluster:
                matrixInd2Cluster[matrixInd] = clusterInd
        self.extractFrames(matrixInd2Cluster, outputPrefix, reps=False)
        
            
    
    def extractFrames(self, matrixIndex2Cluster, outputPrefix, reps=False):
        framesToExtract = {}
        clusterInd2CentFileName = {}
        for thisMatrixIndex in matrixIndex2Cluster:
            thisCluster = matrixIndex2Cluster[thisMatrixIndex]
            npyFileName = self.indexToNpyFile[thisMatrixIndex]
            npyFilePrefix = npyFileName.split('/')[-1].split('frame_')[0]
            frameNum = int(npyFileName.split('/')[-1].split('frame_')[-1].replace('.npy',''))
            prefixMatch = ''
            ## See if this prefix is in our dictionary or trajectories
            for trajPrefix in self.prefixToTrajectory.keys():
                if trajPrefix == npyFilePrefix:
                    if prefixMatch == '':
                        prefixMatch = trajPrefix
                    else: # If a matching prefix has already been found
                        raise Exception('ERROR - file %s matches prefix %s and %s' %(npyFileName, trajPrefix, prefixMatch))
            ## Disabled this block - All prefix-to-trajectory matching should be explicit. This caused an error when POVME was run whith a blank prefix
            #if prefixMatch == '':
            #    trajFileName = npyFilePrefix + '.pdb'
            #else:
            trajFileName = self.prefixToTrajectory[prefixMatch]

            ## Figure out the directory and filename that this frame should be written to    
            outputDir = '%scluster%i' %(outputPrefix, thisCluster)
            if not os.path.exists(outputDir):
                os.system('mkdir %s' %(outputDir))
            if reps == True:
                outputFileName = 'REP_%sframe_%i.pdb' %(prefixMatch, frameNum)
                clusterInd2CentFileName[thisCluster] = outputFileName
            else:
                outputFileName = '%sframe_%i.pdb' %(prefixMatch, frameNum)
            fullOutputFileName = outputDir + '/' + outputFileName
            if not trajFileName in framesToExtract.keys():
                framesToExtract[trajFileName] = {}
                
            framesToExtract[trajFileName][frameNum] = fullOutputFileName
            
        for trajFileName in framesToExtract:
            frameCounter = 1
            frameData = ''
            with open(trajFileName) as fo:
                for line in fo:
                    if frameCounter in framesToExtract[trajFileName]:
                        frameData += line
                    if 'END' == line[:3]:
                        if frameData != '':
                            thisOutputFileName = framesToExtract[trajFileName][frameCounter]
                            with open(thisOutputFileName,'wb') as of:
                                of.write(frameData)
                        frameData = ''
                        frameCounter += 1
                if frameData != '':
                    thisOutputFileName = framesToExtract[trajFileName][frameCounter]
                    with open(thisOutputFileName,'wb') as of:
                        of.write(frameData)
        return clusterInd2CentFileName
        
    '''                
    def extractFrame(self, matrixIndex, outputDir, rep=False):
        npyFileName = self.indexToNpyFile[matrixIndex]
        npyFilePrefix = npyFileName.split('/')[-1].split('frame_')[0]
        frameNum = int(npyFileName.split('/')[-1].split('frame_')[-1].replace('.npy',''))
        prefixMatch = ''
        for trajPrefix in self.prefixToTrajectory.keys():
            if trajPrefix == npyFilePrefix:
                if prefixMatch == '':
                    prefixMatch = trajPrefix
                else: # If a matching prefix has already been found
                    raise Exception('ERROR - file %s matches prefix %s and %s' %(npyFileName, trajPrefix, prefixMatch))
                    
        if prefixMatch == '':
            trajFileName = npyFilePrefix + '.pdb'
        else:
            trajFileName = self.prefixToTrajectory[prefixMatch]
            
        if rep == True:
            outputFileName = '%s/REP_%sframe_%i.pdb' %(outputDir, prefixMatch, frameNum)
        else:
            outputFileName = '%s/%sframe_%i.pdb' %(outputDir, prefixMatch, frameNum)
        frameCounter = 0
        frameData = ''
        with open(trajFileName) as fo:
            for line in fo:
                if frameCounter == frameNum:
                    frameData += line
                if 'END' in line.strip():
                    frameCounter += 1
        
        with open(outputFileName,'wb') as of:
            of.write(frameData)
            '''            
    def generate_difference_maps(self, list_of_clusters, clusterInd2CentFileName, outputPrefix):
        print "Generating difference maps"
        #nFrames = len(frame_assignments)
        nFrames = sum([len(i) for i in list_of_clusters])
        #list_of_clusters = self.separate_clusters(frame_assignments)
        #nClusters = len(list_of_clusters)
        allFrameCounts = {}
        clusterCounts = [] 
        for clusterIndex, cluster in enumerate(list_of_clusters):
            nClusterFrames = len(cluster)
            thisClusterCounts = {}
            for matrixIndex in cluster:
                npyFilename = self.indexToNpyFile[matrixIndex]
                points = numpy.load(npyFilename)
                if len(points) == 0:
                    continue
                # If the list has intensity values
                if points.shape[1]==4:
                    for point in points:
                        tuplePoint = tuple(point[:3])
                        allFrameCounts[tuplePoint] = allFrameCounts.get(tuplePoint,0) + (point[3]/nFrames)
                        thisClusterCounts[tuplePoint] = thisClusterCounts.get(tuplePoint,0) + (point[3]/nClusterFrames)
                else:
                    for point in points:
                        tuplePoint = tuple(point)
                        allFrameCounts[tuplePoint] = allFrameCounts.get(tuplePoint,0) + (1.0/nFrames)
                        thisClusterCounts[tuplePoint] = thisClusterCounts.get(tuplePoint,0) + (1.0/nClusterFrames)
                        
            clusterCounts.append(thisClusterCounts)
            
        allPoints = numpy.array(allFrameCounts.keys())
        allFrameMap = peel.featureMap.fromPovmeList(allPoints, justCoords = True, skinDistance=2.)
        allFrameMap.data[:] = 0.0
        
        for point in allFrameCounts.keys():
            thisIndex = allFrameMap.point_to_nearest_index(point)
            allFrameMap.data[thisIndex] = allFrameCounts[point]
        
        clusterMaps = []
        for thisClusterCounts in clusterCounts:
            thisClusterMap = peel.featureMap.fromPovmeList(allPoints, justCoords = True, skinDistance=2.)
            thisClusterMap.data[:] = 0.0
            for point in thisClusterCounts.keys():
                thisIndex = allFrameMap.point_to_nearest_index(point)
                thisClusterMap.data[thisIndex] = thisClusterCounts[point]
            
            clusterMaps.append(thisClusterMap)

        templateLoadDifference = '''mol new {%s} waitfor all
display projection Orthographic
mol modstyle 0 top Isosurface 0.7500000 0 0 1 1 1
#mol modstyle 0 !MOLID! Isosurface 0.2500000 0 1 1 1 1
#white
mol modcolor 0 top ColorID 8 

mol addfile {%s} waitfor all
mol addrep top
#mol addrep !MOLID!
mol modstyle 1 top Isosurface 0.250000 1 0 1 1 1
#mol modstyle 1 !MOLID! Isosurface 0.250000 1 2 1 1 1
#blue
mol modcolor 1 top ColorID 0 
mol showrep top 1 0 
#mol showrep !MOLID! 1 0 

mol addfile {%s} waitfor all
mol addrep top
#mol addrep !MOLID!
mol modmaterial 2 top Transparent
mol modstyle 2 top Isosurface 0.2500000 2 0 0 1 1
#mol modstyle 2 !MOLID! Isosurface 0.2500000 2 2 1 1 1
#green
mol modcolor 2 top ColorID 12

mol addfile {%s} waitfor all
mol addrep top
#mol addrep !MOLID!
mol modmaterial 3 top Transparent
mol modstyle 3 top Isosurface -0.7500000 3 0 0 1 1
#mol modstyle 3 !MOLID! Isosurface -0.2500000 3 2 1 1 1
#red
mol modcolor 3 top ColorID 1

# Now load the protein
mol addfile {%s} type {pdb} first 0 last -1 step 1 waitfor all
mol modstyle 4 top NewCartoon 0.300000 10.000000 4.100000 0



'''

        plotAll = ''
        
        allFrameDxName = '%saveragePocket.dx' %(outputPrefix)
        allFrameMap.write_dx_file(allFrameDxName)

        for clusterIndex, clusterMap in enumerate(clusterMaps):
            outputDir = '%scluster%i' %(outputPrefix, clusterIndex)
            if not os.path.exists(outputDir):
                os.system('mkdir %s' %(outputDir))
            thisClusterDxName = '%saverage.dx' %(outputPrefix)
            clusterMap.write_dx_file(outputDir+'/'+thisClusterDxName)
            
            differenceMap = copy.deepcopy(clusterMap)
            differenceMap.data = differenceMap.data - allFrameMap.data

            thisDifferenceDxName = '%sdifference.dx' %(outputPrefix)
            differenceMap.write_dx_file(outputDir+'/'+thisDifferenceDxName)
        
            thisCentroidPdbName = clusterInd2CentFileName[clusterIndex]
            
            thisVmdScript = templateLoadDifference %('../'+allFrameDxName, 
                                                     thisClusterDxName, 
                                                     thisDifferenceDxName, 
                                                     thisDifferenceDxName,
                                                     thisCentroidPdbName)
            thisVmdScript = thisVmdScript.replace('!MOLID!', '0')

            with open('%s/visualize.vmd' %(outputDir),'wb') as of:
                of.write(thisVmdScript)
            

            plotAllContrib = templateLoadDifference %(allFrameDxName, 
                                                      outputDir+'/'+thisClusterDxName, 
                                                      outputDir+'/'+thisDifferenceDxName, 
                                                      outputDir+'/'+thisDifferenceDxName,
                                                      outputDir+'/'+thisCentroidPdbName)
            plotAllContrib = plotAllContrib.replace('!MOLID!',str(clusterIndex))

            plotAll += plotAllContrib
        
        with open('%svisualizeAll.vmd' %(outputPrefix),'wb') as of:
            of.write(plotAll)

        ## Write gobstopper view script
        gobstopperViewHeader = '''

display projection Orthographic
color Display Background white
material add copy RTChrome
material change ambient Material23 0.00000
material change diffuse Material23 1.00000
material change specular Material23 0.00000
material change shininess Material23 0.00000
material change mirror Material23 0.00000
material change opacity Material23 0.00000
material change outline Material23 0.00000
material change outlinewidth Material23 0.00000
material change transmode Material23 1.00000


'''
        templateGobstopperViewScript = '''

mol new {!!!CLUSTER AVERAGE FILENAME!!!} waitfor all
mol modstyle 0 top Isosurface 0.7500000 0 0 1 1 1
#mol modstyle 0 !MOLID! Isosurface 0.2500000 0 1 1 1 1
#white
mol modcolor 0 top ColorID 8

mol addrep top
mol modstyle 1 top Isosurface 0.250000 1 0 1 1 1
#blue
mol modcolor 1 top ColorID 0 
mol showrep top 1 0 


'''
            
    # cluster parameter is a list of lists - each list is a cluster
    # Return values:
    # spread 
    # cnum - a value indicating the number of clusters excluding clusters with only 1 frame
    def find_spreads(self,clusters):

#        print cluster
        spreads = []
        cnum = 0
        newWay = True
        #print 'AAAA'
        if newWay == True:
            for cluster in clusters:
                
                thisSpread = numpy.sum(self.one_minus_overlap_values[cluster,:][:,cluster])/2
                spreads.append(thisSpread)
                
        else:
            for current_set in clusters:
                curr_spread = 0
    #            print current_set
                # All combinations of frames in current cluster
                #print current_set
                combinations = itertools.combinations(current_set,2)
        
                #NEED TO USE FIRST AND SECOND VALUES OUT OF ITERATOR!!!!!!!! 
                for frame1, frame2 in combinations:
                    curr_spread += self.one_minus_overlap_values[frame1][frame2]
                
                # Calculate the N(N-1)/2 denominator for the spread of a cluster
                # SET SPREAD TO 1 IF THERE ARE NO ELEMENTS OR ONLY A SINGLE ELEMENT 
                # IN THE CLUSTER
                if len(current_set) <= 1:
                    spreads.append(0)
                else:
                    cnum += 1
                    curr_spread /= (len(current_set)*(len(current_set)-1)/2)
                    spreads.append(curr_spread)
                

        ## Unexpected numbers of clusters are now handled in the Kelley penalty code segment
        #return spread,cnum
        return spreads
    
'''        
    # DO ALL CALCULATIONS IN ONE METHOD OR DO SEPARATE METHODS??? 
    def average_spread(self,spread,cnum):
        avg_value = 0
        for i in spread:
            avg_value += i
        
        avg_value = avg_value / cnum
        return avg_value
    
    def norm_avg_spread(self, list_avg_spread):
        max_avg_spread = list_avg_spread[0]
        min_avg_spread = list_avg_spread[0]
        
        for i in list_avg_spread:
            if i > max_avg_spread:
                max_avg_spread = i
            elif i < min_avg_spread:
                min_avg_spread = i
        
        
        return
'''

#def print_help():
#    print "To run cluster.py: cluster.py [optional -h -k] binding_overlap_file pdb_trajectory_file output_file_names"

class main():

    def __init__(self,argv):
        '''    TEMP INPUT: clsuter.py overlap_file original_traj_file '''

        ''' Pick between running kmeans or hierarchical clustering or both
            Parse the command line inputs.
        '''

        parser = argparse.ArgumentParser(description="Cluster POVME pocket volumes.")
        parser.add_argument('-m', 
                            help='The pocket overlap matrix (generated by binding_site_overlap.py)')
        parser.add_argument('-t', nargs='?', action='append',default=[],
                            help='A mapping between .npy file prefixes and their original pdb file/trajectory, separated by a colon. Can be used repeatedly. Ex: -t 1BYQ:./trajectories/1BYQ.pdb -t 1UYF:./trajectories/1UYF.pdb')
        parser.add_argument('-T', nargs='?', action='append', default=[],
                            help='A file containing a series of -t arguments')
        parser.add_argument('-i', nargs='?', 
                            help='The index file mapping the pocket overlap matrix to frame numbers. Required to return cluster representatives.')
        parser.add_argument('--kmeans', action='store_true', 
                            help='Use kmeans clustering instead of hierarchical.')
        parser.add_argument('-n', nargs='?', 
                            help='Manually set number of clusters. Otherwise the Kelley penalty will calculate the optimal number.')
        parser.add_argument('-N', nargs='?', 
                            help='Set min, min:max, or min:max:skip values for number of clusters that the Kelley penalty can consider.')
        parser.add_argument('-o', nargs='?', default='',
                            help='Specify an output prefix.')
        parser.add_argument('-a', action='store_true',
                            help='Output all frames into cluster subdirectories (not just cluster reps).')
        
        args = parser.parse_args(sys.argv[1:])
        #''' Initial options '''
        #command_input = {}
        #command_input['kmeans'] = False
        #command_input['hierarchical'] = False
        #command_input['csv_file'] = ''
        #command_input['output_name'] = ''
        #command_input['num_clusters'] = None
        #command_input['indexToFrames'] = ''

        #'''Quick and dirty hack for options - need to find more efficient way '''

        #for arg in argv:
        #    print arg
        #    if 'indexMapToFrames.csv' in arg:
        #        command_input['indexToFrames'] = arg
        #    elif '.csv' in arg:
        #        command_input['csv_file'] = arg
        #    elif arg == "-k":
        #        command_input['kmeans']= True
        #    elif arg == "-h":
        #        command_input['hierarchical'] = True
        #    elif arg.isdigit():
        #        command_input['num_clusters'] = int(arg)

        # Print message and exit out of program if missing essential files
        if args.m == '':
            print "Cannot run cluster.py: Need an overlap file from binding_site_overlap.py in order to cluster \n"
            #print_help()
            sys.exit(1)

        
        #if command_input['indexToFrames'] =='':
        #print args.i
 

        
#        if command_input['pdb_file'] == '':
#            print "Cannot run cluster.py: Need the initial trajectory file in order to separate into clusters \n"
#            print_help()
#            sys.exit(1)
                
        ''' Currently only matches correctly if you run script from the folder where all the npy files are located '''
#        for filename in os.listdir('.'):
#            print filename
#            if fnmatch.fnmatch(filename,input_string_file):
#                command_input['traj_file'].append(filename)

        # If both -k and -h weren't specified, then we want to allow both options
        ##if command_input['kmeans'] == command_input['hierarchical']:
        ##    command_input['kmeans'] = True
        ##    command_input['hierarchical'] = True

#        command_input['output_name'] = command_input['pdb_file'].strip('.pdb')

        # Read csv overlap file and parse results
        csv_input = InputReader()
        #csv_input.read_traj(command_input['indexToFrames'])
        #csv_input.read_overlap(command_input['csv_file'])

        csv_input.read_overlap(args.m)
        #if args.i != None:
        csv_input.read_indexFile(args.i)
        csv_input.parse_traj_inputs(args.t, args.T)

        #print csv_input.prefixToTrajectory
        #1/0
        #else:
            #If the user didn't specify an index file, make a dictionary that just returns the input number
            #nFrames = len(csv_input.overlapMatrix)
            #for i in range(nFrames):
            #    csv_input.frameToFileName[i]=i

        #coordinates = Cluster(csv_input.coordinates,csv_input.frames,csv_input.overlapMatrix,csv_input.frameToFileName)
        clustering_obj = Cluster(csv_input)
        #print args.n
        if args.n != None:
            if args.N != None:
                raise Exception('Both -n and -N command line options specified.')
            #If the user manually specified a number of clusters
            if args.kmeans == True:
                list_of_clusters = clustering_obj.kmeans_cluster(int(args.n))
            else:
                list_of_clusters = clustering_obj.hierarchical_cluster(int(args.n))
            #clusters = clustering_obj.separate_clusters(frame_assignments)
        # If the user didn't specify the number of clusters...
        else:
            # ...use the kelley penalty to find the optimal number
            if args.N != None:
                argsNsp = args.N.split(':')
                if len(argsNsp) == 1:
                    maxKPClusters = int(argsNsp[0])
                    userNClusters = range(1,maxKPClusters+1)
                    print "Computing Kelley penalty for nClusters from 1 to %i" %(maxKPClusters)
                elif args.N.count(':') == 1:
                    minKPClusters = int(argsNsp[0])
                    maxKPClusters = int(argsNsp[1])
                    userNClusters = range(minKPClusters,maxKPClusters+1)
                    print "Computing Kelley penalty for nClusters from %i to %i" %(minKPClusters, maxKPClusters)
                elif args.N.count(':') == 2:
                    minKPClusters = int(argsNsp[0])
                    maxKPClusters = int(argsNsp[1])
                    stride = int(argsNsp[2])
                    userNClusters = range(minKPClusters,maxKPClusters+1,stride)
                    print "Computing Kelley penalty for nClusters from %i to %i, taking strides of %i" %(minKPClusters, maxKPClusters, stride)
                    
            else:
                #maxKPClusters = clustering_obj.frames
                maxKPClusters = min(75, clustering_obj.frames)
                userNClusters = range(1,maxKPClusters+1)
                print "Computing Kelley penalty for nClusters from 1 to %i" %(maxKPClusters)
                
                
            ## In order to achieve proper scaling, we must ALWAYS have 1 as a possible cluster number 
            ## in the Kelley penalty computations
            if not 1 in userNClusters:
                potentialNClusters = [1] + userNClusters
            else:
                potentialNClusters = userNClusters
            
            ## We'll be doing numpy-style slicing later so convert it here
            potentialNClusters = numpy.array(potentialNClusters)
                                             
            clustering_results = {}
            #avSpreads = numpy.zeros(maxKPClusters+1)
            avSpreads = []
            # Invalidate "0 clusters" option
            #avSpreads[0] = -1
            
            
            lastProgressDecile = 0
            for index, nClusters in enumerate(potentialNClusters):
                progressDecile = (10*index)/len(potentialNClusters)
                if progressDecile > lastProgressDecile:
                    lastProgressDecile = progressDecile
                    print "Kelley penalty " + str(progressDecile*10) + "% computed"
                if args.kmeans == True:
                    list_of_clusters = clustering_obj.kmeans_cluster(nClusters)
                else:
                    list_of_clusters = clustering_obj.hierarchical_cluster(nClusters)
                
                    
                clustering_results[nClusters] = list_of_clusters
                #clusters = clustering_obj.separate_clusters(frame_assignments)
                if len(list_of_clusters) != nClusters:
                    # If we didn't get as many clusters as we expected, put a placeholder in the array
                    #avSpreads[nClusters] = -1
                    avSpreads.append(-1)
                    # and then skip to the next iteration
                    continue
                
                spreads = clustering_obj.find_spreads(list_of_clusters)
                nSingletonClusters = sum([len(i) == 1 for i in list_of_clusters])
                avSpread = float(sum(spreads)) / (nClusters - nSingletonClusters)
                #avSpreads[nClusters] = avSpread
                avSpreads.append(avSpread)
                
            avSpreads = numpy.array(avSpreads)
            
            ## Remove places where the spread is -1
            # Make boolean index array validSpreads (eg; [0 1 1 0 1] )
            validSpreads = avSpreads != -1
            # find indices of valid spreads (eg; [1 2 4] for above)
            validIndices = numpy.nonzero(validSpreads)[0]

            ## Remove invalid nClusters and avSpread values
            validNClusters = potentialNClusters[validIndices]
            avSpreads = avSpreads[validSpreads]
            ## Now normalize spreads to the range (1, N-1)
            # subtract to bring the minimum value to 0
            avSpreads -= numpy.min(avSpreads)
            # multiply to bring the max value to N-2
            avSpreads *= (clustering_obj.frames-2)/numpy.max(avSpreads)
            # Then add 1 to everything to shift the range to (1, N-1)
            avSpreads += 1
            
            ## and finally compute the penalty value
            kPenalties = avSpreads + (1.0*validNClusters)
            
            #pylab.scatter(validNClusters, kPenalties)
            if 1 in userNClusters:
                pylab.plot(validNClusters, kPenalties, '-o')
            else:
                pylab.plot(validNClusters[1:], kPenalties[1:], '-o')
            pylab.show()
            
            optimal_nClusters = validNClusters[numpy.argsort(kPenalties)[0]]
            list_of_clusters = clustering_results[optimal_nClusters]
            #clusters = clustering_obj.separate_clusters(frame_assignments)
            print "Done computing Kelley penalty. Optimal number of clusters is %i" %(optimal_nClusters)

        clusterInd2CentFileName = clustering_obj.find_centroids(list_of_clusters, args.o)

        ## If desired, output all frames in this cluster
        if args.a == True:
            clustering_obj.outputAllFrames(list_of_clusters, args.o)
        
        ## Generate cluster characteristics
        if (clustering_obj.indexToFrame != {}) and (clustering_obj.prefixToTrajectory != {}):
            clustering_obj.generate_difference_maps(list_of_clusters, clusterInd2CentFileName, args.o)
        


if __name__ == "__main__": main(sys.argv)




















