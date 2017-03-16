#!python

# Calculates the binding site overlap between sets of POVME outputs.
# Started July 9th, 2014
# Celia Wong
# Advised by Jeff Wagner
# Amaro Lab, UCSD

import numpy
import sys
import re
import os
import csv
from optparse import OptionParser

class Trajectory():

    def __init__(self):
        self.coordinates = []
        self.aromatic_coordinates = []
        self.hbondAcceptor_coordinates = []
        self.hbondDonor_coordinates = []
        self.frames = 0
        self.frameToFileName = {}
        
        self.volumetric_filename = []
        self.aromatic_filename = []
        self.hbondAcceptor_filename = []
        self.hbondDonor_filename = []
        

    #traj_file is a list containing all the filenames of npy files found
    #Must run colored vs uncolored from different directories else uncolored will be reading too many files
     
    #2 pass run, first pass for volumentric only
    #Second pass allow for other wanted files assuming the color flag is set (aromatic, hbondAcceptor, hbondDonor)
    def read_traj(self,traj_file, color):
                
        #check regex for volumetric npy files
        expr = re.compile('frame_\d*.npy')
        #frameNo = re.findall(traj_file, 'frame_([0-9]+).npy')
        #frameNo = re.findall(traj_file, 'frame_([0-9]+)_aromatic.npy')
        sortedFramesAndNames = [(int(re.findall('frame_([0-9]+).npy',name)[0]), name) for name in traj_file]
        sortedFramesAndNames.sort(key=lambda x:x[0])
        sortedNames = [i[1] for i in sortedFramesAndNames]
        print sortedNames
        self.frames = len(sortedNames)
        for index, fileName in enumerate(sortedNames):
            self.volumetric_filename.append(fileName)
            self.coordinates.append(set([tuple(dummy_atom) for dummy_atom in numpy.load(fileName)]))
            self.frameToFileName[index] = fileName
            
        #count = 0
        #for i in range(len(traj_file)):
        #    match_value = expr.search(traj_file[i])
        #    if match_value:
        #        self.volumetric_filename.append(traj_file[i])
        #        self.frames += 1
        #        self.coordinates.append(set([tuple(dummy_atom) for dummy_atom in numpy.load(traj_file[i])]))
        #        self.frameToFileName[count] = traj_file[i]
        #        count += 1
                
                
        #print self.volumetric_filename
        
        #Regex to find aromatic, hbondAcceptor, and hbondDonor files
        #This should guarantee that the frames are in the same order for volumetric and all color metrics
        
        if (color):
            for i in range(len(self.volumetric_filename)):
                aromatic_file = self.volumetric_filename[i].strip('.npy')+'_aromatic.npy'
                hbondAcc_file = self.volumetric_filename[i].strip('.npy')+'_hbondAcceptor.npy'
                hbondDon_file = self.volumetric_filename[i].strip('.npy')+'_hbondDonor.npy'
                self.aromatic_filename.append(aromatic_file)
                self.aromatic_coordinates.append(set([tuple(dummy_atom[:3]) for dummy_atom in numpy.load(aromatic_file) if dummy_atom[3] > 0.02]))
                self.hbondAcceptor_filename.append(hbondAcc_file)
                self.hbondAcceptor_coordinates.append(set([tuple(dummy_atom[:3]) for dummy_atom in numpy.load(hbondAcc_file) if dummy_atom[3] > 0.02]))
                self.hbondDonor_filename.append(hbondDon_file)
                self.hbondDonor_coordinates.append(set([tuple(dummy_atom[:3]) for dummy_atom in numpy.load(hbondDon_file) if dummy_atom[3] > 0.02]))
        

class Overlap():

    def __init__(self,coordinates):
        self.coordinates = coordinates
        self.volumes = []
        
        for i in range(len(self.coordinates)):
            self.volumes.append(len(self.coordinates[i]))

    def sum_volumes(self,frame1,frame2):
        total = self.volumes[frame1] + self.volumes[frame2]
        return total

    def number_overlap(self,frame1,frame2):
        # must be a set data structure in order to use
        
        setFrame1 = self.coordinates[frame1]
        setFrame2 = self.coordinates[frame2]
        
        if not isinstance(setFrame1,set) or not isinstance(setFrame2,set):
            print "Coordinates must be contained in a set object"
            sys.exit(1)

        num_overlap_points = len(setFrame1.intersection(setFrame2))
        
        '''Test here for error values '''
        if num_overlap_points > min(len(setFrame1),len(setFrame2)):
            print 'invalid overlap value'
        return num_overlap_points

    # need to calculate the volume that is overlapped
    def volume_overlap(self,frame1,frame2):
        num_overlap_points = self.number_overlap(frame1, frame2)

        overlap_ratio = num_overlap_points / float(len(self.coordinates[frame1]))

        vol_overlap = overlap_ratio * self.volumes[frame1]

        return vol_overlap
    
    def tanimoto_overlap (self,vlap,frame1,frame2):
        #vlap = self.number_overlap(frame1,frame2)
        vtotal = self.sum_volumes(frame1,frame2) - vlap
        tanimoto = float(vlap/float(vtotal))
        return tanimoto

    def tversky_overlap(self,vlap,frame1,frame2):
        ''' Tversky index: overlap/[Va(w/o overlap) + overlap] & overlap/[Vb(w/o overlap) + overlap]
        '''
        #volume_overlap = self.number_overlap(frame1,frame2)
        total_volume = self.volumes[frame1]
        #print 'volume_overlap = {0}, total volume = {1}'.format(volume_overlap,total_volume)
        tversky = float(vlap/float(total_volume))
        return tversky

''' Helper function needed to get all the input files from optparse '''
def vararg_callback(option, opt_str, value, parser):
    assert value is None
    value = []
        
    for arg in parser.rargs:
        if ((arg[:1] == "-" or arg[:2] == "--") and len(arg) > 1):
            break
        value.append(arg)
            
    del parser.rargs[:len(value)]
    setattr(parser.values, option.dest, value)

class main():

    def __init__(self,argv):
        
        if len(argv) == 1:
            print "Cannot run binding_site_overlap.py: Need to specify .npy files to be read"
            sys.exit(1)
        parser = OptionParser()
        
        parser.add_option("-f", dest = "filename", action = "callback", callback = vararg_callback, help = "All files from POVME that you want to run similarity calculations on.")
        parser.add_option("-c", action = "store_true", dest="color", help = "run similarity calculations on colored output from POVME")
        parser.add_option("--csv", action = "store_true", dest="csv", help = "save human-readable CSV distance matrices.")
        (options, args) = parser.parse_args(argv)
        
        #print command_input['traj_file']
        
        file_input = Trajectory()
        #file_input.read_traj(argv)
        file_input.read_traj(options.filename, options.color)



        ''' Saving which index refers to which frame file for use in clustering '''
        
        frames_dict = file_input.frameToFileName
        with open('indexMapToFrames.csv','wb') as csvfile:
            fieldnames = ['index','frame']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            for i in frames_dict:
                writer.writerow({"index": i, "frame" : frames_dict[i]})
        
        overlap_value = Overlap(file_input.coordinates)
        aromatic_overlap = Overlap(file_input.aromatic_coordinates)
        hbondAcc_overlap = Overlap(file_input.hbondAcceptor_coordinates)
        hbondDon_overlap = Overlap(file_input.hbondDonor_coordinates)

        num_frames = len(file_input.coordinates)
        print "The number of frames found was: {0}".format(num_frames)
        
        '''Make a matrix of number of overlapping coordinates first'''
        
        # Always calculate the overlapping points for volumetric
        overlapStyle = 2
        num_overlap = numpy.empty([file_input.frames,file_input.frames], dtype=float)
        if overlapStyle == 1:
            for f1 in range(num_frames):
                for f2 in range(f1,num_frames):
                    num_overlap[f1,f2] = overlap_value.number_overlap(f1, f2)
        elif overlapStyle == 2:
            allPointsSet = set()
            for f1 in range(num_frames):
                for coord in file_input.coordinates[f1]:
                    allPointsSet.add(coord)
            coord2vectPos = dict([(coord, i) for i, coord in enumerate(allPointsSet)])
            nVectPos = len(allPointsSet)
            vecPosMatrix = numpy.zeros((num_frames, nVectPos), dtype=numpy.bool)
            for f1 in range(num_frames):
                for point in file_input.coordinates[f1]:
                    vecPosMatrix[f1,coord2vectPos[point]] = 1
            for f1 in range(num_frames):
                for f2 in range(f1,num_frames):
                    #raise Exception('There\'s a problem here. In comparing two frames in chris condon\'s pockets, values >1 were found in the tanimoto matrix.')
                    num_overlap[f1,f2] = numpy.count_nonzero(vecPosMatrix[f1,:] & vecPosMatrix[f2,:])
                    num_overlap[f2,f1] = num_overlap[f1,f2]
                    
        
        # Only calculate the colored option if the color option was set
        if (options.color):
            aromatic_num_overlap = numpy.empty([file_input.frames,file_input.frames], dtype=float)
            for f1 in range(num_frames):
                for f2 in range(num_frames):
                    aromatic_num_overlap[f1,f2] = aromatic_overlap.number_overlap(f1, f2)
        
        
            hbondAcc_num_overlap = numpy.empty([file_input.frames,file_input.frames], dtype=float)
            for f1 in range(num_frames):
                for f2 in range(num_frames):
                    hbondAcc_num_overlap[f1,f2] = hbondAcc_overlap.number_overlap(f1, f2)
        
        
            hbondDon_num_overlap = numpy.empty([file_input.frames,file_input.frames], dtype=float)
            for f1 in range(num_frames):
                for f2 in range(num_frames):
                    hbondDon_num_overlap[f1,f2] = hbondDon_overlap.number_overlap(f1, f2)
                    
                    
        '''Record the overlap_matrix in a csv file'''

        print "Starting Tanimoto calculations"

        tanimoto_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)
        
        if (options.color):
            aromatic_tanimoto_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)
            hbondAcc_tanimoto_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)
            hbondDon_tanimoto_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)
            colored_tanimoto_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)

        for f1 in range(num_frames):
            for f2 in range(f1,num_frames):
                tanimoto_result = overlap_value.tanimoto_overlap(num_overlap[f1,f2],f1,f2)
                tanimoto_matrix[f1,f2] = tanimoto_result
                tanimoto_matrix[f2,f1] = tanimoto_result
                if (options.color):
                    aromatic_result = aromatic_overlap.tanimoto_overlap(aromatic_num_overlap[f1,f2],f1,f2)
                    aromatic_tanimoto_matrix[f1,f2] = aromatic_result
                    aromatic_tanimoto_matrix[f2,f1] = aromatic_result
                    hbondAcc_result = hbondAcc_overlap.tanimoto_overlap(hbondAcc_num_overlap[f1,f2],f1,f2)
                    hbondAcc_tanimoto_matrix[f1,f2] = hbondAcc_result
                    hbondAcc_tanimoto_matrix[f2,f1] = hbondAcc_result
                    hbondDon_result = hbondDon_overlap.tanimoto_overlap(hbondDon_num_overlap[f1,f2],f1,f2)
                    hbondDon_tanimoto_matrix[f1,f2] = hbondDon_result
                    hbondDon_tanimoto_matrix[f2,f1] = hbondDon_result
                
                    average_similarity = (aromatic_result + hbondAcc_result + hbondDon_result + tanimoto_result)/4
                    colored_tanimoto_matrix[f1,f2] = average_similarity             
                    colored_tanimoto_matrix[f2,f1] = average_similarity
                
        print "Overlap Matrix for Tanimoto calculation"
        print tanimoto_matrix
        
        if (options.color):
            print "Aromatic"
            print aromatic_tanimoto_matrix
            print "Hydrogen Bond Acceptor"
            print hbondAcc_tanimoto_matrix
            print "Hydrogen Bond Donor"
            print hbondDon_tanimoto_matrix
            
        print "\n"
        if (options.csv):
            numpy.savetxt(os.getcwd()+'/POVME_Tanimoto_matrix.csv', tanimoto_matrix, delimiter=',')
        
        if (options.color) and (options.csv):
            numpy.savetxt(os.getcwd()+'/POVME_Tanimoto_matrix_aromatic.csv', aromatic_tanimoto_matrix, delimiter=',')
            numpy.savetxt(os.getcwd()+'/POVME_Tanimoto_matrix_hbondAcceptor.csv', hbondAcc_tanimoto_matrix, delimiter=',')
            numpy.savetxt(os.getcwd()+'/POVME_Tanimoto_matrix_hbondDonor.csv', hbondDon_tanimoto_matrix, delimiter=',')
            numpy.savetxt(os.getcwd()+'/POVME_Tanimoto_matrix_colored.csv', colored_tanimoto_matrix, delimiter=',')
            
        numpy.save('tanimoto_matrix.npy', tanimoto_matrix)
        
        if (options.color):
            numpy.save('aromatic_tanimoto_matrix.npy', aromatic_tanimoto_matrix)
            numpy.save('hbondAcc_tanimoto_matrix.npy', hbondAcc_tanimoto_matrix)
            numpy.save('hbondDon_tanimoto_matrix.npy', hbondDon_tanimoto_matrix)


        print "Starting Tversky calculations"
        tversky_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)
        if (options.color):
            aromatic_tversky_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)
            hbondAcc_tversky_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)
            hbondDon_tversky_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)
            colored_tversky_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)
        

        for f1 in range(num_frames):
            for f2 in range(num_frames):
                #print 'length of frame 1 = {0}'.format(len(overlap_value.coordinates[f1]))
                tversky_matrix[f1,f2] = overlap_value.tversky_overlap(num_overlap[f1,f2],f1, f2)
                if (options.color):
                    aromatic_tversky_matrix[f1,f2] = aromatic_overlap.tversky_overlap(aromatic_num_overlap[f1,f2], f1, f2)
                    hbondAcc_tversky_matrix[f1,f2] = hbondAcc_overlap.tversky_overlap(hbondAcc_num_overlap[f1,f2], f1, f2)
                    hbondDon_tversky_matrix[f1,f2] = hbondDon_overlap.tversky_overlap(hbondDon_num_overlap[f1,f2], f1, f2)
                
                    colored_tversky_matrix[f1,f2] = (hbondAcc_tversky_matrix[f1,f2] + hbondDon_tversky_matrix[f1,f2] + tversky_matrix[f1,f2] + aromatic_tversky_matrix[f1,f2])/4

        print "Overlap Matrix for Tversky calculation"
        print tversky_matrix
        
        if (options.color):
            print "Aromatic"
            print aromatic_tversky_matrix
            print "Hydrogen Bond Acceptor"
            print hbondAcc_tversky_matrix
            print "Hydrogen Bond Donor"
            print hbondDon_tversky_matrix
            
        print "\n"
        
        if (options.csv):
            numpy.savetxt(os.getcwd()+'/POVME_Tversky_matrix.csv',tversky_matrix,delimiter=',')
        
        if (options.color) and (options.csv):
            numpy.savetxt(os.getcwd()+'/POVME_Tversky_matrix_aromatic.csv',aromatic_tversky_matrix,delimiter=',')        
            numpy.savetxt(os.getcwd()+'/POVME_Tversky_matrix_hbondAcceptor.csv',hbondAcc_tversky_matrix,delimiter=',')        
            numpy.savetxt(os.getcwd()+'/POVME_Tversky_matrix_hbondDonor.csv',hbondDon_tversky_matrix,delimiter=',')
            numpy.savetxt(os.getcwd()+'/POVME_Tversky_matrix_colored.csv',colored_tversky_matrix,delimiter=',')
        
        
        numpy.save('tversky_matrix.npy', tversky_matrix)
        
        if (options.color):
            numpy.save('aromatic_tversky_matrix.npy', aromatic_tversky_matrix)
            numpy.save('hbondAcc_tversky_matrix.npy', hbondAcc_tversky_matrix)
            numpy.save('hbondDon_tversky_matrix.npy', hbondDon_tversky_matrix)
                
        #print "Map of index numbers to npy files"
        #print frames_dict

            
if __name__ == "__main__": main(sys.argv)
