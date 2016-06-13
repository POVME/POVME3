# Calculates the binding site overlap between different conformations of the same protein.
# Started July 9th, 2014
# Celia Wong

import numpy
import sys
import re
import os
import csv

class Trajectory():

    def __init__(self):
        self.coordinates = []
        self.frames = 0
        self.frameToFileName = {}

    # Save each frame in the trajectory as a set
    ''' def read_traj(self,traj_file):

        trajectory = open(traj_file,'r')
        frame_coordinates = []
        counter = []
        count = 0
    '''
    ''' Check first line to see if it contains volume count, if not, incorrect pdb file '''

    ''' for line in trajectory:
            if line[0] == 'E':
                self.coordinates.append(set(frame_coordinates))
                #self.coordinates.append(frame_coordinates)
                #numpy.append(self.coordinates,set(frame_coordinates))
                frame_coordinates = []

                #counter.append(count)
                #count = 0

                self.frames += 1

            elif len(line) > 6 and line[7] == 'V':
                frame_vol = float(re.sub("[^0123456789\.]","",line))
                self.volume.append(frame_vol)

            elif line[0] == 'A':
                # If there is an X present there is no viable coordinates for that line
                if line[17] != 'X':
                    #if (float(line[29:37].strip()),float(line[38:45].strip()),float(line[46:54].strip())) in frame_coordinates:
                    #   print (float(line[29:37].strip()),float(line[38:45].strip()),float(line[46:54].strip()))
                    frame_coordinates.append((float(line[29:37].strip()),float(line[38:45].strip()),float(line[46:54].strip())))
                    #count += 1

        #print counter
        print self.coordinates
        trajectory.close()
    '''


    def read_traj(self,traj_file):
        
        i = 0
        for filename in traj_file:
            self.coordinates.append(set([tuple(dummy_atom) for dummy_atom in numpy.load(filename)]))
            self.frames+= 1
            self.frameToFileName[i] = filename
            i += 1
        #print self.coordinates
        #print self.frameToFileName
        #print "\n"
        return

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
        
        #print self.coordinates[frame1]
        #setFrame1 = set([tuple(i) for i in self.coordinates[frame1]])
        #setFrame2 = set([tuple(i) for i in self.coordinates[frame2]])
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

    ''' Incorrect formula do not use
    def schrodinger_overlap(self,frame1,frame2):
        vlap = self.number_overlap(frame1,frame2)
        vtotal = self.sum_volumes(frame1,frame2)
        print vlap, vtotal
        normalized_overlap = float((2*vlap)/float(vtotal))
        return normalized_overlap
    '''
    
    def tanimoto_overlap (self,frame1,frame2):

        vlap = self.number_overlap(frame1,frame2)
        vtotal = self.sum_volumes(frame1,frame2) - vlap
        tanimoto = float(vlap/float(vtotal))
        return tanimoto

    def tversky_overlap(self,frame1,frame2):
        ''' Tversky index: overlap/[Va(w/o overlap) + overlap] & overlap/[Vb(w/o overlap) + overlap]
        '''
        volume_overlap = self.number_overlap(frame1,frame2)
        total_volume = self.volumes[frame1]
        #print 'volume_overlap = {0}, total volume = {1}'.format(volume_overlap,total_volume)
        tversky = float(volume_overlap/float(total_volume))
        return tversky

class main():

    def __init__(self,argv):

        command_input = {}
        command_input['traj_file'] = []
        command_input['output_name'] = ''

#        for arg in argv:
#            if '.pdb' in arg:
#                command_input['traj_file'] = arg
#                command_input['output_name'] = command_input['traj_file'].split('.')[0]

        for arg in argv:
            if ('.npy' in arg) and ('frame' in arg):
                command_input['traj_file'].append(arg)
        print command_input['traj_file']
        
        if command_input == []:
            print "You need to specify .npy files to be read"
            sys.exit(1)

        file_input = Trajectory()
        file_input.read_traj(command_input['traj_file'])

        overlap_value = Overlap(file_input.coordinates)

        num_frames = len(file_input.coordinates)
        print "The number of frames found was: {0}".format(num_frames)
        
        ''' TEST TEST TEST '''
        for i in range(len(file_input.coordinates)):
            print "The number of atoms in {0} is {1}".format(file_input.frameToFileName[i],len(file_input.coordinates[i]))
        #print overlap_value.volumes

        '''TESTING HERE'''
#        print 'intersection value: {0}'.format(len(set.intersection(overlap_value.coordinates[0],overlap_value.coordinates[0])))
    
#        overlap_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)

        ''' 
            Schrodinger calculation
            Find the overlap value for all possible combinations given the number of frames
        for f1 in range(num_frames):
            for f2 in range(f1,num_frames):
                normalized_overlap = overlap_value.schrodinger_overlap(f1, f2)
                overlap_matrix[f1,f2] = normalized_overlap
                overlap_matrix[f2,f1] = normalized_overlap
        print overlap_matrix
        '''
        '''Record the overlap_matrix in a csv file'''
        #print os.getcwd()
        #numpy.savetxt(os.getcwd()+'/POVME_Schrodinger_matrix_'+command_input['output_name']+'.csv', overlap_matrix, delimiter =',')

        print "Starting Tanimoto calculations"

        print "Overlap Matrix for Tanimoto calculation"
        tanimoto_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)

        for f1 in range(num_frames):
            for f2 in range(f1,num_frames):
                tanimoto_result = overlap_value.tanimoto_overlap(f1,f2)
                tanimoto_matrix[f1,f2] = tanimoto_result
                tanimoto_matrix[f2,f1] = tanimoto_result

        print tanimoto_matrix
        #numpy.savetxt(os.getcwd()+'/POVME_Tanimoto_matrix_'+command_input['output_name']+'.csv', tanimoto_matrix, delimiter=',')
        numpy.save(os.getcwd()+'/POVME_Tanimoto_matrix_'+command_input['output_name']+'.npy',tanimoto_matrix)

        print "Starting Tversky calculations"
        tversky_matrix = numpy.empty([file_input.frames,file_input.frames],dtype=float)

        for f1 in range(num_frames):
            for f2 in range(num_frames):
                #print 'length of frame 1 = {0}'.format(len(overlap_value.coordinates[f1]))
                tversky_matrix[f1,f2] = overlap_value.tversky_overlap(f1, f2)

        print tversky_matrix
        #numpy.savetxt(os.getcwd()+'/POVME_Tversky_matrix_'+command_input['output_name']+'.csv',tversky_matrix,delimiter=',')
        numpy.save(os.getcwd()+'/POVME_Tversky_matrix_'+command_input['output_name']+'.npy',tversky_matrix)
        
        
        ''' Saving which index refers to which frame file for use in clustering '''
        
        frames_dict = file_input.frameToFileName
        print frames_dict

        with open('indexMapToFrames.csv','wb') as csvfile:
            fieldnames = ['index','frame']
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            
            for i in frames_dict:
                writer.writerow({"index": i, "frame" : frames_dict[i]})
            
if __name__ == "__main__": main(sys.argv)
