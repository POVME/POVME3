#!/usr/bin/python

#Authors: Shelby Straight
#
#Last Modified: 10 Jan 2015
#
#Description: Python implementation of PCA program to be run on POVME data
#
#usage: ./exec frame1.npy frame2.npy ... frameN.npy 

import numpy as np
from numpy import linalg as LA_num

import scipy
from scipy import spatial
from scipy import linalg as LA_sci

import sys

frame_count=len(sys.argv)

coordinate_set=set([])
common_coordinate_set=set([])

nx3_filename='common_points.xyz'

#first loop is necessary to built the common_coordinate_set; this is the data
#structure which lists all of the points identified by POVME as part of the 
#binding pocket in the N frames given to the PCA program

for x in range(1, frame_count):

	frame_pocket_coordinates=np.load(str(sys.argv[x]))
	frame_pocket_coordinates=set([tuple(i) for i in frame_pocket_coordinates])

	frame_pocket_coordinates_list=list(frame_pocket_coordinates)
	frame_pocket_coordinates_string_list=[[str(i) for i in j] for j in frame_pocket_coordinates_list]
	
	coordinate_set=(frame_pocket_coordinates | coordinate_set)
	print 'atoms in frame', x, len(frame_pocket_coordinates)

	common_coordinate_set=(frame_pocket_coordinates - common_coordinate_set)

print 'Total number of points identified as part of the pocket: ', len(coordinate_set)
total_atoms=len(coordinate_set)

coordinates = np.array([list(i) for i in coordinate_set])

#these common coordinates will need to be print in .xyz format to the file common_pointx.xyz
#this is required to build pdb's of these PCs downstream

write_common_coordinate_filestream = open(nx3_filename, 'w')
print >> write_common_coordinate_filestream, total_atoms
print >> write_common_coordinate_filestream, 'Common points xyz file'

print coordinates[total_atoms-1,0],coordinates[total_atoms-1,1],coordinates[total_atoms-1,2]

for x in range(0,total_atoms):
	print >> write_common_coordinate_filestream, 'O ', coordinates[x,0], coordinates[x,1], coordinates[x,2]

#now that all the common points have been identified, pairwise distances can be
#computed between all points in the data set to identify the resolution (in
#Angstroms) at which POVME was run.
#
#This is necessary to build the 3D arrays/look-up tables necessary to perform 
#"binary presence" analysis 

resolutions = scipy.spatial.distance.pdist(coordinates, 'euclidean')
reso = np.min(resolutions)

#With the resolution calculated, we calculate the shape of the 3D structure we 
#will need to store our binary point presence data in each of the three dimensions

minx = np.min(coordinates[:,0])
maxx = np.max(coordinates[:,0])
xrange = maxx-minx
xshape = xrange/reso

miny = np.min(coordinates[:,1])
maxy = np.max(coordinates[:,1])
yrange = maxy-miny
yshape = yrange/reso

minz = np.min(coordinates[:,2])
maxz = np.max(coordinates[:,2])
zrange = maxz-minz
zshape = zrange/reso

#Initialization of 3D and complimentary 1D binary data structures. Both the
#"average presence" of point i (with coordinates x_i,y_i,z_i) and the immediate
#presence of point i during frame t-eg, presence_of_point(i)=function_of(x,y,z,t)-
#must be used to compute the covariance matrix

average_point_presence = np.zeros((xshape+1, yshape+1, zshape+1))
average_presence = np.zeros((total_atoms+1))
binary_frame_arrays = np.zeros((xshape+1,yshape+1,zshape+1,frame_count+1))
binary_frame_array = np.zeros((total_atoms+1, frame_count+1))

#compute the "average binary presence" of each point in the array, save
#per-frame presence in the binary_frame_arrays data structure

print 'computing average and per-frame presences'

for x in range(1,frame_count):

	frame_pocket_coordinates=np.load(str(sys.argv[x]))
	frame_pocket_coordinates=set([tuple(i) for i in frame_pocket_coordinates])

	frame_coordinates=np.array([list(i) for i in frame_pocket_coordinates])

	for point in frame_coordinates:
		
		#rescale the frame_point with respect to the shape of the data structure
		rel_frame_point = point - np.array([minx,miny,minz])
		rel_frame_point_scaled = rel_frame_point / reso
		rel_frame_point_scaled = np.array(rel_frame_point_scaled, dtype=int)

		#each frame contributes the same to the "average presence" of the point
		#identified as a member of that frame's "frame_coordinates" temp data structure
		contribution = 1.0/frame_count

		#avg_point_presence is initially zeroes, increment it with respect to the avg_contribution
		average_point_presence[rel_frame_point_scaled[0],rel_frame_point_scaled[1],rel_frame_point_scaled[2]] += contribution

		#if the frame's point is in the list, change the binary value of that point in the corresponding
		#3D data structure from zero to one
		binary_frame_arrays[rel_frame_point_scaled[0],rel_frame_point_scaled[1],rel_frame_point_scaled[2],x]=1

print 'converting data from 3D arrays to vectors'

#loop responsible for transferring data from the 3D average presence array to a 1D array 
#This should reduce look-up time when computing the covariance matrix
point_counter = 0
for point in coordinates:

	rel_master_point = point - np.array([minx,miny,minz])
	rel_master_point_scaled = rel_master_point / reso
        rel_master_point_scaled = np.array(rel_master_point_scaled, dtype=int)

	average_presence[point_counter]=average_point_presence[rel_master_point_scaled[0],rel_master_point_scaled[1],rel_master_point_scaled[2]]

	for x in range(1,frame_count):
		binary_frame_array[point_counter,x] = binary_frame_arrays[rel_master_point_scaled[0],rel_master_point_scaled[1],rel_master_point_scaled[2],x]		

	point_counter += 1

#Computation of the covariance matrix
denominator = frame_count - 1		
cov_mat=np.zeros((total_atoms+1,total_atoms+1))

print 'computing covariance matrix'
#This algorithm is take from (REF HERE)
#As is, it takes ~2-3 minutes on a standard linux workstation

for i in range(0,total_atoms):

	avg_point_1 = average_presence[i]

	for j in range(0,total_atoms):

		covariance_entry = 0.0

		avg_point_2 = average_presence[j]
		
		for x in range(1,frame_count):
			
			point_1_variance = binary_frame_array[i,x] - avg_point_1
			point_2_variance = binary_frame_array[j,x] - avg_point_2

			numerator = point_1_variance * point_2_variance

			covariance_entry += (numerator/denominator)

		cov_mat[i, j]=covariance_entry

print 'Covariance matrix computed'
print 'NumPy is now diagonalizing the covariance matrix, which is of dimenstion ', total_atoms 
print 'Please be patient...'

#Using eigh over eig dramatically improves performance: eigenvectors are now computed
#in about 2 minutes because of the symmetry of the matrix, vs 1hr from straight schoolbook
#calculation

eigenvalues, eigenvectors = LA_num.eigh(cov_mat)

#Write the scree_plot data
scree_plot_filename = 'PC_screeplot.dat'
scree_plot_write_stream = open(scree_plot_filename, 'w')

normalization = sum(eigenvalues)

#print to the scree plot file
for count in range(0, 9):
	#it seems as if the eigenvalues are stored in increasing order; thus the index
	#must be reversed during look-up
	print >> scree_plot_write_stream, eigenvalues[(total_atoms-count)]/normalization

#print principal components as a vector, downstream tcl script turns this
#into a .pdb file
for x in range(0,9):

	#here, the eigenvector is stored as a column vector in the matrix "eigenvectors[i,j]"
	#Keep j fixed and loop through i to print each eigenvector to it's column

	#keep in mind, however, that the eigenvectors with the largest eigenvalues are stored
	#furthest "to the right" in the eigenvector data structure

	eigenvector_filename='prin_comp_'
	eigenvector_filename += str(x)
	write_file=open(eigenvector_filename, 'w')
	for y in range(0,total_atoms):
		print >> write_file, np.real(eigenvectors[y,(total_atoms-x)])
	

