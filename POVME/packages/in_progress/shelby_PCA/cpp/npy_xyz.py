#!/usr/bin/python

#usage: ./exec frame1.npy frame2.npy ... frameN.npy 

import numpy as np

import sys

frame_count=len(sys.argv)

coordinate_set=set([])
unique_coordinate_set=set([])
common_coordinate_set=set([])

unique_filename='unique_pocket_coordinate'
nx3_filename='pocket_coordinates'

for x in range(1, frame_count):
	frame_coords_xyz_filename=str(sys.argv[x])
	frame_coords_xyz_filename=frame_coords_xyz_filename[:-4].strip()
	frame_coords_xyz_filename=frame_coords_xyz_filename+'.xyz_arr'
	print(frame_coords_xyz_filename)

	frame_pocket_coordinates=np.load(str(sys.argv[x]))
	frame_pocket_coordinates=set([tuple(i) for i in frame_pocket_coordinates])

	frame_pocket_coordinates_list=list(frame_pocket_coordinates)
	frame_pocket_coordinates_string_list=[[str(i) for i in j] for j in frame_pocket_coordinates_list]
	write_frame_xyz=open(frame_coords_xyz_filename, 'w')
	print('\n'.join(['\t'.join(i) for i in frame_pocket_coordinates_string_list]), file=write_frame_xyz)
	
	coordinate_set=(frame_pocket_coordinates | coordinate_set)
	print('atoms in frame', x, len(frame_pocket_coordinates))

	common_coordinate_set=(frame_pocket_coordinates - common_coordinate_set)

print('total atoms', len(coordinate_set))

total_atom_list=list(coordinate_set)
total_atom_string_list=[[str(i) for i in j] for j in total_atom_list]

write_file=open(nx3_filename, 'w')
print('\n'.join(['\t'.join(i) for i in total_atom_string_list]), file=write_file)

unique_coordinate_set=coordinate_set-common_coordinate_set

print('unique atoms', len(unique_coordinate_set))

unique_atom_list=list(unique_coordinate_set)
unique_atom_string_list=[[str(i) for i in j] for j in unique_atom_list]

write_file_2=open(unique_filename, 'w')
print('\n'.join(['\t'.join(i) for i in unique_atom_string_list]), file=write_file_2)

