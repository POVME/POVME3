# POVME 2.0.1 is released under the GNU General Public License (see http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't hesitate to contact me,
# Jacob Durrant, at jdurrant [at] ucsd [dot] edu.
#
# If you use POVME in your work, please cite Durrant, J. D., C. A. de Oliveira, et al.
#    (2011). "POVME: An algorithm for measuring binding-pocket volumes." J Mol Graph
#    Model 29(5): 773-776.

import math
import sys
import time
import numpy
import pymolecule
import gzip
import os
import shutil
import random
import multiprocessing
import platform

try: from cStringIO import StringIO
except: from StringIO import StringIO

from scipy.spatial.distance import cdist
from scipy.spatial.distance import pdist
from scipy.spatial.distance import squareform

version = "2.0.1"

def log(astr, parameters):
    '''Output POVME statements, either to the screen or to a file

    Arguments:
    astr -- The string to output.
    parameters -- The user-defined parameters.
    
    '''
    
    # Print the output to the screen.
    print astr
    
    # Save it to the output file as well.
    try:
        if parameters['CompressOutput'] == True: f = gzip.open(parameters['OutputFilenamePrefix'] + 'output.txt.gz', 'ab')
        else: f = open(parameters['OutputFilenamePrefix'] + 'output.txt', 'a')
        
        f.write(astr + "\n")
        f.close()
    except: pass

class Multithreading():
    """A class for running calculations on multiple processors"""
    
    results = []
    
    def __init__(self, inputs, num_processors, task_class):
        """Launches a calculation on multiple processors
            
        Arguments:
        inputs -- A list, containing all the input required for the calculation
        num_processors -- An integer, the requested number of processors to use
        task_class -- An class, the class governing what calculations will be run on a given thread
        
        Returns:
        Nothing, though the objects self.results list is populated with the calculation results
        
        """
        
        self.results = []
        
        if num_processors != 1 and (platform.system().upper()[:3] == "WIN" or "NT" in platform.system().upper()): # If it's windows, you can only use one processor.
            print "WARNING: Use of multiple processors is not supported in Windows. Proceeding with one processor..."
            num_processors = 1
        
        if num_processors == 1: # so just running on 1 processor, perhaps under windows
            single_thread = task_class()
            single_thread.total_num_tasks = len(inputs)

            single_thread.results = []
            for item in inputs: single_thread.value_func(item, None)
            
            self.results = single_thread.results
            
        else: # so it actually is running on multiple processors

            cpu_count = 1
            cpu_count = multiprocessing.cpu_count()
    
            # first, if num_processors <= 0, determine the number of processors to use programatically
            if num_processors <= 0: num_processors = cpu_count
    
            # reduce the number of processors if too many have been specified
            if len(inputs) < num_processors: num_processors = len(inputs)
            
            if len(inputs) == 0: # if there are no inputs, there's nothing to do.
                self.results = []
                return
            
            # now, divide the inputs into the appropriate number of processors
            inputs_divided = {}
            for t in range(num_processors): inputs_divided[t] = []
    
            for t in range(0, len(inputs), num_processors):
                for t2 in range(num_processors):
                    index = t + t2
                    if index < len(inputs): inputs_divided[t2].append(inputs[index])

            # now, run each division on its own processor
            running = multiprocessing.Value('i', num_processors)
            mutex = multiprocessing.Lock()
    
            arrays = []
            threads = []
            for i in range(num_processors):
                athread = task_class()
                athread.total_num_tasks = len(inputs)
                
                threads.append(athread)
                arrays.append(multiprocessing.Array('i',[0, 1]))
    
            results_queue = multiprocessing.Queue() # to keep track of the results
            
            processes = []
            for i in range(num_processors):
                p = multiprocessing.Process(target=threads[i].runit, args=(running, mutex, results_queue, inputs_divided[i]))
                p.start()
                processes.append(p)
    
            while running.value > 0: is_running = 0 # wait for everything to finish

            # compile all results into one list
            for thread in threads:
                chunk = results_queue.get()
                self.results.extend(chunk)
                
class MultithreadingTaskGeneral:
    """A parent class of others that governs what calculations are run on each thread"""
    
    results = []
    
    def runit(self, running, mutex, results_queue, items):
        """Launches the calculations on this thread
            
        Arguments:
        running -- A multiprocessing.Value object
        mutex -- A multiprocessing.Lock object
        results_queue -- A multiprocessing.Queue() object for storing the calculation output
        items -- A list, the input data required for the calculation
        
        """

        for item in items: self.value_func(item, results_queue)
        
        mutex.acquire()
        running.value -= 1
        mutex.release()
        results_queue.put(self.results)

    def value_func(self, item, results_queue): # so overwriting this function
        """The definition that actually does the work.

        Arguments:
        item -- A list or tuple, the input data required for the calculation
        results_queue -- A multiprocessing.Queue() object for storing the calculation output
        
        """

        # input1 = item[0]
        # input2 = item[1]
        # input3 = item[2]
        # input4 = item[3]
        # input5 = item[4]
        # input6 = item[5]
        
        # use inputs to come up with a result, some_result
        
        #self.results.append(some_result)
        
        pass

class ConvexHull():
    """A class to handle convex-hull calculations"""

    def get_seg_dict_num(self, seg_dict, seg_index):
        """seg_dict is a dictionary object that contains information about segments within the convex hull. The keys are 2x3 tuples, which represent two ends of a segment in space. The values of seg_dict are the number of times a segment has been part of a triangle, either 1 or 2. (Zero times would mean that the segment doesn't exist in the dictionary yet). This function looks up and returns the value of a seg_index from seg_dict
            
        Arguments:
        seg_dict -- the dictionary of segment 2x3 tuples as keys, integers as values
        seg_index -- the key of the dictionary member we are going to retrieve
        
        Returns:
        if seg_index exists in the keys of seg_dict, return the value. Otherwise, return 0
        
        """

        if seg_index[0][0] > seg_index[1][0]: # we want the index with the greater x-value, so we don't get identical segments in the dictionary more than once
            index = seg_index
        else:
            index = seg_index[::-1]
        
        if index in seg_dict:
            return seg_dict[index]
        else:
            return 0
    
    def increment_seg_dict(self, seg_dict, seg_index):
        """seg_dict is a dictionary object that contains information about segments within the convex hull. The keys are 2x3 tuples, which represent two ends of a segment in space. The values of seg_dict are the number of times a segment has been part of a triangle, either 1 or 2. (Zero times would mean that the segment doesn't exist in the dictionary yet). This function increments the values within seg_dict, or initiates them if they dont exist yet.
            
        Arguments:
        seg_dict -- the dictionary of segment 2x3 tuples as keys, integers as values
        seg_index -- the key of the dictionary member we are going to increment
        
        Returns:
        None: the values of seg_dict are received and modified by reference
        """

        if seg_index[0][0] > seg_index[1][0]: # we want the index with the greater x-value, so we don't get identical segments in the dictionary more than once
            index = seg_index
        else:
            index = seg_index[::-1]
        
        #"putting index:", index, "into seg_dict because", index[0][0], ">", index[1][0]  
        
        if index in seg_dict: # if the entry already exists in seg_dict
            seg_dict[index] += 1 # increment
        else:
            seg_dict[index] = 1 # initiate with a value of 1 because it now exists on a triangle
        return
    
    def gift_wrapping_3d(self, raw_points):
        """Gift wrapping for 3d convex hull
            
        Arguments:
        raw_points -- A nx3 array of points, where each row corresponds to an x,y,z point coordinate
        
        Returns:
        A convex hull represented by a list of triangles. Each triangle is a 3x3 array, where each row is an x,y,z coordinate in space. The 3 rows describe the location of the 3 corners of the triangle. Each of the 3 points are arranged so that a cross product will point outwards from the hull
        
        
        """

        n = numpy.shape(raw_points)[0] # number of points
        point1 = raw_points[0] # take the first point
        xaxis = numpy.array([1,0,0]) # create a ref vector pointing along x axis
        maxx = raw_points[0][0] # initiate highest x value
        points = [] # a list of tuples for easy dictionary lookup
        seg_dict = {} # a dictionary that contains the number of triangles a seg is in
        
        begintime = time.time()
        for i in range(n): # find the n with the largest x value
            point = tuple(raw_points[i])
            points.append(point)
            if point[0] > maxx:
                maxx = point[0]
                point1 = raw_points[i]
        #print "find max x:", time.time() - begintime
        
        best_dot = -1.0 # initiate dot relative to x-axis
        point2 = numpy.array(raw_points[1]) # initiate best segment
        
        # find first/best segment
        begintime = time.time()
        for i in range(n):
            pointi = raw_points[i]
            if numpy.array_equal(pointi, point1): continue
            diff_vec = pointi - point1
            diff_len = numpy.linalg.norm(diff_vec)
            
            test_dot = numpy.dot(diff_vec/diff_len,xaxis)
            if test_dot > best_dot:
                best_dot = test_dot
                point2 = pointi
        
        #print "find first segment:", time.time() - begintime
        point1 = tuple(point1)
        point2 = tuple(point2)
        ref_vec = xaxis

        # now find the best triangle
        triangles = []
        
        seg_list = set([(point1, point2),])
        norm_dict = {(point1,point2):xaxis}
        self.increment_seg_dict( seg_dict, (point1,point2) )
        
        counter = 0
        first_time = True
        
        begintime = time.time()
        section1 = 0.0
        section2 = 0.0
        section3 = 0.0
        while seg_list: # as long as there are unexplored edges of triangles in the hull...
            
            counter += 1
            seg = seg_list.pop() # take a segment out of the seg_list
            tuple1 = seg[0] # the two ends of the segment
            tuple2 = seg[1]
            point1 = numpy.array(seg[0])
            point2 = numpy.array(seg[1])
            result = self.get_seg_dict_num( seg_dict, (seg[0],seg[1]) ) 

            if result >= 2: # then we already have 2 triangles on this segment
                continue # forget about drawing a triangle for this seg

            ref_vec = norm_dict[(seg[0],seg[1])] # get the norm for a triangle that the segment is part of
            
            best_dot_cross = -1.0
            best_point = None
            
            for i in range(n): # look at each point
                
                pointi = raw_points[i]
                #if numpy.array_equal(pointi, point1) or numpy.array_equal(pointi, point2): continue # if we are trying one of the points that are point1 or point2
                diff_vec1 = point2 - point1
                #diff_len1 = numpy.linalg.norm(diff_vec1)
                diff_vec2 = pointi - point2
                #diff_len2 = numpy.linalg.norm(diff_vec2)
                
                #test_cross = numpy.cross(diff_vec1/diff_len1,diff_vec2/diff_len2)
                #test_cross = numpy.cross(diff_vec1,diff_vec2)
                test_cross = numpy.array([diff_vec1[1]*diff_vec2[2]-diff_vec1[2]*diff_vec2[1], diff_vec1[2]*diff_vec2[0]-diff_vec1[0]*diff_vec2[2], diff_vec1[0]*diff_vec2[1]-diff_vec1[1]*diff_vec2[0]]) # cross product
                
                test_cross_len = numpy.sqrt(test_cross[0]*test_cross[0] + test_cross[1]*test_cross[1] + test_cross[2]*test_cross[2]) #numpy.linalg.norm(test_cross) # get the norm of the cross product
                
                if test_cross_len <= 0.0: continue
                #test_cross_len_inv = 1 / test_cross_len
                test_cross = test_cross / test_cross_len
                dot_cross = numpy.dot(test_cross, ref_vec)
                #dot_cross = test_cross[0]*ref_vec[0] + test_cross[1]*ref_vec[1] + test_cross[2]*ref_vec[2]
                if dot_cross > best_dot_cross:
                    best_cross = test_cross
                    best_dot_cross = dot_cross
                    best_point = pointi
                    tuple3 = points[i]
                
            
            point3 = best_point
            
            if self.get_seg_dict_num( seg_dict, (tuple2,tuple1) ) > 2: continue
            if self.get_seg_dict_num( seg_dict, (tuple3,tuple2) ) > 2: continue
            if self.get_seg_dict_num( seg_dict, (tuple1,tuple3) ) > 2: continue
            
            # now we have a triangle from point1 -> point2 -> point3
            # must test each edge
            if first_time:
                self.increment_seg_dict( seg_dict, (tuple2,tuple1) )
                seg_list.add((tuple2, tuple1))
                norm_dict[(tuple2,tuple1)] = best_cross
            
            self.increment_seg_dict( seg_dict, (tuple3,tuple2) )
            seg_list.add((tuple3, tuple2))
            norm_dict[(tuple3,tuple2)] = best_cross
            
            self.increment_seg_dict( seg_dict, (tuple1,tuple3) )
            seg_list.add((tuple1, tuple3))
            norm_dict[(tuple1,tuple3)] = best_cross
            
            triangles.append((numpy.array(tuple1),numpy.array(tuple2),numpy.array(tuple3)))
            
            first_time = False
            
        #print "find all triangles:", time.time() - begintime
        
        #print "section1:", section1
        #print "section2:", section2
        #print "section3:", section3
        return triangles
    
    def akl_toussaint(self, points):
        """The Akl-Toussaint Heuristic. Given a set of points, this definition will create an octahedron whose corners are the extremes in x, y, and z directions. Every point within this octahedron will be removed because they are not part of the convex hull. This causes any expected running time for a convex hull algorithm to be reduced to linear time.
        
        Arguments:
        points -- An nx3 array of x,y,z coordinates
        
        Returns:
        All members of original set of points that fall outside the Akl-Toussaint octahedron
        
        """

        x_high = (-1e99,0,0); x_low = (1e99,0,0); y_high = (0,-1e99,0); y_low = (0,1e99,0); z_high = (0,0,-1e99); z_low = (0,0,1e99)
        
        
        for point in points: # find the corners of the octahedron
            if point[0] > x_high[0]: x_high = point
            if point[0] < x_low[0]: x_low = point
            if point[1] > y_high[1]: y_high = point
            if point[1] < y_low[1]: y_low = point
            if point[2] > z_high[2]: z_high = point
            if point[2] < z_low[2]: z_low = point
          
        octahedron = [ # define the triangles of the surfaces of the octahedron
        numpy.array((x_high,y_high,z_high)),
        numpy.array((x_high,z_low,y_high)),
        numpy.array((x_high,y_low,z_low)),
        numpy.array((x_high,z_high,y_low)),
        numpy.array((x_low,y_low,z_high)),
        numpy.array((x_low,z_low,y_low)),
        numpy.array((x_low,y_high,z_low)),
        numpy.array((x_low,z_high,y_high)),
        ]
        new_points = [] # everything outside of the octahedron
        for point in points: # now check to see if a point is inside or outside the octahedron
            outside = self.outside_hull(point, octahedron, epsilon=-1.0e-5)
            if outside:
                new_points.append(point)
            
        return numpy.array(new_points) # convert back to an array
    
    def outside_hull(self, our_point, triangles, epsilon=1.0e-5):
        """Given the hull as defined by a list of triangles, this definition will return whether a point is within these or not.
            
        Arguments:
        our_point -- an x,y,z array that is being tested to see whether it exists inside the hull or not
        triangles -- a list of triangles that define the hull
        epsilon -- needed for imprecisions in the floating-point operations.
        
        Returns:
        True if our_point exists outside of the hull, False otherwise
        
        """

        our_point = numpy.array(our_point) # convert it to an numpy.array
        for triangle in triangles:
            rel_point = our_point - triangle[0] # vector from triangle corner 0 to point
            vec1 = triangle[1] - triangle[0] # vector from triangle corner 0 to corner 1
            vec2 = triangle[2] - triangle[1] # vector from triangle corner 1 to corner 2
            our_cross = numpy.cross(vec1, vec2) # cross product between vec1 and vec2
            our_dot = numpy.dot(rel_point,our_cross) # dot product to determine whether cross is point inward or outward
            if numpy.dot(rel_point,our_cross) > epsilon: # if the dot is greater than 0, then its outside
                return True
            
        return False

def unique_rows(a):
    """Identifies unique points (rows) in an array of points.
    
    Arguments:
    a -- A nx3 numpy.array representing 3D points.

    Returns:
    A nx2 numpy.array containing the 3D points that are unique.

    """
    
    a[a == -0.0] = 0.0
    b = numpy.ascontiguousarray(a).view(numpy.dtype((numpy.void, a.dtype.itemsize * a.shape[1])))
    unique_a = numpy.unique(b).view(a.dtype).reshape(-1, a.shape[1])

    return unique_a

def create_pdb_line(numpy_array, index, resname, letter):
    """Create a string formatted according to the PDB standard.

    Arguments:
    numpy_array -- A 1x3 numpy.array representing a 3D point.
    letter -- A string, the atom name/chain/etc to use for the output.

    Returns:
    A string, formatted according to the PDB standard.

    """

    if len(numpy_array) == 2: numpy_array = numpy.array([numpy_array[0], numpy_array[1], 0.0])
    if numpy_array.shape == (1, 3): numpy_array = numpy_array[0]

    output = "ATOM "
    output = output + str(index % 999999).rjust(6) + letter.rjust(5) + resname.rjust(4) + letter.rjust(2) + str(index % 9999).rjust(4)
    output = output + ("%.3f" % numpy_array[0]).rjust(12)
    output = output + ("%.3f" % numpy_array[1]).rjust(8)
    output = output + ("%.3f" % numpy_array[2]).rjust(8)
    output = output + letter.rjust(24)
    
    return output

def numpy_to_pdb(narray, letter, resname=""):
    """Create a string formatted according to the PDB standard.

    Arguments:
    narray -- A nx3 numpy.array representing a 3D point.
    letter -- A string, the atom name/chain/etc to use for the output.

    Returns:
    (Optionally) A string, formatted according to the PDB standard.

    """
    
    if len(narray.flatten()) == 3:
        return create_pdb_line(narray, 1, "AAA", letter) + "\n"
    else:
        if resname == "":
            letters = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "W", "X", "Y", "Z"]
            resnames = []
            for l1 in letters:
                for l2 in letters:
                    for l3 in letters:
                        resnames.append(l1+l2+l3)
            resnames.remove("XXX") # because this is reserved for empty atoms
        else:
            resnames = [resname]

        t = ""
        for i, item in enumerate(narray): t = t + create_pdb_line(item, i+1, resnames[i % len(resnames)], letter) + "\n"
        return t

def dx_freq(freq_mat, parameters):
    '''
    Generates a DX file that records the frequency that a volume element is open
    
    Arguments:
    freq_mat -- a Nx4 matrix, where the first 3 columns are the x,y,z coords of the point, and the 4th column is the frequency of emptiness for that point in space

    '''
    
    header_template = """# Data from POVME 2.0.1
# 
# FREQUENCY (unitless)
# 
object 1 class gridpositions counts %d %d %d
origin %8.6e %8.6e %8.6e
delta %8.6e 0.000000e+00 0.000000e+00
delta 0.000000e+00 %8.6e 0.000000e+00
delta 0.000000e+00 0.000000e+00 %8.6e
object 2 class gridconnections counts %d %d %d
object 3 class array type double rank 0 items %d data follows
"""
    
    footer_template = """
attribute "dep" string "positions"
object "regular positions regular connections" class field
component "positions" value 1
component "connections" value 2
component "data" value 3"""
    
    # 1. Sort the points into the proper order for a dx file
    
    # already sorted correctly
    
    # 2. Obtain key information about the grid
    N = freq_mat.shape[0] # number of data points
    
    minx = min(freq_mat[:,0])
    miny = min(freq_mat[:,1])
    minz = min(freq_mat[:,2]) # find the upper and lower corners of the grid
    maxx = max(freq_mat[:,0])
    maxy = max(freq_mat[:,1])
    maxz = max(freq_mat[:,2])
    
    widthx = maxx - minx # find the widths of the grid
    widthy = maxy - miny
    widthz = maxz - minz
    
    xs = numpy.unique(freq_mat[:,0])
    ys = numpy.unique(freq_mat[:,1])
    zs = numpy.unique(freq_mat[:,2])
    
    resx = xs[1]- xs[0]
    resy = ys[1]- ys[0]
    resz = zs[1]- zs[0]
    
    #resx = freq_mat[(widthz+1)*(widthy+1),0] - freq_mat[0,0]
    #resy = freq_mat[widthz+1,1] - freq_mat[0,1] # find the resolution of the grid
    #resz = freq_mat[1,2] - freq_mat[0,2]

    nx = (widthx) / resx + 1 # number of grid points in each dimension
    ny = (widthy) / resy + 1 # need to add one because the subtraction leaves out an entire row
    nz = (widthz) / resz + 1

    # test to make sure all is well with the size of the grid and its dimensions
    assert (nx * ny * nz) == N, "Something is wrong with the freq_mat array: it is not a prismatic shape"
    
    # 3. write the header and footer
    if parameters['SaveVolumetricDensityMap'] == True:
        if parameters['CompressOutput'] == True: dx_file = gzip.open(parameters['OutputFilenamePrefix'] + "volumetric_density.dx.gz",'wb')
        else: dx_file = open(parameters['OutputFilenamePrefix'] + "volumetric_density.dx",'w')
        
        header = header_template % (nx, ny, nz, minx, miny, minz, resx, resy, resz, nx, ny, nz, N) # format the header
        footer = footer_template # the footer needs no formatting
        dx_file.write(header)
        newline_counter = 1
        for i in range(N): # write the data to the DX file
            dx_file.write("%8.6e" % freq_mat[i,3])
            if newline_counter == 3:
                newline_counter = 0
                dx_file.write("\n")
            else:
                dx_file.write(" ")
            newline_counter += 1
        dx_file.write(footer)
        dx_file.close
    return

class MultithreadingCalcVolumeTask(MultithreadingTaskGeneral):
    '''A class for calculating the volume.'''
    
    def value_func(self, item, results_queue):
        """Calculate the volume.

        Arguments:
        item -- A list or tuple, the input data required for the calculation
        results_queue -- A multiprocessing.Queue() object for storing the calculation output
        
        """

        frame_indx = item[0]
        pdb = item[1]
        parameters = item[2]
        # pts = parameters['pts_orig'].copy() # this works
        pts = parameters['pts_orig'] # also works, so keep because faster
        
        # if the user wants to save empty points (points that are removed), then we need a copy of the original
        if parameters['OutputEqualNumPointsPerFrame'] == True:
            pts_deleted = pts.copy()
        
        # you may need to load it from disk if the user so specified
        if parameters['UseDiskNotMemory'] == True: # so you need to load it from disk
            pym_filename = pdb
            pdb = pymolecule.Molecule()
            pdb.fileio.load_pym_into(pym_filename)
        
        # remove the points that are far from the points region anyway
        min_pts = numpy.min(pts,0) - parameters['DistanceCutoff'] - 1
        max_pts = numpy.max(pts,0) + parameters['DistanceCutoff'] + 1
        
        # identify atoms that are so far away from points that they can be ignored
        index_to_keep1 = numpy.nonzero((pdb.information.coordinates[:,0] > min_pts[0]))[0] # x's too small
        index_to_keep2 = numpy.nonzero((pdb.information.coordinates[:,0] < max_pts[0]))[0] # x's too large
        index_to_keep3 = numpy.nonzero((pdb.information.coordinates[:,1] > min_pts[1]))[0] # y's too small
        index_to_keep4 = numpy.nonzero((pdb.information.coordinates[:,1] < max_pts[1]))[0] # y's too large
        index_to_keep5 = numpy.nonzero((pdb.information.coordinates[:,2] > min_pts[2]))[0] # z's too small
        index_to_keep6 = numpy.nonzero((pdb.information.coordinates[:,2] < max_pts[2]))[0] # z's too large
        
        index_to_keep = numpy.intersect1d(index_to_keep1, index_to_keep2, assume_unique=True)
        index_to_keep = numpy.intersect1d(index_to_keep, index_to_keep3, assume_unique=True)
        index_to_keep = numpy.intersect1d(index_to_keep, index_to_keep4, assume_unique=True)
        index_to_keep = numpy.intersect1d(index_to_keep, index_to_keep5, assume_unique=True)
        index_to_keep = numpy.intersect1d(index_to_keep, index_to_keep6, assume_unique=True)
        
        # keep only relevant atoms
        if len(index_to_keep) > 0: pdb = pdb.selections.create_molecule_from_selection(index_to_keep)
    
        # get the vdw radii of each protein atom
        vdw = numpy.ones(len(pdb.information.coordinates)) # so the default vdw is 1.0
        
        # get vdw... you might want to fill this out with additional vdw values
        vdw[numpy.nonzero(pdb.information.atom_information['element_stripped'] == "H")[0]] = 1.2
        vdw[numpy.nonzero(pdb.information.atom_information['element_stripped'] == "C")[0]] = 1.7
        vdw[numpy.nonzero(pdb.information.atom_information['element_stripped'] == "N")[0]] = 1.55
        vdw[numpy.nonzero(pdb.information.atom_information['element_stripped'] == "O")[0]] = 1.52
        vdw[numpy.nonzero(pdb.information.atom_information['element_stripped'] == "F")[0]] = 1.47
        vdw[numpy.nonzero(pdb.information.atom_information['element_stripped'] == "P")[0]] = 1.8
        vdw[numpy.nonzero(pdb.information.atom_information['element_stripped'] == "S")[0]] = 1.8
        vdw = numpy.repeat(numpy.array([vdw]).T, len(pts), axis=1)
    
        # now identify the points that are close to the protein atoms
        dists = cdist(pdb.information.coordinates, pts)
        close_pt_index = numpy.nonzero((dists < (vdw + parameters['DistanceCutoff'])))[1]
        
        # now keep the appropriate points
        pts = numpy.delete(pts, close_pt_index, axis=0)
    
        # exclude points outside convex hull
        if parameters['ConvexHullExclusion'] == True:
            convex_hull_3d = ConvexHull()
            
            # get the coordinates of the non-hydrogen atoms (faster to discard hydrogens)
            hydros = pdb.selections.select_atoms({'element_stripped':['H']})
            not_hydros = pdb.selections.invert_selection(hydros)
            not_hydros_coors = pdb.information.coordinates[not_hydros]

            #not_hydros = pdb.selections.select_atoms({'name_stripped':['CA']})
            #not_hydros_coors = pdb.information.coordinates[not_hydros]

            # modify pts here.
            # note that the atoms of the pdb frame are in pdb.information.coordinates
            #begintime = time.time() # measure execution time
            akl_toussaint_pts = convex_hull_3d.akl_toussaint(not_hydros_coors) # quickly reduces input size
            #print "akl Toussaint:", time.time() - begintime
            begintime = time.time() # measure execution time
            hull = convex_hull_3d.gift_wrapping_3d(akl_toussaint_pts) # calculate convex hull using gift wrapping algorithm
            #print "gift_wrapping:", time.time() - begintime
            
            old_pts = pts # we will need to regenerate the pts list, disregarding those outside the hull
            pts = []
            for pt in old_pts:
              pt_outside = convex_hull_3d.outside_hull(pt, hull) # check if pt is outside hull
              if not pt_outside:
                pts.append(pt) # if its not outside the hull, then include it in the volume measurement
            pts = numpy.array(pts)
            
        # Now, enforce contiguity if needed
        if len(parameters['ContiguousPocketSeedRegions']) > 0 and len(pts) > 0:
            # first, for each point, determine how many neighbors it has
            cutoff_dist = parameters['GridSpacing'] * 1.01 * math.sqrt(3) # to count kiddy-corner points too
            pts_dists = squareform(pdist(pts))
            neighbor_counts = numpy.sum(pts_dists < cutoff_dist,axis=0) - 1 # minus 1 because an atom shouldn't be considered its own neighor
            
            # remove all the points that don't have enough neighbors
            pts = pts[numpy.nonzero(neighbor_counts >= parameters['ContiguousPointsCriteria'])[0]]
            
            # get all the points in the defined parameters['ContiguousPocket'] seed regions
            contig_pts = parameters['ContiguousPocketSeedRegions'][0].points_set(parameters['GridSpacing'])
            for Contig in parameters['ContiguousPocketSeedRegions'][1:]: contig_pts = numpy.vstack((contig_pts, Contig.points_set(parameters['GridSpacing'])))
            contig_pts = unique_rows(contig_pts)
            
            try: # error here if there are no points of contiguous seed region outside of protein volume.
                # now just get the ones that are not near the protein
                contig_pts = pts[numpy.nonzero(cdist(contig_pts, pts) < 1e-7)[1]]
                
                last_size_of_contig_pts = 0
                while last_size_of_contig_pts != len(contig_pts):
                    last_size_of_contig_pts = len(contig_pts)
                    
                    # now get the indecies of all points that are close to the contig_pts
                    all_pts_close_to_contig_pts_boolean = (cdist(pts, contig_pts) < cutoff_dist)
                    index_all_pts_close_to_contig_pts = numpy.unique(numpy.nonzero(all_pts_close_to_contig_pts_boolean)[0])
                    contig_pts = pts[index_all_pts_close_to_contig_pts]
                    
                pts = contig_pts
            except:
                log("\tFrame " + str(frame_indx) + ": None of the points in the contiguous-pocket seed region\n\t\tare outside the volume of the protein! Assuming a pocket\n\t\tvolume of 0.0 A.", parameters)
                pts = numpy.array([])

        # now write the pdb and calculate the volume
        volume = len(pts) * math.pow(parameters['GridSpacing'],3)

        log("\tFrame " + str(frame_indx) + ": " + repr(volume) + " A^3", parameters)
        if parameters['SaveIndividualPocketVolumes'] == True:
            frame_text = ""
            frame_text = frame_text + "REMARK Frame " + str(frame_indx) + "\n"
            frame_text = frame_text + "REMARK Volume = " + repr(volume) + " Cubic Angtroms\n"
            frame_text = frame_text + numpy_to_pdb(pts,'X')

            if parameters['OutputEqualNumPointsPerFrame'] == True:
                # you need to find the points that are in pts_deleted but not in pts
                tmp = reduce(lambda x, y: x |  numpy.all(pts_deleted == y, axis=-1), pts, numpy.zeros(pts_deleted.shape[:1], dtype=numpy.bool))
                indices = numpy.where(tmp)[0]
                pts_deleted = numpy.delete(pts_deleted, indices, axis=0)

                pts_deleted = numpy.zeros(pts_deleted.shape) # So extra points will always be at the origin. These can be easily hidden with your visualization software.
                frame_text = frame_text + numpy_to_pdb(pts_deleted,'X',"XXX")
                
            frame_text = frame_text + "END\n"

            if parameters['CompressOutput'] == True: fl = gzip.open(parameters['OutputFilenamePrefix'] + "frame_" + str(frame_indx) + ".pdb.gz", 'wb')
            else: fl = open(parameters['OutputFilenamePrefix'] + "frame_" + str(frame_indx) + ".pdb", 'w')
            fl.write(frame_text)
            fl.close()

        extra_data_to_add = {}
        if parameters['SaveVolumetricDensityMap'] == True: extra_data_to_add['SaveVolumetricDensityMap'] = pts

        self.results.append((frame_indx, volume, extra_data_to_add))
        
        #if len(extra_data_to_add.keys()) != 0: 
        #else: self.results.append((frame_indx, volume))

class MultithreadingStringToMoleculeTask(MultithreadingTaskGeneral):
    '''A class for loading PDB frames (as strings) into pymolecule.Molecule objects.'''
    
    def value_func(self, item, results_queue):
        """Convert a PDB string into a pymolecule.Molecule object

        Arguments:
        item -- A list or tuple, the input data required for the calculation
        results_queue -- A multiprocessing.Queue() object for storing the calculation output
        
        """
        
        pdb_string = item[0]
        index = item[1]
        parameters = item[2]
        
        # make the pdb object
        str_obj = StringIO(pdb_string)
        tmp = pymolecule.Molecule()
        tmp.fileio.load_pdb_into_using_file_object(str_obj, False, False, False)

        log("\tFurther processing frame " + str(index), parameters)
        
        if parameters['UseDiskNotMemory'] == False: # so load the whole trajectory into memory
            self.results.append((index, tmp))
        else: # save to disk, record filename
            pym_filename = "." + os.sep + ".povme_tmp" + os.sep + "frame_" + str(index) + ".pym"
            tmp.fileio.save_pym(pym_filename, False, False, False, False, False)
            self.results.append((index, pym_filename))

class Region:
    '''A class for defining regions that will be filled with points.'''
    
    def __init__(self):
        '''Initialize some variables.'''
        
        self.center = numpy.array([9999.9, 9999.9, 9999.9])
        self.radius = 9999.9 # in case the region is a sphere
        self.box_dimen = numpy.array([9999.9, 9999.9, 9999.9]) # in case the region is a box
        
        self.region_type = "SPHERE" # could also be BOX

    def __str__(self):
        '''Returns a string representation of the region.'''
        
        if self.region_type == "SPHERE": return "sphere at (" + str(self.center[0]) + ", " + str(self.center[1]) + ", " + str(self.center[2]) + "), radius = " + str(self.radius)
        if self.region_type == "BOX": return "box centered at (" + str(self.center[0]) + ", " + str(self.center[1]) + ", " + str(self.center[2]) + ") with x,y,z dimensions of (" + str(self.box_dimen[0]) + ", " + str(self.box_dimen[1]) + ", " + str(self.box_dimen[2]) + ")"
        return ''

    def __snap(self, pts, reso):
        """Snaps a set of points to a fixed grid.
    
        Arguments:
        pts -- A nx3 numpy.array representing 3D points.
        reso -- A float, the resolution of the grid.
    
        Returns:
        A nx3 numpy.array with the 3D points snapped to the nearest grid point.
    
        """
        
        # unfortunately, numpy.around rounds evenly, so 0.5 rounds to 0.0 and 1.5 rounds to 2.0.
        # very annoying, I'll just add a tiny amount to 0.5 => 0.500001
        # this should work, since user is unlikely to select region center or radius with such
        # precision
        
        pts = pts + 1e-10
        return numpy.around(pts/reso) * reso

    def points_set(self, reso):
        """Generates a point field by filling the region with equally spaced points.
    
        Arguments:
        reso -- A float, the resolution of the grid on which the points will be placed.
    
        Returns:
        A nx3 numpy.array with the 3D points filling the region.
    
        """
        
        total_pts = None
        
        if self.region_type == "BOX":
            xs = numpy.arange(self.center[0] - self.box_dimen[0] / 2, self.center[0] + self.box_dimen[0] / 2, reso)
            ys = numpy.arange(self.center[1] - self.box_dimen[1] / 2, self.center[1] + self.box_dimen[1] / 2, reso)
            zs = numpy.arange(self.center[2] - self.box_dimen[2] / 2, self.center[2] + self.box_dimen[2] / 2, reso)
            
            total_pts = numpy.empty((len(xs) * len(ys) * len(zs), 3))

            i = 0
            for x in xs:
                for y in ys:
                    for z in zs:
                        total_pts[i][0] = x
                        total_pts[i][1] = y
                        total_pts[i][2] = z
                        
                        i = i + 1

            total_pts = self.__snap(total_pts, reso)

        elif self.region_type == "SPHERE":
            xs = numpy.arange(self.center[0] - self.radius, self.center[0] + self.radius, reso)
            ys = numpy.arange(self.center[1] - self.radius, self.center[1] + self.radius, reso)
            zs = numpy.arange(self.center[2] - self.radius, self.center[2] + self.radius, reso)
            
            total_pts = numpy.empty((len(xs) * len(ys) * len(zs), 3))

            i = 0
            for x in xs:
                for y in ys:
                    for z in zs:
                        total_pts[i][0] = x
                        total_pts[i][1] = y
                        total_pts[i][2] = z
                        
                        i = i + 1

            total_pts = self.__snap(total_pts, reso)
            
            # now remove all the points outside of this sphere
            index_inside_sphere = numpy.nonzero(cdist(total_pts, numpy.array([self.center])) < self.radius)[0]
            total_pts = total_pts[index_inside_sphere]
                        
        return total_pts

class ConfigFile:
    '''A class for processing the user-provided configuration file.'''
    
    entities = []
    
    def __init__ (self, FileName):
        """Generates a point field by filling the region with equally spaced points.
    
        Arguments:
        FileName -- A string, the filename of the configuration file.
    
        """

        f = open(FileName,'r')
        lines = f.readlines()
        f.close()
        
        for line in lines:
            # remove comments
            line = line.split("#",1)[0]
            # line = line.split("//",1)[0] # We can't have these kinds of comments any more because of Windows filenames.

            line = line.strip()

            if line != "":

                # replace ; and , and : with space
                # line = line.replace(',',' ')
                # line = line.replace(';',' ')
                # line = line.replace(':',' ') # this messes up Windows filenames
                line = line.replace("\t",' ')
                
                # now strip string
                line = line.strip()
    
                # now, replace double spaces with one space
                while '  ' in line: line = line.replace('  ',' ')
                
                # Now split the thing
                line = line.split(' ',1)
    
                # now, make it upper case
                line[0] = line[0].upper()
                
                # If there's QUIT, EXIT, or STOP, then don't continue.
                if line[0] in ['QUIT','EXIT','STOP']: break
                
                self.entities.append(line)
    
class runit():
    '''The main class to run POVME.'''

    def reference(self, parameters, before=""):
        '''Print out a message regarding terms of use.'''
        
        log("", parameters)
        log(before + "If you use POVME in your research, please cite the following reference:", parameters)
        log(before + "  Durrant, J. D., C. A. de Oliveira, et al. (2011). \"POVME: An algorithm", parameters)
        log(before + "  for measuring binding-pocket volumes.\" J Mol Graph Model 29(5): 773-776.", parameters)

    def load_multi_frame_pdb(self, filename, parameters):
        """Load a multi-frame PDB into memory or into separate files (depending on user specifications).
    
        Arguments:
        filename -- A string, the filename of the multi-frame PDB
        parameters -- A python dictionary, where the keys are the user-defined parameter names and the values are the corresponding parameter values.
    
        Returns:
        If the user has requested that the disk be used to save memory, this function returns a list of tuples, where the first item in each tuple is the frame index, and the second is a filename containing the individual frame.
        If memory is to be used instead of the disk, this function returns a list of tuples, where the first item in each tuple is the frame index, and the second is a pymolecule.Molecule object representing the frame.
    
        """
        
        pdb_strings = []
        growing_string = ''
        
        log("", parameters)
        log("Reading frames from " + filename, parameters)
        
        f = open(filename, 'rb')
        while True:
            
            if parameters['NumFrames'] != -1:
                if len(pdb_strings) >= parameters['NumFrames']: break
            
            line = f.readline()
            
            if len(line) == 0:
                pdb_strings.append(growing_string)
                break
            if line[:3] == "END":
                pdb_strings.append(growing_string)
                growing_string = ''
            else:
                growing_string = growing_string + line
            
        f.close()
        
        while '' in pdb_strings: pdb_strings.remove('')
        
        # now convert each pdb string into a pymolecule.Molecule object
        molecules = Multithreading([(pdb_strings[idx], idx + 1, parameters) for idx in range(len(pdb_strings))], parameters['NumProcessors'], MultithreadingStringToMoleculeTask)
        molecules = molecules.results

        return molecules
    
    def __init__(self, argv):
        '''Start POVME

        Arguments:
        argv -- A list of the command-line arguments.
        
        '''
        
        start_time = time.time()

        # Load the configuration file
        if len(argv) == 1:
                print "\nPOVME " + version
                print "\nPlease specify the input file from the command line!\n\nExample: python POVME.py input_file.ini"
                self.reference({})
                print
                sys.exit()

        config = ConfigFile(argv[1])
        
        # Process the config file
        parameters = {}
        
        parameters['GridSpacing'] = 1.0 # default
        parameters['PointsIncludeRegions'] = []
        parameters['PointsExcludeRegions'] = []
        parameters['SavePoints'] = False # default
        parameters['LoadPointsFilename'] = '' # default
        parameters['PDBFileName'] = "" # default
        parameters['DistanceCutoff'] = 1.09 # default is VDW radius of hydrogen
        parameters['ConvexHullExclusion'] = True
        parameters['ContiguousPocketSeedRegions'] = []
        parameters['ContiguousPointsCriteria'] = 4
        parameters['NumProcessors'] = 4
        parameters['UseDiskNotMemory'] = False
        parameters['OutputFilenamePrefix'] = "POVME_output." + time.strftime("%m-%d-%y") + "." + time.strftime("%H-%M-%S") + os.sep
        parameters['SaveIndividualPocketVolumes'] = False
        parameters['SavePocketVolumesTrajectory'] = False
        parameters['OutputEqualNumPointsPerFrame'] = False
        parameters['SaveTabbedVolumeFile'] = False
        parameters['SaveVolumetricDensityMap'] = False
        parameters['CompressOutput'] = False
        parameters['NumFrames'] = -1 # This is a parameter for debugging purposes only.
        
        float_parameters = ["GridSpacing", "DistanceCutoff"]
        boolean_parameters = ["SavePoints", "ConvexHullExclusion", "CompressOutput", "UseDiskNotMemory", "SaveVolumetricDensityMap", "OutputEqualNumPointsPerFrame", "SaveIndividualPocketVolumes", "SaveTabbedVolumeFile", "SavePocketVolumesTrajectory"]
        int_parameters = ["NumFrames", "ContiguousPointsCriteria", "NumProcessors"]
        string_parameters = ["OutputFilenamePrefix", "PDBFileName", "LoadPointsFilename"]
        
        print config.entities
        
        for entity in config.entities:
            try:
                index = [p.upper() for p in float_parameters].index(entity[0])
                parameters[float_parameters[index]] = float(entity[1])
            except: pass
            
            try: 
                index = [p.upper() for p in boolean_parameters].index(entity[0])
                if entity[1].upper() in ["YES", "TRUE"]: parameters[boolean_parameters[index]] = True
                else: parameters[boolean_parameters[index]] = False
            except: pass

            try:
                index = [p.upper() for p in int_parameters].index(entity[0])
                parameters[int_parameters[index]] = int(entity[1])
            except: pass

            try:
                index = [p.upper() for p in string_parameters].index(entity[0])
                parameters[string_parameters[index]] = entity[1].strip()
            except: pass

            # Regions are handled separately for each parameter...
            if entity[0] == "POINTSINCLUSIONSPHERE":
                Include = Region()
                items = entity[1].split(' ')
                Include.center[0] = float(items[0])
                Include.center[1] = float(items[1])
                Include.center[2] = float(items[2])
                Include.radius = float(items[3])
                Include.region_type = "SPHERE"
                parameters['PointsIncludeRegions'].append(Include)
            elif entity[0] == "POINTSINCLUSIONBOX":
                Include = Region()
                items = entity[1].split(' ')
                Include.center[0] = float(items[0])
                Include.center[1] = float(items[1])
                Include.center[2] = float(items[2])
                Include.box_dimen[0] = float(items[3])
                Include.box_dimen[1] = float(items[4])
                Include.box_dimen[2] = float(items[5])
                Include.region_type = "BOX"
                parameters['PointsIncludeRegions'].append(Include)
            if entity[0] == "CONTIGUOUSPOCKETSEEDSPHERE":
                Contig = Region()
                items = entity[1].split(' ')
                Contig.center[0] = float(items[0])
                Contig.center[1] = float(items[1])
                Contig.center[2] = float(items[2])
                Contig.radius = float(items[3])
                Contig.region_type = "SPHERE"
                parameters['ContiguousPocketSeedRegions'].append(Contig)
            elif entity[0] == "CONTIGUOUSPOCKETSEEDBOX":
                Contig = Region()
                items = entity[1].split(' ')
                Contig.center[0] = float(items[0])
                Contig.center[1] = float(items[1])
                Contig.center[2] = float(items[2])
                Contig.box_dimen[0] = float(items[3])
                Contig.box_dimen[1] = float(items[4])
                Contig.box_dimen[2] = float(items[5])
                Contig.region_type = "BOX"
                parameters['ContiguousPocketSeedRegions'].append(Contig)
            elif entity[0] == "POINTSEXCLUSIONSPHERE":
                Exclude = Region()
                items = entity[1].split(' ')
                Exclude.center[0] = float(items[0])
                Exclude.center[1] = float(items[1])
                Exclude.center[2] = float(items[2])
                Exclude.radius = float(items[3])
                Exclude.region_type = "SPHERE"
                parameters['PointsExcludeRegions'].append(Exclude)
            elif entity[0] == "POINTSEXCLUSIONBOX":
                Exclude = Region()
                items = entity[1].split(' ')
                Exclude.center[0] = float(items[0])
                Exclude.center[1] = float(items[1])
                Exclude.center[2] = float(items[2])
                Exclude.box_dimen[0] = float(items[3])
                Exclude.box_dimen[1] = float(items[4])
                Exclude.box_dimen[2] = float(items[5])
                Exclude.region_type = "BOX"
                parameters['PointsExcludeRegions'].append(Exclude)
        
        # If the output prefix includes a directory, create that directory if necessary
        if os.sep in parameters['OutputFilenamePrefix']:
            output_dirname = os.path.dirname(parameters['OutputFilenamePrefix'])
            #if os.path.exists(output_dirname): shutil.rmtree(output_dirname) # So delete the directory if it already exists.
            try: os.mkdir(output_dirname)
            except: pass

        # print out the header
        self.reference(parameters, "")
        log('', parameters)
        
        # create temp swap directory if needed
        if parameters['UseDiskNotMemory'] == True:
            if os.path.exists('.' + os.sep + '.povme_tmp'): shutil.rmtree('.' + os.sep + '.povme_tmp')
            os.mkdir('.' + os.sep + '.povme_tmp')

        # print out parameters
        log("Parameters:", parameters)
        for i in parameters.keys():
            
            if i == 'NumFrames' and parameters['NumFrames'] == -1: continue # So only show this parameter if it's value is not the default.
            
            if type(parameters[i]) is list:
                for i2 in parameters[i]:
                    if i2 != "": log("\t" + str(i) + ": " + str(i2), parameters)
            else:
                if parameters[i] != "": log("\t" + str(i) + ": " + str(parameters[i]), parameters)

        pts = None
        if len(parameters['PointsIncludeRegions']) > 0: # so create the point file
        
            log("\nGenerating the pocket-encompassing point field", parameters)
        
            # get all the points of the inclusion regions
            pts = parameters['PointsIncludeRegions'][0].points_set(parameters['GridSpacing'])
            for Included in parameters['PointsIncludeRegions'][1:]: pts = numpy.vstack((pts, Included.points_set(parameters['GridSpacing'])))
            pts = unique_rows(pts)
            
            # get all the points of the exclusion regions
            if len(parameters['PointsExcludeRegions']) > 0:
                pts_exclusion = parameters['PointsExcludeRegions'][0].points_set(parameters['GridSpacing'])
                for Excluded in parameters['PointsExcludeRegions'][1:]: pts_exclusion = numpy.vstack((pts_exclusion, Excluded.points_set(parameters['GridSpacing'])))
                pts_exclusion = unique_rows(pts_exclusion)
            
                # remove the exclusion points from the inclusion points
                # I think there ought to be a set-based way of doing this,
                # but I'm going to go for the pairwise comparison.
                # consider rewriting later
                index_to_remove = numpy.nonzero(cdist(pts, pts_exclusion) < 1e-7)[0]
                pts = numpy.delete(pts, index_to_remove, axis=0)
                
            # save the points as PDB
            if parameters['SavePoints'] == True:
                
                # First, save the point field itself
                
                log("\nSaving the point field as a PDB and NPY file", parameters)
                
                points_filename = parameters['OutputFilenamePrefix'] + "point_field.pdb"
                
                if parameters['CompressOutput'] == True: afile = gzip.open(points_filename + ".gz", 'wb')
                else: afile = open(points_filename,'w')
                
                afile.write(numpy_to_pdb(pts, "X"))
                afile.close()
                
                # save the points as npy
                numpy.save(points_filename + ".npy", pts)
            
                log("\tPoint field saved to " + points_filename + " to permit visualization", parameters)
                log("\tPoint field saved to " + points_filename + ".npy to optionally load for the volume calculation", parameters)
                log("", parameters)
                
                # Now, save the contiguous seed points as well, if specified.
                if len(parameters['ContiguousPocketSeedRegions']) > 0:
                    # get all the contiguous points
                    contig_pts = parameters['ContiguousPocketSeedRegions'][0].points_set(parameters['GridSpacing'])
                    for Contig in parameters['ContiguousPocketSeedRegions'][1:]: contig_pts = numpy.vstack((contig_pts, Contig.points_set(parameters['GridSpacing'])))
                    contig_pts = unique_rows(contig_pts)
                    
                    log("\nSaving the contiguous-pocket seed points as a PDB file", parameters)
                    
                    points_filename = parameters['OutputFilenamePrefix'] + "contiguous_pocket_seed_points.pdb"
                    
                    if parameters['CompressOutput'] == True: afile = gzip.open(points_filename + ".gz", 'wb')
                    else: afile = open(points_filename,'w')
                    
                    afile.write(numpy_to_pdb(contig_pts, "X"))
                    afile.close()
                    
                    log("\tContiguous-pocket seed points saved to " + points_filename + " to permit visualization", parameters)
                    log("", parameters)
                
        if parameters['PDBFileName'] != '': # so there's a PDB point specified for calculating the volume.
            
            # load the points in they aren't already present
            if pts is None:
                log("\nLoading the point-field NPY file...", parameters)
                parameters['pts_orig'] = numpy.load(parameters['LoadPointsFilename'])
            else: parameters['pts_orig'] = pts

            # load the PDB frames
            index_and_pdbs = self.load_multi_frame_pdb(parameters['PDBFileName'], parameters)
            
            # calculate all the volumes
            log("", parameters)
            log("Calculating the pocket volume of each frame", parameters)
            tmp = Multithreading([(index, pdb_object, parameters) for index, pdb_object in index_and_pdbs], parameters['NumProcessors'], MultithreadingCalcVolumeTask)

            # delete the temp swap directory if necessary
            if parameters['UseDiskNotMemory'] == True:
                if os.path.exists('.' + os.sep + '.povme_tmp'): shutil.rmtree('.' + os.sep + '.povme_tmp')

            # display the results
            results_dic = {}
            for result in tmp.results: results_dic[result[0]] = result[1]
            log("", parameters)
            log("FRAME        | VOLUME (A^3)", parameters)
            log("-------------+-------------", parameters)
            for i in sorted(results_dic.keys()): log(str(i).ljust(13) + "| " + str(results_dic[i]), parameters)
            
            log("", parameters)
            log("Execution time = " + str(time.time()-start_time) + " sec", parameters)
            log("", parameters)
            
            # if the user requested a separate volume file, save that as well
            if parameters['SaveTabbedVolumeFile'] == True:
                if parameters['CompressOutput'] == True: f = gzip.open(parameters['OutputFilenamePrefix'] + "volumes.tabbed.txt.gz", 'wb')
                else: f = open(parameters['OutputFilenamePrefix'] + "volumes.tabbed.txt", 'w')
                
                for i in sorted(results_dic.keys()): f.write(str(i) + "\t" + str(results_dic[i]) + "\n")
                f.close()
            
            # if the user wanted a single trajectory containing all the volumes, generate that here.
            if parameters['SavePocketVolumesTrajectory'] == True:
                if parameters['CompressOutput'] == True: traj_file = gzip.open(parameters['OutputFilenamePrefix'] + "volume_trajectory.pdb.gz", 'wb')
                else: traj_file = open(parameters['OutputFilenamePrefix'] + "volume_trajectory.pdb", 'w')

                for frame_index in range(1,len(results_dic.keys())+1):
                    if parameters['CompressOutput'] == True: frame_file = gzip.open(parameters['OutputFilenamePrefix'] + "frame_" + str(frame_index) + ".pdb.gz", 'rb')
                    else: frame_file = open(parameters['OutputFilenamePrefix'] + "frame_" + str(frame_index) + ".pdb", 'r')
                    
                    traj_file.write(frame_file.read())
                    frame_file.close()

                traj_file.close()            

            # if the user requested a volumetric density map, then generate it here
            if parameters['SaveVolumetricDensityMap'] == True:
                unique_points = {}
                
                overall_min = numpy.ones(3) * 1e100
                overall_max = numpy.ones(3) * -1e100

                for result in tmp.results:
                    pts = result[2]['SaveVolumetricDensityMap']

                    if len(pts) > 0:
                        amin = numpy.min(pts,axis=0)
                        amax = numpy.max(pts,axis=0)
                        
                        overall_min = numpy.min(numpy.vstack((overall_min, amin)), axis=0)
                        overall_max = numpy.max(numpy.vstack((overall_max, amax)), axis=0)
                        
                        for pt in pts:
                            pt_key = str(pt[0]) + ";" + str(pt[1]) + ";" + str(pt[2])
                            try: unique_points[pt_key] = unique_points[pt_key] + 1
                            except: unique_points[pt_key] = 1
                if overall_min[0] == 1e100:
                    log("ERROR! Cannont save volumetric density file because no volumes present in any frame.", parameters)
                else:
                    xpts = numpy.arange(overall_min[0], overall_max[0] + parameters['GridSpacing'], parameters['GridSpacing'])
                    ypts = numpy.arange(overall_min[1], overall_max[1] + parameters['GridSpacing'], parameters['GridSpacing'])
                    zpts = numpy.arange(overall_min[2], overall_max[2] + parameters['GridSpacing'], parameters['GridSpacing'])
                    
                    all_pts = numpy.zeros((len(xpts)*len(ypts)*len(zpts), 4))
                    
                    i = 0
                    for x in xpts:
                        for y in ypts:
                            for z in zpts:
                                key = str(x) + ";" + str(y) + ";" + str(z)
                                all_pts[i][0] = x
                                all_pts[i][1] = y
                                all_pts[i][2] = z
                                
                                try: all_pts[i][3] = unique_points[key]
                                except: pass
                                
                                i = i + 1
                    
                    # convert the counts in the fourth column into frequencies
                    all_pts[:,3] = all_pts[:,3] / len(tmp.results)
                    dx_freq(all_pts, parameters) # save the dx file
                    
                    #print "To turn into a DX file:"
                    #print all_pts
                    #import cPickle as pickle
                    #pickle.dump(all_pts, open('dill.pickle', 'w'))
                
if __name__ == "__main__": dorun = runit(sys.argv)
