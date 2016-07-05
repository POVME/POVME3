# POVME Pocket ID 1.0 is released under the GNU General Public License 
# (see http://www.gnu.org/licenses/gpl.html).

# If you have any questions, comments, or suggestions, please don't hesitate to contact me,
# Jacob Durrant, at jdurrant [at] ucsd [dot] edu.

import sys
import numpy
from scipy import spatial
from scipy.cluster.vq import kmeans2
from scipy.spatial.distance import cdist
import textwrap
import getopt
from numpy.lib.recfunctions import append_fields
import multiprocessing

# POVME Pocket ID 1.0 is a program for identifying protein pockets and generating
# appropriate pocket-encompassing inclusion spheres. These spheres, modified as required,
# can then be used as POVME input.

# Some classes are required to support the loading and manipulation of 3D molecular information

class Information():
    """A class for storing and accessing information about the elements of a Molecule object"""
    
    def __init__(self, parent_molecule_object):
        """Initializes the Information class.
                
            Arguments:
            parent_molecule_object -- The Molecule object associated with this class.
            
            """
        
        self.__parent_molecule = parent_molecule_object

        self.__constants = {}
        self.__constants['i8_fields'] = ['serial','resseq']
        self.__constants['f8_fields']= ['x','y','z','occupancy','tempfactor']

        self.__atom_information = None
        self.__coordinates = None

    def get_atom_information(self): return self.__atom_information
    def get_coordinates(self): return self.__coordinates
    def get_constants(self): return self.__constants
    
    def set_atom_information(self,atom_information): self.__atom_information = atom_information
    def set_coordinates(self,coordinates): self.__coordinates = coordinates

    def get_bounding_box(self, selection = None, padding=0.0):
        """Calculates a box that bounds (encompasses) a set of atoms.
            
            Arguments:
            selection -- An optional numpy.array containing the indices of the atoms to consider. If ommitted, all atoms of the Molecule object will be considered.
            padding -- An optional float. The bounding box will extend this many angstroms beyond the atoms being considered.
            
            Returns:
            A numpy array representing two 3D points, (min_x, min_y, min_z) and (max_x, max_y, max_z), that bound the molecule.
            
            """
        
        if selection is None: selection = self.__parent_molecule.select_all()
        
        return numpy.vstack((numpy.min(self.__coordinates[selection],0), numpy.max(self.__coordinates[selection],0)))
    
class FileIO():
    """A class for saving and loading molecular data into a Molecule object"""
    
    def __init__(self, parent_molecule_object):
        """Initializes the FileIO class.
                
            Arguments:
            parent_molecule_object -- The Molecule object associated with this class.
            
            """
        
        self.__parent_molecule = parent_molecule_object
        
    def load_pdb_into(self, filename):
        """Loads the molecular data contained in a pdb file into the current Molecule object.
                
            Arguments:
            filename -- A string, the filename of the pdb file.
            
            """

        # open/read the file
        afile = open(filename,"r")
        self.load_pdb_into_using_file_object(afile)
        afile.close()

    def load_pdb_into_using_file_object(self, file_obj):
        """Loads molecular data from a python file object (pdb formatted) into the current Molecule object. Note that most users will want to use the load_pdb_into() function instead, which is identical except that it accepts a filename string instead of a python file object.
                
            Arguments:
            file_obj -- A python file object, containing pdb-formatted data.
            
            """

        #source_data = numpy.genfromtxt(file_obj, dtype="S6,S5,S5,S4,S2,S4,S4,S8,S8,S8,S6,S6,S10,S2,S2", names=['record_name', 'serial', 'name', 'resname', 'chainid', 'resseq', 'empty', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'empty2', 'element', 'charge'], delimiter=[6, 5, 5, 4, 2, 4, 4, 8, 8, 8, 6, 6, 10, 2, 2])
        source_data = numpy.genfromtxt(file_obj, dtype="S6,S5,S5,S5,S1,S4,S4,S8,S8,S8,S6,S6,S10,S2,S3", names=['record_name', 'serial', 'name', 'resname', 'chainid', 'resseq', 'empty', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'empty2', 'element', 'charge'], delimiter=[6, 5, 5, 5, 1, 4, 4, 8, 8, 8, 6, 6, 10, 2, 3])
        
        if source_data.ndim == 0: source_data = source_data.reshape(1, -1) # in case the pdb file has only one line
        
        # get the ones that are ATOM or HETATOM in the record_name
        or_matrix = numpy.logical_or((source_data['record_name'] == "ATOM  "), (source_data['record_name'] == "HETATM"))
        indices_of_atom_or_hetatom = numpy.nonzero(or_matrix)[0]
        self.__parent_molecule.set_atom_information(source_data[indices_of_atom_or_hetatom])

        # now, some of the data needs to change types
        # first, fields that should be numbers cannot be empty strings
        for field in self.__parent_molecule.get_constants()['i8_fields'] + self.__parent_molecule.get_constants()['f8_fields']:
            check_fields = self.__parent_molecule.get_atom_information()[field]
            check_fields = numpy.core.defchararray.strip(check_fields)
            indices_of_empty = numpy.nonzero(check_fields == '')[0]
            self.__parent_molecule.get_atom_information()[field][indices_of_empty] = '0'
            
        # now actually change the type
        old_types = self.__parent_molecule.get_atom_information().dtype
        descr = old_types.descr
        for field in self.__parent_molecule.get_constants()['i8_fields']:
            index = self.__parent_molecule.get_atom_information().dtype.names.index(field)
            descr[index] = (descr[index][0], 'i8')
        for field in self.__parent_molecule.get_constants()['f8_fields']:
            index = self.__parent_molecule.get_atom_information().dtype.names.index(field)
            descr[index] = (descr[index][0], 'f8')
        new_types = numpy.dtype(descr)
        self.__parent_molecule.set_atom_information(self.__parent_molecule.get_atom_information().astype(new_types))
        
        # remove some of the fields that just contain empty data
        self.__parent_molecule.set_atom_information(self.__parent_molecule.numpy_structured_array_remove_field(self.__parent_molecule.get_atom_information(), ['empty', 'empty2']))
        
        # the coordinates need to be placed in their own special numpy array to facilitate later manipulation
        self.__parent_molecule.set_coordinates(numpy.vstack([self.__parent_molecule.get_atom_information()['x'], self.__parent_molecule.get_atom_information()['y'], self.__parent_molecule.get_atom_information()['z']]).T)
        self.__parent_molecule.set_atom_information(self.__parent_molecule.numpy_structured_array_remove_field(self.__parent_molecule.get_atom_information(), ['x', 'y', 'z'])) # now remove the coordinates from the atom_information object to save memory
        
        # string values in self.__parent_molecule.information.get_atom_information() should also be provided in stripped format for easier comparison
        fields_to_strip = ['name', 'resname', 'chainid', 'element']
        for f in fields_to_strip: self.__parent_molecule.set_atom_information(append_fields(self.__parent_molecule.get_atom_information(), f + '_stripped', data=numpy.core.defchararray.strip(self.__parent_molecule.get_atom_information()[f])))
        
class Selections():
    """A class for selecting atoms"""

    ######## selections ########
    def __init__(self, parent_molecule_object):
        """Initializes the Selections class.
                
            Arguments:
            parent_molecule_object -- The Molecule object associated with this class.
            
            """
        
        self.__parent_molecule = parent_molecule_object

    def select_atoms(self, selection_criteria):
        """Select a set of atoms based on user-specified criteria.
        
            Arguments:
            selection_criteria -- An dictionary, where the keys correspond to keys in the self.__parent_molecule.information.get_atom_information() structured numpy array, and the values are lists of acceptable matches.
                The selection is a logical "AND" between dictionary entries, but "OR" within the value lists themselves.
                For example: {'atom':['CA','O'], 'chain':'A', 'resname':'PRO'} would select all atoms with the names CA or O that are located in the PRO residues of chain A.
            
            Returns:
            A numpy.array containing the indices of the atoms of the selection.
            
            """
        
        try:
            selection = numpy.ones(len(self.__parent_molecule.get_atom_information()), dtype=bool) # start assuming everything is selected
            
            for key in selection_criteria.keys():
                
                vals = selection_criteria[key]
                
                # make sure the vals are in a list
                if not type(vals) is list and not type(vals) is tuple: vals = [vals] # if it's a single value, put it in a list
                
                # make sure the vals are in the right format
                if key in self.__parent_molecule.get_constants()['f8_fields']: vals = [float(v) for v in vals]
                elif key in self.__parent_molecule.get_constants()['i8_fields']: vals = [int(v) for v in vals]
                else: vals = [v.strip() for v in vals]
                
                # "or" all the vals together
                subselection = numpy.zeros(len(self.__parent_molecule.get_atom_information()), dtype=bool) # start assuming nothing is selected
                for val in vals: subselection = numpy.logical_or(subselection, (self.__parent_molecule.get_atom_information()[key] == val))
                
                # now "and" that with everything else
                selection = numpy.logical_and(selection, subselection)
            
            # now get the indices of the selection
            return numpy.nonzero(selection)[0]
        except:
            print "ERROR: Could not make the selection. Existing fields:"
            print "\t" + ", ".join(self.__parent_molecule.get_atom_information().dtype.names)
            sys.exit(0)
            
    def invert_selection(self, selection):
        """Inverts a user-defined selection (i.e., identifies all atoms that are not in the seleciton).
        
            Arguments:
            selection -- A numpy.array containing the indices of the user-defined selection.
            
            Returns:
            A numpy.array containing the indices of all atoms that are not in the user-defined seleciton.
            
            """

        # selection is a list of atom indices
        all_atoms = numpy.arange(0,len(self.__parent_molecule.get_atom_information()), 1, dtype=int)
        remaining_indicies = numpy.delete(all_atoms, selection)
        return remaining_indicies
    
    def select_all(self):
        """Selects all the atoms in a Molecule object.
        
            Returns:
            A numpy.array containing the indices of all atoms in the Molecule object.
            
            """

        return self.select_atoms({})

    def get_molecule_from_selection(self, selection):
        """Creates a Molecule from a user-defined atom selection.
        
            Arguments
            selection -- A numpy.array containing the indices of the atoms in the user-defined selection.
        
            Returns:
            A Molecule object containing the atoms of the user-defined selection.
        
            """

        new_mol = Molecule()
        new_mol.set_coordinates(self.__parent_molecule.get_coordinates()[selection])
        new_mol.set_atom_information(self.__parent_molecule.get_atom_information()[selection])
        
        # note that hierarchy will have to be recalculated
        
        return new_mol
    
# here's the actual Molecule class
class Molecule:
    """Loads, saves, and manupulates molecuar models. The main pymolecule class."""
    
    def __init__ (self):
        """Initializes the variables of the Molecule class."""
        
        self.fileio = FileIO(self)
        self.selections = Selections(self)
        self.information = Information(self)
  
    # Information methods
    def get_coordinates(self): return self.information.get_coordinates()
    def get_atom_information(self): return self.information.get_atom_information()
    def get_constants(self): return self.information.get_constants()
    def get_bounding_box(self, selection=None, padding=0.0): return self.information.get_bounding_box(selection, padding)
    def set_atom_information(self,atom_information): self.information.set_atom_information(atom_information)
    def set_coordinates(self,coordinates): self.information.set_coordinates(coordinates)

    # File I/O class methods
    def load_pdb_into(self, filename): self.fileio.load_pdb_into(filename)
    def load_pdb_into_using_file_object(self, file_obj): self.fileio.load_pdb_into_using_file_object(file_obj)
    
    # Selections class
    def get_molecule_from_selection(self, selection): return self.selections.get_molecule_from_selection(selection)
    def select_atoms(self, selection_criteria): return self.selections.select_atoms(selection_criteria)
    def invert_selection(self, selection): return self.selections.invert_selection(selection)
    def select_all(self): return self.selections.select_all()
    
    ######## supporting functions ########
    
    def numpy_structured_array_remove_field(self, narray, field_names): # surprised this doesn't come with numpy
        """Removes a specific field name from a structured numpy array.
                
            Arguments:
            narray -- A structured numpy array.
            field_names -- A list of strings, where each string is one of the field names of narray.
            
            Returns:
            A structured numpy array identical to narray, but with the field names in field_names removed.
            
            """
        
        names = list(narray.dtype.names) # now remove the coordinates from the atom_information object to save memory
        for f in field_names: names.remove(f)
        return narray[names]

# Some classes are required for calculating convex hulls

class ConvexHull():
    """A class to handle convex-hull calculations"""
    
    def __init__(self, pts):
        """Initializes the ConvexHull class."""

        akl_toussaint_pts = self.akl_toussaint(pts) # quickly reduces input size
        self.hull = self.gift_wrapping_3d(akl_toussaint_pts) # calculate convex hull using gift wrapping algorithm
    
    def inside_hull(self, our_point):
        """Determines if a point is inside the hull
            
            Arguments:
            our_point -- An x,y,z array
            
            Returns:
            A boolean, True if the point is inside the hull, False otherwise
            
            """
            
        return not self.outside_hull(our_point, self.hull)
    
    def outside_hull(self, our_point, triangles, epsilon=1.0e-5): # this one used internally
        """Given the hull as defined by a list of triangles, this definition will return whether a point is within these or not.
            
            Arguments:
            our_point -- an x,y,z array
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
        
        for i in range(n): # find the n with the largest x value
            point = tuple(raw_points[i])
            points.append(point)
            if point[0] > maxx:
                maxx = point[0]
                point1 = raw_points[i]
        
        best_dot = -1.0 # initiate dot relative to x-axis
        point2 = numpy.array(raw_points[1]) # initiate best segment
        
        # find first/best segment
        for i in range(n):
            pointi = raw_points[i]
            if numpy.array_equal(pointi, point1): continue
            diff_vec = pointi - point1
            diff_len = numpy.linalg.norm(diff_vec)
            
            test_dot = numpy.dot(diff_vec/diff_len,xaxis)
            if test_dot > best_dot:
                best_dot = test_dot
                point2 = pointi
        
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
                diff_vec1 = point2 - point1
                diff_vec2 = pointi - point2
                
                test_cross = numpy.array([diff_vec1[1]*diff_vec2[2]-diff_vec1[2]*diff_vec2[1], diff_vec1[2]*diff_vec2[0]-diff_vec1[0]*diff_vec2[2], diff_vec1[0]*diff_vec2[1]-diff_vec1[1]*diff_vec2[0]]) # cross product
                
                test_cross_len = numpy.sqrt(test_cross[0]*test_cross[0] + test_cross[1]*test_cross[1] + test_cross[2]*test_cross[2]) #numpy.linalg.norm(test_cross) # get the norm of the cross product
                
                if test_cross_len <= 0.0: continue
                test_cross = test_cross / test_cross_len
                dot_cross = numpy.dot(test_cross, ref_vec)
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

# Some classes are required for multiprocessing

class MultiThreading():
    """A class for multi-processor support."""

    results = []
    
    def __init__(self, inputs, num_processors, task_class_name):
        """Initializes the MultiThreading class."""

        self.results = []

        # first, if num_processors <= 0, determine the number of processors to use programatically
        if num_processors <= 0: num_processors = multiprocessing.cpu_count()

        # reduce the number of processors if too many have been specified
        if len(inputs) < num_processors: num_processors = len(inputs)

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
	    threads.append(task_class_name())
	    arrays.append(multiprocessing.Array('i',[0, 1]))

	results_queue = multiprocessing.Queue() # to keep track of the results

	processes = []
	for i in range(num_processors):
	    p = multiprocessing.Process(target=threads[i].runit, args=(running, mutex, results_queue, inputs_divided[i]))
	    p.start()
	    #p.join()
	    processes.append(p)

	while running.value > 0: is_running = 0 # wait for everything to finish
	
	# compile all results
	for thread in threads:
	    chunk =  results_queue.get()
	    self.results.extend(chunk)

class GeneralTask:
    """A class that determines the specific calculations that will be performed when multi-processor support is used. Other, more specific classes will inherit this one."""

    results = []
    
    def runit(self, running, mutex, results_queue, items):
        for item in items: self.value_func(item, results_queue)
        mutex.acquire()
        running.value -= 1
        mutex.release()
        results_queue.put(self.results)

    def value_func(self, item, results_queue): # this is the function that changes through inheritance
        print item # here's where you do something
        self.results.append(item) # here save the results for later compilation

# You'll also need a class representing a box of points, with associated definitions

class BoxOfPoints():
    """A class representing a box of equidistant points"""

    def __init__(self, box, reso):
        """Initialize the class.
    
            Arguments:
            box -- A numpy array representing two 3D points, (min_x, min_y, min_z) and (max_x, max_y, max_z), that define a box.
            reso -- The space between the points of the box, in the X, Y, and Z direction.

            """
        
        self.write_pdbs = write_pdbs()

        min_x = self.__snap_float(box[0][0], reso)
        min_y = self.__snap_float(box[0][1], reso)
        min_z = self.__snap_float(box[0][2], reso)
        max_x = self.__snap_float(box[1][0], reso) + 1.1 * reso
        max_y = self.__snap_float(box[1][1], reso) + 1.1 * reso
        max_z = self.__snap_float(box[1][2], reso) + 1.1 * reso

        x, y, z = numpy.mgrid[min_x:max_x:reso, min_y:max_y:reso, min_z:max_z:reso]
        self.points = numpy.array(zip(x.ravel(), y.ravel(), z.ravel()))
        
    def __snap_float(self, val, reso):
        """Snaps an arbitrary point to the nearest grid point.
    
            Arguments:
            val -- A numpy array corresponding to a 3D point.
            reso -- The resolution (distance in the X, Y, and Z directions between adjacent points) of the grid.
        
            Returns:
            A numpy array corresponding to a 3D point near val that is on a nearby grid point.
        
            """

        return numpy.floor(val / reso) * reso
    
    def remove_points_outside_convex_hull(self, hull):
        """Removes box points that are outside a convex hull.
    
            Arguments:
            hull -- The convex hull.
        
            """

        chunks = [(hull, t) for t in numpy.array_split(self.points, params['processors'])]
        tmp = MultiThreading(chunks, params['processors'], self.__MultiIdHullPts)
        self.points = numpy.vstack(tmp.results)

    class __MultiIdHullPts(GeneralTask):
        """A class to remove points outside a convex hull using multiple processors."""

        def value_func(self, items, results_queue): # so overwriting this function
            """The calculations that will run on a single processor to remove points outside a convex hull."""

            hull = items[0]
            some_points = items[1]
            
            # Note this would be much faster if it were matrix-based intead of point-by-point based.
            new_pts = [] # Can preallocate numpy array size because I don't know beforehand how many points will be in the hull
            for pt in some_points: 
                if hull.inside_hull(pt) == True: new_pts.append(pt)
                
            if len(new_pts) == 0: pass # here save the results for later compilation
            else: self.results.append(numpy.array(new_pts))

    def remove_all_points_close_to_other_points(self, other_points, dist_cutoff):
        """Removes all points in this box that come within the points specified in a numpy array
    
            Arguments:
            other_points -- A numpy array containing the other points.
            dist_cutoff -- A float, the cutoff distance to use in determining whether or not box points will be removed.
        
            """

        box_of_pts_distance_tree = spatial.KDTree(self.points) # note, in newer versions of scipy use cKDTree
        chunks = [(box_of_pts_distance_tree, dist_cutoff, t) for t in numpy.array_split(other_points, params['processors'])]
        tmp = MultiThreading(chunks, params['processors'], self.__MultiGetClosePoints)
        indicies_of_box_pts_close_to_molecule_points = numpy.unique(numpy.hstack(tmp.results))
        
        self.points = numpy.delete(self.points, indicies_of_box_pts_close_to_molecule_points, axis=0) # remove the ones that are too close to molecule atoms

    class __MultiGetClosePoints(GeneralTask):
        """A class to remove box points that are near other, user-specified points, using multiple processors."""

        def value_func(self, items, results_queue): # so overwriting this function
            """The calculations that will run on a single processor."""

            box_of_pts_distance_tree = items[0]
            dist_cutoff = items[1]
            other_points = items[2]
            
            other_points_distance_tree = spatial.KDTree(other_points) # note, in newer versions of scipy use cKDTree
            sparce_distance_matrix = other_points_distance_tree.sparse_distance_matrix(box_of_pts_distance_tree, dist_cutoff)
            indicies_of_box_pts_close_to_molecule_points = numpy.unique(sparce_distance_matrix.tocsr().indices) #tocsr()
            
            self.results.append(indicies_of_box_pts_close_to_molecule_points)
            
    def to_pdb(self, let="X"):
        """Converts the points in this box into a PDB representation.
    
            Arguments:
            let -- An optional string, the chain ID to use. "X" by default.

            Returns:
            A PDB-formatted string.

            """

        return self.write_pdbs.numpy_to_pdb(self.points, let)

    def expand_around_existing_points(self, num_pts, reso):
        """Add points to the current box that surround existing points, essentially increasing the resolution of the box.
    
            Arguments:
            num_pts -- An int, the number of points to place on each side of the existing points, in the X, Y, and Z directions.
            reso -- The distance between adjacent added points.

            """

        new_pts = []
        
        i = numpy.arange(-num_pts * reso, num_pts * reso + reso*0.01, reso)
        for xi in i:
            for yi in i:
                for zi in i:
                    vec = numpy.array([xi, yi, zi])
                    new_pts.append(self.points + vec)
        self.points = numpy.vstack(new_pts)
        
        self.__unique_points()
        
    def __unique_points(self):
        """Identifies unique points (rows) in an array of points.
        
        Arguments:
        a -- A nx3 numpy.array representing 3D points.
        
        Returns:
        A nx2 numpy.array containing the 3D points that are unique.
        
        """
        
        b = numpy.ascontiguousarray(self.points).view(numpy.dtype((numpy.void, self.points.dtype.itemsize * self.points.shape[1])))
        unique_points = numpy.unique(b).view(self.points.dtype).reshape(-1, self.points.shape[1])
        
        self.points = unique_points
        
    def filter_isolated_points_until_no_change(self, reso, number_of_neighbors):
        """Keep removing points that don't have enough neighbors, until no such points exist.
        
        Arguments:
        reso -- The distance between adjacent points.
        number_of_neighbors -- The minimum number of permissible neighbors. 
        
        """
        
        # calculate the pairwise distances between all box points
        
        box_of_pts_distance_tree = spatial.KDTree(self.points) # note, in newer versions of scipy use cKDTree
        print self.points
        self.dist_matrix = box_of_pts_distance_tree.sparse_distance_matrix(box_of_pts_distance_tree, reso * numpy.sqrt(3.0) * 1.1).todense() # so kiddy-corner counted as a neighbor
        
        # note that the diagnol of self.dist_matrix is zero, as expected, but ones with dist > reso * numpy.sqrt(3.0) * 1.1 are also 0. Pretty convenient.
        
        num_pts = 0
        while num_pts != len(self.points): # keep running the pass until there are no changes (points are stable)
            
            num_pts = len(self.points)
            
            # identify the points that have enough neighbors
            columns_nonzero_count = numpy.array((self.dist_matrix != 0).sum(0))[0]
            columns_nonzero_count_match_criteria = (columns_nonzero_count >= number_of_neighbors)
            columns_nonzero_count_match_criteria_index = numpy.nonzero(columns_nonzero_count_match_criteria)
            
            self.__keep_limited_points(columns_nonzero_count_match_criteria_index)
            
    def __keep_limited_points(self, pt_indices):
        """A support function"""
        
        # keep only those points
        self.points = self.points[pt_indices]

        # update the distance matrix so it doesn't need to be recalculated
        self.dist_matrix = self.dist_matrix[pt_indices,:][0]
        self.dist_matrix = self.dist_matrix.T
        self.dist_matrix = self.dist_matrix[pt_indices,:][0]
        #self.dist_matrix = self.dist_matrix.T # not necessary because it's a symetrical matrix
    
    def separate_out_pockets(self):
        """Separate the points according to the pocket they belong to. Determined by looking at patches of contiguous points.
        
        Returns:
        A list of point arrays, each array corresponding to the points of a separate pocket.
        
        """
        
        all_pockets = []
        
        while len(self.points) != 0:
        
            pocket_indexes = numpy.array([0])
            
            num_pts_in_pocket = 0
            
            while num_pts_in_pocket != len(pocket_indexes):
                
                num_pts_in_pocket = len(pocket_indexes)
                
                # get all the adjacent points
                pocket_indexes = numpy.hstack((pocket_indexes,numpy.array(numpy.nonzero(self.dist_matrix[pocket_indexes, :])[1])[0]))
                pocket_indexes = numpy.unique(pocket_indexes)
            
            pocket = self.points[pocket_indexes,:]
            all_pockets.append(pocket)
            
            self.__delete_limited_points(pocket_indexes)
        
        # sort the pockets by size
        all_pockets = sorted(all_pockets, key=lambda pts: -len(pts))
        
        return all_pockets
        
    def __delete_limited_points(self, pt_indices):
        """A support function"""
        
        # keep only those points
        self.points = numpy.delete(self.points, pt_indices, axis=0)

        # update the distance matrix so it doesn't need to be recalculated
        self.dist_matrix = numpy.delete(self.dist_matrix,pt_indices, axis=0)
        self.dist_matrix = self.dist_matrix.T
        self.dist_matrix = numpy.delete(self.dist_matrix,pt_indices, axis=0)

# Also, you need a class to save numpy arrays as PDB files

class write_pdbs():
    """A class for converting numpy arrays into PDB-formatted strings"""

    def __create_pdb_line(self, numpy_array, index, resname, letter):
        """Create a string formatted according to the PDB standard.
    
        Arguments:
        numpy_array -- A 1x3 numpy.array representing a 3D point.
        index -- An integer, the atom index to use in the string.
        resname -- A string, the RESNAME to use.
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
    
    def numpy_to_pdb(self, narray, letter, resname=""):
        """Create a string formatted according to the PDB standard.
    
        Arguments:
        narray -- A nx3 numpy.array representing a 3D point.
        letter -- A string, the atom name/chain/etc to use for the output.
        resname -- An optional string, the RESNAME to use for the output.
    
        Returns:
        A string, formatted according to the PDB standard.
    
        """
        
        if len(narray.flatten()) == 3:
            return self.__create_pdb_line(narray, 1, "AAA", letter) + "\n"
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
            for i, item in enumerate(narray): t = t + self.__create_pdb_line(item, i+1, resnames[i % len(resnames)], letter) + "\n"
            return t

####### Now the meat of the program ########

# First, show a brief help file describing the command-line arguments.

help_lines = []
help_lines.append('')
help_lines.append('POVME Pocket ID 1.0')
help_lines.append('===================')
help_lines.append('')
help_lines.append('Required command-line parameters:')
help_lines.append('')
help_lines.append('--filename: The PDB filename to be analyzed.')
help_lines.append('')
help_lines.append('Optional command-line parameters:')
help_lines.append('')
help_lines.append('--pocket_detection_resolution: The distance between probe points used to initially find the pockets (4.0 by default).')
help_lines.append('--pocket_measuring_resolution: The distance between probe points used to measure identified pockets in greater detail. Should divide --pocket_detection_resolution evenly. (1.0 by default).')
help_lines.append('--clashing_cutoff: In measuring the pockets, any points closer than this cutoff to receptor atoms will be removed. (3.0 by default).')
help_lines.append('--number_of_neighbors: In measuring the pockets, any points with fewer than this number of neighbors will be deleted. These are usually just stray points that don\'t belong to any real pocket. (4 by default).')
help_lines.append('--processors: The number of processors to use. (1 by default).')
help_lines.append('--number_of_spheres: The number of inclusion spheres to generate for each pocket. (5 by default).')
help_lines.append('--sphere_padding: How much larger the radius of the inclusion spheres should be, beyond what is required to encompass the identified pockets. (5.0 by default).')
help_lines.append('')
help_lines.append('Example:')
help_lines.append('')
help_lines.append('python pocket_id.py --filename rel1_example.pdb --pocket_detection_resolution 4.0 --pocket_measuring_resolution 1.0 --clashing_cutoff 3.0 --number_of_neighbors 4 --processors 1 --number_of_spheres 5 --sphere_padding 5.0 ')
help_lines.append('')

def printit(text): print textwrap.fill(text, initial_indent='', subsequent_indent='     ')

for line in help_lines: printit(line) 
if len(sys.argv[1:]) == 0: sys.exit(0)

# Now, parse the command-line arguments

params = {
    'filename': '',
    'pocket_detection_resolution': 4.0,
    'pocket_measuring_resolution': 1.0,
    'clashing_cutoff': 3.0,
    'number_of_neighbors': 4,
    'processors': 1,
    'number_of_spheres': 5,
    'sphere_padding': 5.0    
}

for item in getopt.getopt(sys.argv[1:], '', [ 'filename=', 'pocket_detection_resolution=', 'pocket_measuring_resolution=', 'clashing_cutoff=', 'number_of_neighbors=', 'processors=', 'number_of_spheres=', 'sphere_padding=' ])[0]: params[item[0].replace('--','')] = item[1]

if params['filename'] == '':
    print "ERROR: Must specify the --filename parameter!"
    print
    sys.exit(0)

for key in ['number_of_neighbors', 'processors', 'number_of_spheres']: params[key] = int(params[key])
for key in ['pocket_detection_resolution', 'pocket_measuring_resolution', 'clashing_cutoff', 'sphere_padding']: params[key] = float(params[key])

print 'Specified command-line arguments:'
print
for key in params: print "     --" + key + ': ' + str(params[key])
print

# Step 1: Load in the protein

printit("Step 1. Loading the PDB file " + params['filename'] + "...")
molecule = Molecule()
molecule.load_pdb_into(params['filename'])

# Step 2: Get rid of hydogen atoms. They just slow stuff down.

print "Step 2. Removing hydrogen atoms..."
sel = molecule.selections.select_atoms({'element_stripped':'H'})
sel = molecule.selections.invert_selection(sel)
molecule = molecule.selections.get_molecule_from_selection(sel)

# Step 3: Calculate the convex hull of the protein alpha carbons.
print "Step 3. Calculating the convex hull of the PDB file..."

molecule_alpha_carbons = molecule.selections.get_molecule_from_selection(molecule.selections.select_atoms({'name_stripped':'CA'})) # Get a version of the protein with just the alpha carbons. In my experience, that's better for convex hull identification. Otherwise the program identifies shallow contors in the protein surface as pockets.
convex_hull_3d = ConvexHull(molecule_alpha_carbons.get_coordinates())

# Step 4. Get a box of equispaced points that surround the protein, snapped to reso. I'm putting a whole bunch of other functions in this class as well to manipulate the points of this box.

printit("Step 4. Making a box of points spaced " + str(params['pocket_detection_resolution']) + " A apart that entirely encompasses the protein...")
    
box_pts = BoxOfPoints(molecule.get_bounding_box(), params['pocket_detection_resolution'] * 4) # note that the initial box is low resolution (* 4) so convex hull will be very fast

# Step 5. Remove points outside the convex hull. Gradually fill in protein-occupying region with denser point fields. Faster this way, I think.
printit("Step 5. Removing points that fall outside the protein's convex hull...")
box_pts.remove_points_outside_convex_hull(convex_hull_3d)
box_pts.expand_around_existing_points(2, params['pocket_detection_resolution'] * 2)
box_pts.remove_points_outside_convex_hull(convex_hull_3d)
box_pts.expand_around_existing_points(2, params['pocket_detection_resolution'])
box_pts.remove_points_outside_convex_hull(convex_hull_3d)

# Step 6. Remove the points in this box that are too close to protein atoms.
# For simplicity's sake, don't worry about atomic radii. Just a simple cutoff.
printit("Step 6. Removing points that come within " + str(params['clashing_cutoff']) + " A of any protein atom...")
box_pts.remove_all_points_close_to_other_points(molecule.get_coordinates(), params['clashing_cutoff'])

# Step 7. Now surround each of these points with higher density points that in the same regions. This is for getting a more detailed view of the identified pockets.
if params['pocket_measuring_resolution'] != params['pocket_detection_resolution']:
    printit("Step 7. Flooding the identified pockets with points spaced " + str(params['pocket_measuring_resolution']) + " A apart for a more detailed measurement of the pocket volume...")
    print "\tAdding points..."
    box_pts.expand_around_existing_points(params['pocket_detection_resolution']/params['pocket_measuring_resolution'], params['pocket_measuring_resolution'])
    printit("\tRemoving points that fall outside the convex hull...")
    box_pts.remove_points_outside_convex_hull(convex_hull_3d)
    printit("\tRemoving points within " + str(params['clashing_cutoff']) + " A of any protein atom...")
    box_pts.remove_all_points_close_to_other_points(molecule.get_coordinates(), params['clashing_cutoff'])
    
# Step 8. Now start doing a repeated pass filter (keep repeating until no change). Don't know if this is a high pass or low pass filter. I've heard these terms, though, and they sound cool.
printit("Step 8. Removing points until all points have at least " + str(params['number_of_neighbors']) + " neighbors...")
box_pts.filter_isolated_points_until_no_change(params['pocket_measuring_resolution'], params['number_of_neighbors'])

# Step 9. Separate out the pockets so they can be considered in isolation.
printit("Step 9. Partitioning the remaining points by pocket...")
all_pockets = box_pts.separate_out_pockets()

# Step 10. Get povme spheres that encompass each pocket, write pockets to seprate pdb files
printit("Step 10. Saving the points of each pocket...")
let_ids = ['A','B','C','D','E','F','G','H','I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W', 'X', 'Y', 'Z']
write_some_pdbs = write_pdbs()

for i,pts in enumerate(all_pockets):
    filename = 'pocket' + str(i+1) + '.pdb'
    printit("\tSaving " + filename + "...")
    f = open(filename,'w')
    f.write("REMARK Pocket #" + str(i+1) + "\n")
    
    # do I need to whiten stuff here? not sure what whitening is.
    
    centroids, idx = kmeans2(pts, params['number_of_spheres'])

    pts_string = ""    
    for cluster_num in range(params['number_of_spheres']):
        indexes_for_this_cluster = numpy.nonzero(idx == cluster_num)[0]
        cluster_pts = pts[indexes_for_this_cluster]
        cluster_center = numpy.mean(cluster_pts, axis=0)
        try:
            cluster_radius = numpy.max(cdist(numpy.array([cluster_center]), cluster_pts))
            f.write("REMARK CHAIN " + let_ids[cluster_num] + ": PointsInclusionSphere " + str(numpy.round(cluster_center[0],2)) + ' ' + str(numpy.round(cluster_center[1],2)) + ' ' + str(numpy.round(cluster_center[2],2)) + ' ' + str(numpy.round(cluster_radius + params['sphere_padding'],2)) + "\n")
            pts_string = pts_string + write_some_pdbs.numpy_to_pdb(cluster_pts, let_ids[cluster_num])
        except:
            print
            printit("There was an error, but I don't think it was catastrophic. Could be that one of the pocket clusters was empty.")
            print
            
    f.write(pts_string)
    f.close()

print
printit("Done. See the pocket{n}.pdb files. Using a visualization program like VMD, identify which of these files includes the pocket you wish to measure. POVME Pocket ID has divided each pocket volume into " + str(params['number_of_spheres']) + " sections (i.e., PDB chains). In some cases, the pocket you're interested in might be included in a larger identified pocket, so feel free to use only certain sections of a given pocket as well.")
printit("The POVME PointsInclusionSphere commands are located in the header of each pocket{n}.pdb file. A text editor can be used to copy and paste these commands into a POVME input file.")
print
