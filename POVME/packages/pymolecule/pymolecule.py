'''pymolecule is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    pymolecule is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
    
    Copyright 2011 Jacob D. Durrant. If you have any questions, comments, or
    suggestions, please don't hesitate to contact me at jdurrant [at] ucsd [dot] edu.
    
    The latest version of pymolecule can be downloaded from 
    http://sourceforge.net/projects/autoclickchem/
    
    If you use pymolecule in your work, please cite [REFERENCE HERE]'''

import os
import itertools
import numpy
import numpy.linalg
from numpy.lib.recfunctions import append_fields
from scipy.spatial.distance import pdist
from scipy.spatial.distance import cdist 
import scipy
import copy
import sys
import cPickle as pickle
import shutil
import warnings
warnings.filterwarnings("ignore", message="changing the sparsity structure of a csr_matrix is expensive. lil_matrix is more efficient.")

version="2.0"

if __name__ == '__main__': print "\npymolecule " + version + "\n"

class Information():
    '''A class for storing and accessing information about the elements of a pymolecule.Molecule object'''
    
    def __init__(self, parent_molecule_object):
        '''Initializes the pymolecule.Information class.
                
            Arguments:
            parent_molecule_object -- The pymolecule.Molecule object associated with this class.
            
            '''
        
        self.__parent_molecule = parent_molecule_object

        self.__constants = {}
        #self.__constants['element_names_with_two_letters'] = ['BR', 'CL', 'BI', 'AS', 'AG', 'LI', 'HG', 'MG', 'RH', 'ZN', 'MN']
        #Removed HG from this list to avoid capturing gamma hydrogens
        self.__constants['element_names_with_two_letters'] = ['BR', 'CL', 'BI', 'AS', 'AG', 'LI', 'MG', 'RH', 'ZN', 'MN'] 
        
        #SHORTEN LENGTH OF BOND_LENGTH_DICT
        
        self.__constants['bond_length_dict'] = {'C-C': 1.53, 'N-N': 1.425, 'O-O': 1.469, 'S-S': 2.048, 'C-H': 1.059, 'H-C': 1.059, 'C-N': 1.469, 'N-C': 1.469, 'C-O': 1.413, 'O-C': 1.413, 'C-S': 1.819, 'S-C': 1.819, 'N-H': 1.009, 'H-N': 1.009, 'N-O': 1.463, 'O-N': 1.463, 'O-S': 1.577, 'S-O': 1.577, 'O-H': 0.967, 'H-O': 0.967, 'S-H': 1.35, 'H-S': 1.35, 'S-N': 1.633, 'N-S': 1.633, 'C-F': 1.399, 'F-C': 1.399, 'C-CL': 1.790, 'CL-C': 1.790, 'C-BR': 1.910, 'BR-C': 1.910, 'C-I':2.162, 'I-C':2.162, 'S-BR': 2.321, 'BR-S': 2.321, 'S-CL': 2.283, 'CL-S': 2.283, 'S-F': 1.640, 'F-S': 1.640, 'S-I': 2.687, 'I-S': 2.687, 'P-BR': 2.366, 'BR-P': 2.366, 'P-CL': 2.008, 'CL-P': 2.008, 'P-F': 1.495, 'F-P': 1.495, 'P-I': 2.490, 'I-P': 2.490, 'P-C': 1.841, 'C-P': 1.841, 'P-N': 1.730, 'N-P': 1.730, 'P-O': 1.662, 'O-P': 1.662, 'P-S': 1.954, 'S-P': 1.954, 'N-BR': 1.843, 'BR-N': 1.843, 'N-CL': 1.743, 'CL-N': 1.743, 'N-F': 1.406, 'F-N': 1.406, 'N-I': 2.2, 'I-N': 2.2, 'SI-BR': 2.284, 'BR-SI': 2.284, 'SI-CL': 2.072, 'CL-SI': 2.072, 'SI-F': 1.636, 'F-SI': 1.636, 'SI-P': 2.264, 'P-SI': 2.264, 'SI-S': 2.145, 'S-SI': 2.145, 'SI-SI': 2.359, 'SI-SI': 2.359, 'SI-C': 1.888, 'C-SI': 1.888, 'SI-N': 1.743, 'N-SI': 1.743, 'SI-O': 1.631, 'O-SI': 1.631, 'X-X': 1.53, 'X-C': 1.53, 'C-X': 1.53, 'X-H': 1.059, 'H-X': 1.059, 'X-N': 1.469, 'N-X': 1.469, 'X-O': 1.413, 'O-X': 1.413, 'X-S': 1.819, 'S-X': 1.819, 'X-F': 1.399, 'F-X': 1.399, 'X-CL': 1.790, 'CL-X': 1.790, 'X-BR': 1.910, 'BR-X': 1.910, 'X-I':2.162, 'I-X':2.162, 'SI-X': 1.888, 'X-SI': 1.888, 'H-H':0.74}
        
        #INCLUDE ALL N AND C TERMINAL
        
        self.__constants['protein_residues'] = ["ALA", "ARG", "ASH", "ASN", "ASP", "CYM", "CYS", "CYX", "GLN", "GLH", "GLU", "GLY", "HID", "HIE", "HIP", "HIS", "ILE", "LEU", "LYN", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL","MSE","TPO","PTR","SEP", 'CALA', 'CARG', 'CASN', 'CASP', 'CCYS', 'CCYX', 'CGLN', 'CGLU', 'CGLY', 'CHID', 'CHIE', 'CHIP', 'CHIS', 'CILE', 'CLEU', 'CLYS', 'CMET', 'CPHE', 'CPRO', 'CSER', 'CTHR', 'CTRP', 'CTYR', 'CVAL', 'NALA', 'NARG', 'NASN', 'NASP', 'NCYS', 'NCYX', 'NGLN', 'NGLU', 'NGLY', 'NHID', 'NHIE', 'NHIP', 'NHIS', 'NILE', 'NLEU', 'NLYS', 'NMET', 'NPHE', 'NPRO', 'NSER', 'NTHR', 'NTRP', 'NTYR', 'NVAL']
        self.__constants['dna_residues'] = ["A", "C", "G", "T", "DA", "DA3", "DA5", "DAN", "DC", "DC3", "DC4", "DC5", "DCN", "DG", "DG3", "DG5", "DGN", "DT", "DT3", "DT5", "DTN"]
        self.__constants['rna_residues'] = ["A", "C", "G", "U", "RA", "RA3", "RA5", "RAN", "RC", "RC3", "RC4", "RC5", "RCN", "RG", "RG3", "RG5", "RGN", "RU", "RU3", "RU5", "RUN"]
        self.__constants['mass_dict'] = {'H':  1.00794, 'C':  12.0107, 'N':  14.0067, 'O':  15.9994, 'S':  32.065, 'P':  30.973762, 'NA':  22.9897, 'MG':  24.3050, 'F':  18.9984032, 'CL':  35.453, 'K':  39.0983, 'CA':  40.078, 'I':  126.90447, 'LI':  6.941, 'BE':  9.0122, 'B':  10.811, 'AL':  26.9815, 'MN':  54.938, 'FE':  55.845, 'CO':  58.9332, 'CU':  63.9332, 'ZN':  65.38, 'AS':  74.9216, 'BR':  79.904, 'MO':  95.94, 'RH':  102.9055, 'AG':  107.8682, 'AU':  196.9655, 'HG':  200.59, 'PB':  207.2, 'BI':  208.98040}
        self.__constants['vdw_dict'] = {'H':  1.2, 'C':  1.7, 'N':  1.55, 'O':  1.52, 'F':  1.47, 'P':  1.8, 'S':  1.8, 'B':  2.0, 'LI':  1.82, 'NA':  2.27, 'MG':  1.73, 'AL':  2.00, 'CL':  1.75, 'CA':  2.00, 'MN':  2.00, 'FE':  2.00, 'CO':  2.00, 'CU':  1.40, 'ZN':  1.39, 'AS':  1.85, 'BR':  1.85, 'MO':  2.00, 'RH':  2.00, 'AG':  1.72, 'AU':  1.66, 'PB':  2.02, 'BI':  2.00, 'K':  2.75, 'I':  1.98}
        self.__constants['i8_fields'] = ['serial','resseq']
        self.__constants['f8_fields']= ['x','y','z','occupancy','tempfactor']
        self.__constants['max_number_of_bonds_permitted'] = {"C": 4, "N": 4, "O": 2, "H": 1, "F": 1, "Cl": 1, "BR": 1, "CL": 1, "I": 1, "P": 5, "S": 6}

        self.__filename = ""
        self.__remarks = []
        self.__atom_information = None
        self.__coordinates = None
        self.__coordinates_undo_point = None
        self.__bonds = None
        self.__hierarchy = {}
        self.__max_ring_size = 50

    def get_filename(self): return self.__filename
    def get_remarks(self): return self.__remarks
    def get_atom_information(self): return self.__atom_information
    def get_coordinates(self): return self.__coordinates
    def get_coordinates_undo_point(self): return self.__coordinates_undo_point
    def get_bonds(self): return self.__bonds
    def get_hierarchy(self): return self.__hierarchy
    def get_constants(self): return self.__constants
    
    def set_filename(self,filename): self.__filename = filename
    def set_remarks(self,remarks): self.__remarks = remarks
    def set_atom_information(self,atom_information): self.__atom_information = atom_information
    def set_coordinates(self,coordinates): self.__coordinates = coordinates
    def set_coordinates_undo_point(self,coordinates_undo_point): self.__coordinates_undo_point = coordinates_undo_point
    def set_bonds(self,bonds): 
        if bonds is None:
            self.__bonds = None
        else:
            self.__bonds = scipy.sparse.csr_matrix(bonds)
    def set_hierarchy(self,hierarchy): self.__hierarchy = hierarchy

    def belongs_to_protein(self, atom_index):
        '''Checks if the atom is part of a protein.  Taken primarily from Amber residue names.
            
            Arguments:
            atom_index -- An int, the index of the atom to consider.

            Returns:
            A boolean. True if part of protein, False if not.
            
            '''
        
        # this function is retained for legacy reasons. past versions of pymolecule had
        # this functionality.

        if self.__atom_information['resname_stripped'][atom_index] in self.__constants['protein_residues']: return True
        return False

    def belongs_to_dna(self, atom_index):
        '''Checks if the atom is part of DNA.
            
            Arguments:
            atom_index -- An int, the index of the atom to consider.

            Returns:
            A boolean. True if part of dna, False if not.
            
            '''

        # this function is retained for legacy reasons. past versions of pymolecule had
        # this functionality.
        
        if self.__atom_information['resname_stripped'][atom_index] in self.__constants['dna_residues']: return True
        return False

    def belongs_to_rna(self, atom_index):
        '''Checks if the atom is part of RNA.
            
            Arguments:
            atom_index -- An int, the index of the atom to consider.

            Returns:
            A boolean. True if part of rna, False if not.
            
            '''

        # this function is retained for legacy reasons. past versions of pymolecule had
        # this functionality.
        
        if self.__atom_information['resname_stripped'][atom_index] in self.__constants['rna_residues']: return True
        return False

    def assign_masses(self):
        '''Assigns masses to the atoms of the pymolecule.Molecule object.'''

        if not "mass" in self.__atom_information.dtype.names: # only assign if not been assigned previously
            masses = numpy.empty((len(self.__atom_information['element_stripped'])))
            for i in range(len(self.__atom_information['element_stripped'])):
                element = self.__atom_information['element_stripped'][i]
                mass = self.__constants['mass_dict'][element]
                masses[i] = mass
            
            self.__atom_information = append_fields(self.__atom_information, 'mass', data=masses)
            
    def assign_elements_from_atom_names(self, selection=None):
        '''Determines the elements of all atoms from the atom names. Note that this will overwrite any existing element assignments, including those explicitly specified in loaded files. Note that this doesn't populate elements_stripped.

            Arguments:
            selection -- An optional numpy.array containing the indices of the atoms to consider when calculating the center of mass. If ommitted, all atoms of the pymolecule.Molecule object will be considered.

        '''

        if selection is None: selection = self.__parent_molecule.select_all()

        if len(selection) == 0: return

        # get the atom names
        fix_element_names = numpy.core.defchararray.upper(self.__atom_information['name'][selection])
        fix_element_names = numpy.core.defchararray.strip(fix_element_names)
        
        # first remove any numbers at the begining of these names
        fix_element_names = numpy.core.defchararray.lstrip(fix_element_names,'0123456789')
        
        # remove any thing, letters or numbers, that follows a number, including the number itself. so C2L becomes C, not CL.
        for num in ['0', '1', '2', '3', '4', '5', '6', '7', '8', '9']: # I wish there was a more numpified way of doing this. :(
            tmp = numpy.core.defchararray.split(fix_element_names, num)
            fix_element_names = numpy.empty(len(fix_element_names), dtype="S5")
            for i, item in enumerate(tmp): fix_element_names[i] = tmp[i][0]
        
        # take just first two letters of each item
        fix_element_names = numpy.array(fix_element_names, dtype="|S2")

        # identify ones that are two-letter elements and one-letter elements
        one_that_should_be_two_letters = (fix_element_names == self.__constants['element_names_with_two_letters'][0])
        for other_two_letter in self.__constants['element_names_with_two_letters'][1:]: one_that_should_be_two_letters = numpy.logical_or(one_that_should_be_two_letters, (fix_element_names == other_two_letter))
        indices_of_two_letter_elements = numpy.nonzero(one_that_should_be_two_letters)[0]
        indices_of_one_letter_elements = numpy.nonzero(numpy.logical_not(one_that_should_be_two_letters))[0]

        # get ones that are one-letter elements
        fix_element_names[indices_of_one_letter_elements] = numpy.core.defchararray.rjust(numpy.array(fix_element_names[indices_of_one_letter_elements], dtype="|S1"),2)
        
        # they should be capitalized for consistency
        fix_element_names = numpy.core.defchararray.upper(fix_element_names)
        
        # now map missing element names back
        self.__atom_information['element'][selection] = fix_element_names  #This line throws an imminent-deprecation warning on 12/18/2014
        ### This doesn't fix it
        #self.__atom_information['element'] = fix_element_names  #This line throws an imminent-deprecation warning on 12/18/2014
        ### This also doesn't fix it
        #for index in selection:
        #    self.__atom_information['element'][index] = fix_element_names[index]
        
        
        
        # element_stripped also needs to be updated
        #try: self.__parent_molecule.information.get_atom_information()['element_stripped'][selection] = numpy.core.defchararray.strip(fix_element_names)
        #except: # so element_stripped hasn't been defined yet
        #    self.__parent_molecule.information.get_atom_information() = append_fields(self.__parent_molecule.information.get_atom_information(), 'element_stripped', data=numpy.core.defchararray.strip(self.__parent_molecule.information.get_atom_information()['element']))
    
    def get_center_of_mass(self, selection=None):
        '''Determines the center of mass.

            Arguments:
            selection -- An optional numpy.array containing the indices of the atoms to consider when calculating the center of mass. If ommitted, all atoms of the pymolecule.Molecule object will be considered.

            Returns:
            A numpy.array containing to the x, y, and z coordinates of the center of mass.
        
            '''
        
        if selection is None: selection = self.__parent_molecule.select_all()
        
        # make sure the masses have been asigned
        self.assign_masses()
            
        # calculate the center of mass
        
        # multiply each coordinate by its mass
        center_of_mass = self.__coordinates[selection] * numpy.vstack((self.__atom_information['mass'][selection], self.__atom_information['mass'][selection], self.__atom_information['mass'][selection])).T
        
        # now sum all that
        center_of_mass = numpy.sum(center_of_mass,0)
        
        # now divide by the total mass
        center_of_mass = center_of_mass / self.get_total_mass(selection)
        
        return center_of_mass
        
    def get_geometric_center(self, selection=None):
        '''Determines the geometric center.
            
            Arguments:
            selection -- An optional numpy.array containing the indices of the atoms to consider when calculating the geometric center. If ommitted, all atoms of the pymolecule.Molecule object will be considered.

            Returns:
            A numpy.array containing to the x, y, and z coordinates of the geometric center.

            '''
        
        if selection is None: selection = self.__parent_molecule.select_all()
        
        return numpy.sum(self.__coordinates[selection],0) / self.get_total_number_of_atoms(selection)
    
    def get_total_mass(self, selection=None):
        '''Calculates the total atomic mass.
            
            Arguments:
            selection -- An optional numpy.array containing the indices of the atoms to consider when calculating the total mass. If ommitted, all atoms of the pymolecule.Molecule object will be considered.

            Returns:
            A double, the total mass.

            '''

        if selection is None: selection = self.__parent_molecule.select_all()

        # assign masses if necessary
        self.assign_masses()
        
        # return total mass
        return numpy.sum(self.__atom_information['mass'][selection])

    def get_total_number_of_atoms(self, selection=None):
        '''Counts the number of atoms.
            
            Arguments:
            selection -- An optional numpy.array containing the indices of the atoms to count. If ommitted, all atoms of the pymolecule.Molecule object will be considered.

            Returns:
            An int, the total number of atoms.

            '''

        if selection is None: selection = self.__parent_molecule.select_all()

        if self.__coordinates is None: return 0
        else: return len(self.__coordinates[selection])

    def get_total_number_of_heavy_atoms(self):
        '''Counts the number of heavy atoms (i.e., atoms that are not hydrogens).
            
            Arguments:
            selection -- An optional numpy.array containing the indices of the atoms to count. If ommitted, all atoms of the pymolecule.Molecule object will be considered.

            Returns:
            An int, the total number of heavy (non-hydrogen) atoms.

            '''

        if self.__coordinates is None: return 0
        
        all_hydrogens = self.__parent_molecule.select_atoms({'element_stripped':'H'})
        
        return self.get_total_number_of_atoms() - len(all_hydrogens)
    
    def get_bounding_box(self, selection = None, padding=0.0):
        '''Calculates a box that bounds (encompasses) a set of atoms.
            
            Arguments:
            selection -- An optional numpy.array containing the indices of the atoms to consider. If ommitted, all atoms of the pymolecule.Molecule object will be considered.
            padding -- An optional float. The bounding box will extend this many angstroms beyond the atoms being considered.
            
            Returns:
            A numpy array representing two 3D points, (min_x, min_y, min_z) and (max_x, max_y, max_z), that bound the molecule.
            
            '''
        
        if selection is None: selection = self.__parent_molecule.select_all()
        
        return numpy.vstack((numpy.min(self.__coordinates[selection],0), numpy.max(self.__coordinates[selection],0)))
    
    def get_bounding_sphere(self, selection=None, padding=0.0):
        '''Calculates a sphere that bounds (encompasses) a set of atoms.
            
            Arguments:
            selection -- An optional numpy.array containing the indices of the atoms to consider. If ommitted, all atoms of the pymolecule.Molecule object will be considered.
            padding -- An optional float. The bounding sphere will extend this many angstroms beyond the atoms being considered.
            
            Returns:
            A tuple containing two elements. The first is a numpy.array representing a 3D point, the (x, y, z) center of the sphere. The second is a float, the radius of the sphere.
            
            '''
        
        if selection is None: selection = self.__parent_molecule.select_all()
        
        # get center
        center_of_selection = numpy.array([self.get_geometric_center(selection)])
        
        # get distance to farthest point in selection
        return (center_of_selection[0], numpy.max(cdist(center_of_selection, self.__coordinates[selection])[0]))

    def define_molecule_chain_residue_spherical_boundaries(self):
        '''Identifies spheres that bound (encompass) the entire molecule, the chains, and the residues. This information is stored in pymolecule.Molecule.information.hierarchy.'''

        # first, check to see if it's already been defined
        if 'spheres' in self.__hierarchy.keys(): return
        
        # set up the new structure
        self.__hierarchy['spheres'] = {}
        self.__hierarchy['spheres']['molecule'] = {}
        self.__hierarchy['spheres']['chains'] = {}
        self.__hierarchy['spheres']['residues'] = {}

        # get all the chains and residues
        chains = self.__parent_molecule.selections_of_chains()
        residues = self.__parent_molecule.selections_of_residues()
        
        # do calcs for the whole molcules
        whole_mol_calc = self.get_bounding_sphere()
        self.__hierarchy['spheres']['molecule']['center'] = numpy.array([whole_mol_calc[0]])
        self.__hierarchy['spheres']['molecule']['radius'] = whole_mol_calc[1]
        
        # do calcs for the chains
        self.__hierarchy['spheres']['chains']['keys'] = numpy.array(self.__hierarchy['chains']['indices'].keys()) #numpy string array e.g. ['a','b','c']
        self.__hierarchy['spheres']['chains']['centers'] = numpy.empty((len(self.__hierarchy['spheres']['chains']['keys']),3))
        self.__hierarchy['spheres']['chains']['radii'] = numpy.empty(len(self.__hierarchy['spheres']['chains']['keys']))
        
        for index, chainid in enumerate(self.__hierarchy['spheres']['chains']['keys']):
            asphere  = self.get_bounding_sphere(selection=self.__hierarchy['chains']['indices'][chainid])
            self.__hierarchy['spheres']['chains']['centers'][index][0] = asphere[0][0]
            self.__hierarchy['spheres']['chains']['centers'][index][1] = asphere[0][1]
            self.__hierarchy['spheres']['chains']['centers'][index][2] = asphere[0][2]
            self.__hierarchy['spheres']['chains']['radii'][index] = asphere[1]
            
        # do calcs for the residues
        self.__hierarchy['spheres']['residues']['keys'] = numpy.array(self.__hierarchy['residues']['indices'].keys())
        self.__hierarchy['spheres']['residues']['centers'] = numpy.empty((len(self.__hierarchy['spheres']['residues']['keys']),3))
        self.__hierarchy['spheres']['residues']['radii'] = numpy.empty(len(self.__hierarchy['spheres']['residues']['keys']))
        
        for index, resid in enumerate(self.__hierarchy['spheres']['residues']['keys']):
            asphere  = self.get_bounding_sphere(selection=self.__hierarchy['residues']['indices'][resid])
            self.__hierarchy['spheres']['residues']['centers'][index][0] = asphere[0][0]
            self.__hierarchy['spheres']['residues']['centers'][index][1] = asphere[0][1]
            self.__hierarchy['spheres']['residues']['centers'][index][2] = asphere[0][2]
            self.__hierarchy['spheres']['residues']['radii'][index] = asphere[1]

    def serial_reindex(self):
        '''Reindexes the serial field of the atoms in the molecule, starting with 1'''
        
        for i in range(len(self.__atom_information['serial'])): self.__atom_information['serial'][i] = i + 1

    def resseq_reindex(self):
        '''Reindexes the resseq field of the atoms in the molecule, starting with 1'''
        
        keys = numpy.core.defchararray.add(self.__atom_information['resname_stripped'], '-')
        keys = numpy.core.defchararray.add(keys, numpy.array([str(t) for t in self.__atom_information['resseq']]))
        keys = numpy.core.defchararray.add(keys, '-')
        keys = numpy.core.defchararray.add(keys, self.__atom_information['chainid_stripped'])
        
        keys2 = numpy.insert(keys,0, '')[:-1]
        index_of_change = numpy.nonzero(numpy.logical_not(keys == keys2))[0]
        index_of_change = numpy.append(index_of_change, len(self.__atom_information))
        
        count = 1
        for t in range(len(index_of_change[:-1])):
            start = index_of_change[t]
            end = index_of_change[t+1]
            self.__atom_information['resseq'][numpy.arange(start,end,1,dtype='int')] = count
            count = count + 1

class FileIO():
    '''A class for saving and loading molecular data into a pymolecule.Molecule object'''
    
    def __init__(self, parent_molecule_object):
        '''Initializes the pymolecule.FileIO class.
                
            Arguments:
            parent_molecule_object -- The pymolecule.Molecule object associated with this class.
            
            '''
        
        self.__parent_molecule = parent_molecule_object
        
    def load_pym_into(self, filename):
        '''Loads the molecular data contained in a pym file into the current pymolecule.Molecule object.
                
            Arguments:
            filename -- A string, the filename of the pym file.
            
            '''

        if filename[-1:] != os.sep: filename = filename + os.sep
        
        # first, get the files that must exist
        self.__parent_molecule.set_atom_information(pickle.load(open(filename + 'atom_information', "rb" )))
        self.__parent_molecule.set_coordinates(numpy.load(filename + "coordinates.npz")['arr_0'])
        # now look for other possible files (optional output)
        if os.path.exists(filename + 'remarks'): self.__parent_molecule.set_remarks(pickle.load( open( filename + 'remarks', "rb" ) ) ) 
        if os.path.exists(filename + 'hierarchy'): self.__parent_molecule.set_hierarchy(pickle.load( open( filename + 'hierarchy', "rb" ) ) )
        if os.path.exists(filename + 'filename'): self.__parent_molecule.set_filename(pickle.load( open( filename + 'filename', "rb" ) ) )
        if os.path.exists(filename + "bonds.npz"):
            # Some extra work here to load sparse bond matrix
            bonds_raw = numpy.load(filename + "bonds.npz") #['arr_0']
            bonds_sparse = scipy.sparse.csr_matrix((bonds_raw['data'],
                                                    bonds_raw['indices'],
                                                    bonds_raw['indptr']),
                                                   shape=bonds_raw['shape'])
            self.__parent_molecule.set_bonds(bonds_sparse)
        if os.path.exists(filename + "coordinates_undo_point.npz"): self.__parent_molecule.set_coordinates_undo_point(numpy.load(filename + "coordinates_undo_point.npz")['arr_0'])
    
    def load_pdb_into(self, filename, bonds_by_distance=True, serial_reindex=True, resseq_reindex=False):
        '''Loads the molecular data contained in a pdb file into the current pymolecule.Molecule object.
                
            Arguments:
            filename -- A string, the filename of the pdb file.
            bonds_by_distance -- An optional boolean, whether or not to determine atomic bonds based on atom proximity. True by default.
            serial_reindex -- An optional boolean, whether or not to reindex the pdb serial field. True by default.
            resseq_reindex -- An optional boolean, whether or not to reindex the pdb resseq field. False by default.
            
            '''

        self.__parent_molecule.set_filename(filename)
        
        # open/read the file
        afile = open(filename,"r")
        self.load_pdb_into_using_file_object(afile, bonds_by_distance, serial_reindex, resseq_reindex)
        afile.close()

    def load_pdb_into_using_file_object(self, file_obj, bonds_by_distance=True, serial_reindex=True, resseq_reindex=False):
        '''Loads molecular data from a python file object (pdb formatted) into the current pymolecule.Molecule object. Note that most users will want to use the load_pdb_into() function instead, which is identical except that it accepts a filename string instead of a python file object.
                
            Arguments:
            file_obj -- A python file object, containing pdb-formatted data.
            bonds_by_distance -- An optional boolean, whether or not to determine atomic bonds based on atom proximity. True by default.
            serial_reindex -- An optional boolean, whether or not to reindex the pdb serial field. True by default.
            resseq_reindex -- An optional boolean, whether or not to reindex the pdb resseq field. False by default.
            
            '''

        #source_data = numpy.genfromtxt(file_obj, dtype="S6,S5,S5,S4,S2,S4,S4,S8,S8,S8,S6,S6,S10,S2,S2", names=['record_name', 'serial', 'name', 'resname', 'chainid', 'resseq', 'empty', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'empty2', 'element', 'charge'], delimiter=[6, 5, 5, 4, 2, 4, 4, 8, 8, 8, 6, 6, 10, 2, 2])
        source_data = numpy.genfromtxt(file_obj, dtype="S6,S5,S5,S5,S1,S4,S4,S8,S8,S8,S6,S6,S10,S2,S3", names=['record_name', 'serial', 'name', 'resname', 'chainid', 'resseq', 'empty', 'x', 'y', 'z', 'occupancy', 'tempfactor', 'empty2', 'element', 'charge'], delimiter=[6, 5, 5, 5, 1, 4, 4, 8, 8, 8, 6, 6, 10, 2, 3])
        
        # get the remarks, if any. good to hold on to this because some of my programs might retain info via remarks
        remark_indices = numpy.nonzero(source_data['record_name'] == "REMARK")[0]
        remarks = []
        for index in remark_indices:
            astr = ""
            for name in source_data.dtype.names[1:]: astr = astr + source_data[name][index]
            remarks.append(astr.rstrip())
        self.__parent_molecule.set_remarks(remarks)
        
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
        
        # now determine element from atom name for those entries where it's not given
        # note that the molecule.information.assign_elements_from_atom_names function can be used to overwrite this and
        # assign elements based on the atom name only.
        indicies_where_element_is_not_defined = numpy.nonzero(numpy.core.defchararray.strip(self.__parent_molecule.get_atom_information()['element']) == '')[0]
        
        self.__parent_molecule.assign_elements_from_atom_names(indicies_where_element_is_not_defined)
        
        # string values in self.__parent_molecule.information.get_atom_information() should also be provided in stripped format for easier comparison
        fields_to_strip = ['name', 'resname', 'chainid', 'element']
        for f in fields_to_strip: self.__parent_molecule.set_atom_information(append_fields(self.__parent_molecule.get_atom_information(), f + '_stripped', data=numpy.core.defchararray.strip(self.__parent_molecule.get_atom_information()[f])))
        
        # now, if there's conect data, load it. this part of the code is not that "numpyic"
        conect_indices = numpy.nonzero(source_data['record_name'] == "CONECT")[0]
        if len(conect_indices) > 0:
            
            self.__parent_molecule.set_bonds(numpy.zeros((len(self.__parent_molecule.get_atom_information()), len(self.__parent_molecule.get_atom_information()))))
            
            # build serial to index mapping
            serial_to_index = {}
            for index, inf in enumerate(self.__parent_molecule.get_atom_information()['serial']): serial_to_index[inf] = index # is there a faster way?
            
            # get the connect data
            for index in conect_indices:
                astr = ""
                for name in source_data.dtype.names[1:]: astr = astr + source_data[name][index]
                astr= astr.rstrip()
                
                indices = []
                for i in xrange(0, len(astr), 5): indices.append(serial_to_index[int(astr[i:i+5])])
                
                for partner_index in indices[1:]:
                    self.__parent_molecule.add_bond(indices[0],partner_index)

        #else: # create empty bond array
        #    self.__parent_molecule.information.get_bonds() = numpy.zeros((len(self.__parent_molecule.information.get_atom_information()), len(self.__parent_molecule.information.get_atom_information())))

        if bonds_by_distance == True: self.__parent_molecule.create_bonds_by_distance(True)
        if serial_reindex == True: self.__parent_molecule.serial_reindex()
        if resseq_reindex == True: self.__parent_molecule.resseq_reindex()
        
    def save_pym(self, filename, save_bonds=False, save_filename=False, save_remarks=False, save_hierarchy=False, save_coordinates_undo_point=False):
        '''Saves the molecular data contained in a pymolecule.Molecule object to a pym file.
                
            Arguments:
            filename -- An string, the filename to use for saving. (Note that this is actually a directory, not a file.)
            save_bonds -- An optional boolean, whether or not to save information about atomic bonds. False by default.
            save_filename -- An optional boolean, whether or not to save the original (pdb) filename. False by default.
            save_remarks -- An optional boolean, whether or not to save remarks associated with the molecule. False by default.
            save_hierarchy -- An optional boolean, whether or not to save information about spheres the bound (encompass) the whole molecule, the chains, and the residues. False by default.
            save_coordinates_undo_point -- An optional boolean, whether or not to save the last coordinate undo point. False by default.
            
            '''
            
        # Why not just pickle self.parent.information? Because it's a huge file,
        # can't selectively not save bonds, for example, 
        # and numpy.save is faster than cPickle protocol 2 on numpy arrays
        
        # if the directory already exists, first delete it
        if os.path.exists(filename):
            try: shutil.rmtree(filename)
            except: pass
            
            # it could be a file, not a directory
            try: os.remove(filename)
            except: pass
        
        # filename is actually a directory, so append separator if needed
        if filename[-1:] != os.sep: filename = filename + os.sep

        # make directory
        os.mkdir(filename)
        
        # save components
        
        # python objects must be pickled
        if save_hierarchy == True: pickle.dump(self.__parent_molecule.get_hierarchy(), open(filename + 'hierarchy','wb'), -1) # note this is a combo of python objects and numpy arrays, so must be pickled.
        if save_remarks == True: pickle.dump(self.__parent_molecule.get_remarks(), open(filename + 'remarks','wb'), -1) # using the latest protocol
        if save_filename == True: pickle.dump(self.__parent_molecule.get_filename(), open(filename + 'filename','wb'), -1)
        
        # unfortunately, the speedy numpy.save doesn't work on masked arrays
        # masked arrays have a dump method, but it just uses cPickle
        # so we're just going to cPickle masked arrays. Could be so much faster if numpy were up to speed... :(
        # not clear that numpy.ma.dump accepts protocol parameter, so let's just use cPickle directly
        pickle.dump(self.__parent_molecule.get_atom_information(), open(filename + 'atom_information','wb'), -1)

        # fortunately, coordinates and bonds are regular numpy arrays
        # they can be saved with numpy's speedy numpy.save function
        # note that I'm compressing them here. benchmarking suggests
        # this takes longer to save, but is much faster to load.
        # so I'm prioritizing load times over save times
        # note also that numpy.savez can save multiple arrays to a single file,
        # probably speeding up load.
        
        numpy.savez(filename + "coordinates.npz", self.__parent_molecule.get_coordinates())
        if save_bonds == True:
            bonds = self.__parent_molecule.get_bonds()
            numpy.savez(filename + "bonds.npz",
                        data=bonds.data,
                        indices=bonds.indices,
                        indptr=bonds.indptr,
                        shape=bonds.shape)
        if save_coordinates_undo_point == True: numpy.savez(filename + "coordinates_undo_point.npz", self.__parent_molecule.get_coordinates_undo_point())
        
    def save_pdb(self, filename="", serial_reindex=True, resseq_reindex=False, return_text = False):
        '''Saves the molecular data contained in a pymolecule.Molecule object to a pdb file.
                
            Arguments:
            filename -- A string, the filename to use for saving.
            serial_reindex -- An optional boolean, whether or not to reindex the pdb serial field. True by default.
            resseq_reindex -- An optional boolean, whether or not to reindex the pdb resseq field. False by default.
            return_text -- An optional boolean, whether or not to return text instead of writing to a file. If True, the filename variable is ignored.

            Returns:
            If return_text is True, a PDB-formatted string. Otherwise, returns nothing.
            
            '''

        if len(self.__parent_molecule.get_atom_information()) > 0: # so the pdb is not empty (if it is empty, don't save)
            
            if serial_reindex == True: self.__parent_molecule.serial_reindex()
            if resseq_reindex == True: self.__parent_molecule.resseq_reindex()
            
            if return_text == False: afile = open(filename,"w")
            else: return_string = ""

            # print out remarks
            for line in self.__parent_molecule.get_remarks():
                remark = "REMARK" + line + "\n"
                
                if return_text == False: afile.write(remark)
                else: return_string = return_string + remark
                
            # print out coordinates
            atom_information = self.__parent_molecule.get_atom_information()
            coordinates = self.__parent_molecule.get_coordinates()
            
            printout = numpy.core.defchararray.add(atom_information['record_name'], numpy.core.defchararray.rjust(atom_information['serial'].astype('|S5'),5))
            printout = numpy.core.defchararray.add(printout, atom_information['name'])
            printout = numpy.core.defchararray.add(printout, atom_information['resname'])
            printout = numpy.core.defchararray.add(printout, atom_information['chainid'])
            printout = numpy.core.defchararray.add(printout, numpy.core.defchararray.rjust(atom_information['resseq'].astype('|S4'),4))
            printout = numpy.core.defchararray.add(printout, '    ')
            printout = numpy.core.defchararray.add(printout,  numpy.core.defchararray.rjust(numpy.array(["%.3f" % t for t in coordinates[:,0]]), 8))
            printout = numpy.core.defchararray.add(printout,  numpy.core.defchararray.rjust(numpy.array(["%.3f" % t for t in coordinates[:,1]]), 8))
            printout = numpy.core.defchararray.add(printout,  numpy.core.defchararray.rjust(numpy.array(["%.3f" % t for t in coordinates[:,2]]), 8))
            printout = numpy.core.defchararray.add(printout,  numpy.core.defchararray.rjust(numpy.array(["%.2f" % t for t in atom_information['occupancy']]), 6))
            printout = numpy.core.defchararray.add(printout,  numpy.core.defchararray.rjust(numpy.array(["%.2f" % t for t in atom_information['tempfactor']]), 6))
            printout = numpy.core.defchararray.add(printout, '          ')
            printout = numpy.core.defchararray.add(printout, atom_information['element'])
            printout = numpy.core.defchararray.add(printout, atom_information['charge'])

            if return_text == False: 
                if printout[0][-1:] == "\n": afile.write("".join(printout) + "\n")
                else: afile.write("\n".join(printout) + "\n")
            else:
                if printout[0][-1:] == "\n": return_string = return_string + ("".join(printout) + "\n")
                else: return_string = return_string + ("\n".join(printout) + "\n")
            
            # print out connect
            if not self.__parent_molecule.get_bonds() is None:
                for index in range(self.__parent_molecule.get_bonds().shape[0]):
                    indices_of_bond_partners = self.__parent_molecule.select_all_atoms_bound_to_selection(numpy.array([index]))
                    if len(indices_of_bond_partners) > 0:

                        if return_text == False: afile.write("CONECT" + str(self.__parent_molecule.get_atom_information()["serial"][index]).rjust(5) + "".join([str(self.__parent_molecule.get_atom_information()["serial"][t]).rjust(5) for t in indices_of_bond_partners]) + "\n")
                        else: return_string = return_string + ("CONECT" + str(self.__parent_molecule.get_atom_information()["serial"][index]).rjust(5) + "".join([str(self.__parent_molecule.get_atom_information()["serial"][t]).rjust(5) for t in indices_of_bond_partners]) + "\n")

            if return_text == False: afile.close()
            else: return return_string
            
        else: print "ERROR: Cannot save a Molecule with no atoms (file name \"" + filename + "\")"

class AtomsAndBonds():
    '''A class for adding and deleting atoms and bonds'''
    
    def __init__(self, parent_molecule_object):
        '''Initializes the pymolecule.AtomsAndBonds class.
                
            Arguments:
            parent_molecule_object -- The pymolecule.Molecule object associated with this class.
            
            '''
        
        self.__parent_molecule = parent_molecule_object

    def create_bonds_by_distance(self, remove_old_bond_data=True, delete_excessive_bonds=True):
        '''Determines which atoms are bound to each other based on their proximity.
                
            Arguments:
            remove_old_bond_data -- An optional boolean, whether or not to discard old bond data before adding in bonds determined by distance. True by default.
            delete_excessive_bonds -- An optional boolean, whether or not to check for and delete excessive bonds. True by default.
            
            '''

        # create/recreate the bond array if needed
        if remove_old_bond_data == True or self.__parent_molecule.get_bonds() is None: self.__parent_molecule.set_bonds(numpy.zeros((len(self.__parent_molecule.information.get_atom_information()), len(self.__parent_molecule.information.get_atom_information()))))
        
        # get the longest bond length on record
        max_bond_length = numpy.max([self.__parent_molecule.get_constants()['bond_length_dict'][key] for key in self.__parent_molecule.get_constants()['bond_length_dict'].keys()])
        
        # which ones could possibly be bound (less than the max_bond_length)
        distances = scipy.spatial.distance.squareform(pdist(self.__parent_molecule.get_coordinates()))
        ones_to_consider = numpy.nonzero(distances < max_bond_length * 1.2)
        
        for index in range(len(ones_to_consider[0])):
            index1 = ones_to_consider[0][index]
            index2 = ones_to_consider[1][index]
            
            if index1 != index2: # so an atom is not bound to itself.__parent_molecule
                key = self.__parent_molecule.get_atom_information()['element_stripped'][index1] +'-' + self.__parent_molecule.get_atom_information()['element_stripped'][index2]
                
                try: bond_dist = self.__parent_molecule.get_constants()['bond_length_dict'][key]
                except:
                    print "ERROR: Unknown bond distance between elements " + self.__parent_molecule.get_atom_information()['element_stripped'][index1] + ' and ' + self.__parent_molecule.get_atom_information()['element_stripped'][index2] + '. Assuming ' + str(max_bond_length) + '.'
                    bond_dist = max_bond_length
                    
                if distances[index1][index2] < bond_dist * 1.2 and distances[index1][index2] > bond_dist * 0.5: # so they should be bonded
                    self.__parent_molecule.add_bond(index1,index2)
                
        if delete_excessive_bonds == True:
            # now do a sanity check. C cannot have more than 4 bonds, O cannot have more than 2, and N cannot have more than 2
            # if more, than use ones closest to ideal bond length
            for index in range(len(self.__parent_molecule.get_atom_information())):
                # get the info of the index atom
                element = self.__parent_molecule.get_atom_information()['element_stripped'][index]
                
                bond_partner_indices = self.__parent_molecule.select_all_atoms_bound_to_selection(numpy.array([index]))
                number_of_bonds = len(bond_partner_indices)
                
                try:
                    if number_of_bonds > self.__parent_molecule.get_constants()['max_number_of_bonds_permitted'][element]: # so this atom has too many bonds
                        # get the distances of this atoms bonds
                        dists = distances[index][bond_partner_indices]
                        
                        # get the ideal distances of those bonds
                        ideal_dists = numpy.empty(len(dists)) # initialize the vector
                        
                        for t in range(len(bond_partner_indices)): # populate the ideal-bond-length vector
                            index_partner = bond_partner_indices[t]
                            element_partner = self.__parent_molecule.get_atom_information()['element_stripped'][index_partner]
                            ideal_dists[t] = self.__parent_molecule.get_constants()['bond_length_dict'][element + '-' + element_partner]
                        
                        diff = numpy.absolute(dists - ideal_dists) # get the distance
                        
                        # identify the bonds to discard
                        indices_in_order = diff.argsort()
                        indicies_to_throw_out = indices_in_order[self.__parent_molecule.get_constants()['max_number_of_bonds_permitted'][element]:]
                        indicies_to_throw_out = bond_partner_indices[indicies_to_throw_out]
                        
                        # discard the extra bonds
                        for throw_out_index in indicies_to_throw_out:
                            self.__parent_molecule.delete_bond(index, throw_out_index)
                        
                except: pass # element probably wasn't in the dictionary

    def get_number_of_bond_partners_of_element(self, atom_index, the_element):
        '''Counts the number of atoms of a given element bonded to a specified atom of interest.
            
            Arguments:
            atom_index -- An int, the index of the atom of interest.
            the_element -- A string describing the element of the neighbors to be counted.
            
            Returns:
            An int, the number of neighboring atoms of the specified element.
            
            '''
        
        # this function is really here for historical reasons. it's similar to the old
        # number_of_neighors_of_element function. it could be done pretty easily with
        # numpy
        
        the_element = the_element.strip()
        bond_partners_selection = self.__parent_molecule.select_all_atoms_bound_to_selection(numpy.array([atom_index]))
        elements = self.__parent_molecule.get_atom_information()['element_stripped'][bond_partners_selection]
        return len(numpy.nonzero(elements == the_element)[0])
        
    def get_index_of_first_bond_partner_of_element(self, atom_index, the_element): 
        '''For a given atom of interest, returns the index of the first neighbor of a specified element.
        
            Arguments:
            atom_index -- An int, the index of the atom of interest.
            the_element -- A string specifying the desired element of the neighbor.
            
            Returns:
            An int, the index of the first neighbor atom of the specified element. If no such neighbor exists, returns -1.
            
            '''

        # this function is really here for historical reasons. it's similar to the old
        # index_of_neighbor_of_element function. it could be done pretty easily with
        # numpy
        
        the_element = the_element.strip()
        bond_partners_selection = self.__parent_molecule.select_all_atoms_bound_to_selection(numpy.array([atom_index]))
        elements = self.__parent_molecule.get_atom_information()['element_stripped'][bond_partners_selection]
        return bond_partners_selection[numpy.nonzero(elements == the_element)[0]][0]

    def delete_bond(self, index1, index2):
        '''Deletes a bond.
        
            Arguments:
            index1 -- An int, the index of the first atom of the bonded pair.
            index2 -- An int, the index of the second atom of the bonded pair.
            
            '''
        bonds = self.__parent_molecule.get_bonds()
        try:
            bonds[index1,index2] = 0
            bonds[index2,index1] = 0
        except: print "Could not delete bond between " + str(index1) + " and " + str(index2) + "."
        
    def add_bond(self, index1, index2, order=1):
        '''Adds a bond.
        
            Arguments:
            index1 -- An int, the index of the first atom of the bonded pair.
            index2 -- An int, the index of the second atom of the bonded pair.
            order -- An optional int, the order of the bond. 1 by default.
            
            '''
        
        bonds = self.__parent_molecule.get_bonds()
        bonds[index1,index2] = order
        bonds[index2,index1] = order

        
    def delete_atom(self, index):
        '''Deletes an atom.
        
            Arguments:
            index -- An int, the index of the atom to delete.
            
            '''

        # remove the atom information
        self.__parent_molecule.set_atom_information(numpy.delete(self.__parent_molecule.get_atom_information(), index))
        
        # remove the coordinates
        self.__parent_molecule.set_coordinates(numpy.delete(self.__parent_molecule.get_coordinates(), index, axis=0))
        try: self.__parent_molecule.set_coordinates_undo_point(numpy.delete(self.__parent_molecule.get_coordinates_undo_point(), index, axis=0))
        except: pass
        
        # remove the relevant bonds
        self.__parent_molecule.set_bonds(numpy.delete(self.__parent_molecule.information.get_bonds().todense(), index, 0))
        self.__parent_molecule.set_bonds(numpy.delete(self.__parent_molecule.information.get_bonds().todense(), index, 1))

        ## The above conversion to and from dense matrices is quite inefficient. This is the beginning of an attempt to fix it with an efficient sparse equivalent.
        #self.__parent_molecule.set_bonds(delete_row_and_col_csr(self.__parent_molecule.information.get_bonds(), index))

        
        # the hierarchy will have to be recomputed
        self.__hierarchy = {}
    
    def add_atom(self, record_name="ATOM", serial=1, name="X", resname="XXX", chainid="X", resseq=1, occupancy=0.0, tempfactor=0.0, charge='', element="X", coordinates=numpy.array([0.0, 0.0, 0.0]), autoindex = True):
        '''Adds an atom.
        
            Arguments:
            record_name -- An optional string, the record name of the atom. "ATOM" is the default.
            serial -- An optional int, the serial field of the atom. 1 is the default.
            name -- An optional string, the name of the atom. "X" is the default.
            resname -- An optional string, the resname of the atom. "XXX" is the default.
            chainid -- An optional string, chainid of the atom. "X" is the default.
            resseq -- An optional int, the resseq field of the atom. 1 is the default.
            occupancy -- An optional float, the occupancy of the atom. 0.0 is the default.
            tempfactor -- An optional float, the tempfactor of the atom. 0.0 is the default.
            charge -- An optional string, the charge of the atom. "" is the default.
            element -- An optional string, the element of the atom. "X" is the default.
            coordinates -- An optional numpy.array, the (x, y, z) coordinates of the atom. numpy.array([0.0, 0.0, 0.0]) is the default.
            
            '''

        # add the atom information
        
        if len(record_name) < 6: record_name = record_name.ljust(6)
        if len(name) < 5:
            if len(name) < 4: name = name.rjust(4) + ' '
            else: name = name.rjust(5)
        if len(resname) < 4: resname = resname.rjust(4)
        if len(chainid) < 2: chainid = chainid.rjust(2)
        if len(charge) < 2: charge = charge.ljust(2)
        if len(element) < 2: element = element.rjust(2)
        
        name_stripped = name.strip()
        resname_stripped = resname.strip()
        chainid_stripped = chainid.strip()
        element_stripped = element.strip()
        
        try: mass = self.__parent_molecule.get_constants()['mass_dict'][element_stripped]
        except: mass = 0.0
        
        # if there is no atom_information, you need to create it.
        if self.__parent_molecule.get_atom_information() is None:
            self.__parent_molecule.information.set_atom_information(numpy.zeros((1,), dtype=[('record_name', '|S6'), ('serial', '<i8'), ('name', '|S5'), ('resname', '|S4'), ('chainid', '|S2'), ('resseq', '<i8'), ('occupancy', '<f8'), ('tempfactor', '<f8'), ('element', '|S2'), ('charge', '|S2'), ('name_stripped', '|S5'), ('resname_stripped', '|S4'), ('chainid_stripped', '|S2'), ('element_stripped', '|S2')]))

        # ********

        atom_information = numpy.ma.resize(self.__parent_molecule.information.get_atom_information(), self.__parent_molecule.information.get_total_number_of_atoms()+1)
        atom_information['record_name'][-1] = record_name
        atom_information['name'][-1] = name
        atom_information['resname'][-1] = resname
        atom_information['chainid'][-1] = chainid
        atom_information['charge'][-1] = charge
        atom_information['element'][-1] = element
        atom_information['name_stripped'][-1] = name_stripped
        atom_information['resname_stripped'][-1] = resname_stripped
        atom_information['chainid_stripped'][-1] = chainid_stripped
        atom_information['element_stripped'][-1] = element_stripped
        atom_information['serial'][-1] = serial
        atom_information['resseq'][-1] = resseq
        atom_information['occupancy'][-1] = occupancy
        atom_information['tempfactor'][-1] = tempfactor
        
        self.__parent_molecule.set_atom_information(atom_information)
        #self.__parent_molecule.assign_masses()
        
        if 'mass' in self.__parent_molecule.information.get_atom_information().dtype.names:
            self.__parent_molecule.information.get_atom_information()['mass'][-1] = mass

        # now add the coordinates
        if self.__parent_molecule.get_coordinates() is None: self.__parent_molecule.set_coordinates(numpy.array([coordinates]))
        else: self.__parent_molecule.set_coordinates(numpy.vstack((self.__parent_molecule.get_coordinates(), coordinates)))
        
        # now add places for bonds, though bonds will only be added if done explicitly, not here
        if self.__parent_molecule.get_bonds() is None:
            self.__parent_molecule.set_bonds(numpy.array([[0]]))
        else:
            self.__parent_molecule.set_bonds(numpy.vstack((self.__parent_molecule.information.get_bonds().todense(), numpy.zeros(self.__parent_molecule.get_total_number_of_atoms()-1))))
            
            self.__parent_molecule.set_bonds(numpy.hstack((
                self.__parent_molecule.get_bonds().todense(),
                numpy.zeros((1,self.__parent_molecule.information.get_total_number_of_atoms())).T
            )))
            
class Selections():
    '''A class for selecting atoms'''

    ######## selections ########
    def __init__(self, parent_molecule_object):
        '''Initializes the pymolecule.Selections class.
                
            Arguments:
            parent_molecule_object -- The pymolecule.Molecule object associated with this class.
            
            '''
        
        self.__parent_molecule = parent_molecule_object

    def select_atoms(self, selection_criteria):
        '''Select a set of atoms based on user-specified criteria.
        
            Arguments:
            selection_criteria -- An dictionary, where the keys correspond to keys in the self.__parent_molecule.information.get_atom_information() structured numpy array, and the values are lists of acceptable matches.
                The selection is a logical "AND" between dictionary entries, but "OR" within the value lists themselves.
                For example: {'atom':['CA','O'], 'chain':'A', 'resname':'PRO'} would select all atoms with the names CA or O that are located in the PRO residues of chain A.
            
            Returns:
            A numpy.array containing the indices of the atoms of the selection.
            
            '''
        
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
            
    def select_atoms_in_bounding_box(self, bounding_box):
        '''Selects all the atoms that are within a bounding box.
        
            Arguments:
            bounding_box -- A 2x3 numpy.array containing the minimum and maximum points of the bounding box. Example: numpy.array([[min_x, min_y, min_z], [max_x, max_y, max_z]]).
            
            Returns:
            A numpy.array containing the indices of the atoms that are within the bounding box.
            
            '''
        
        min_pt = bounding_box[0]
        max_pt = bounding_box[1]
        coordinates = self.__parent_molecule.get_coordinates()
        sel1 = numpy.nonzero((coordinates[:,0] > min_pt[0]))[0]
        sel2 = numpy.nonzero((coordinates[:,0] < max_pt[0]))[0]
        sel3 = numpy.nonzero((coordinates[:,1] > min_pt[1]))[0]
        sel4 = numpy.nonzero((coordinates[:,1] < max_pt[1]))[0]
        sel5 = numpy.nonzero((coordinates[:,2] > min_pt[2]))[0]
        sel6 = numpy.nonzero((coordinates[:,2] < max_pt[2]))[0]
        sel = numpy.intersect1d(sel1, sel2)
        sel = numpy.intersect1d(sel, sel3)
        sel = numpy.intersect1d(sel, sel4)
        sel = numpy.intersect1d(sel, sel5)
        sel = numpy.intersect1d(sel, sel6)
        
        return sel

    def select_all_atoms_bound_to_selection(self, selection):
        '''Selects all the atoms that are bound to a user-specified selection.
        
            Arguments:
            selection -- A numpy.array containing the indices of the user-specified selection.
            
            Returns:
            A numpy.array containing the indices of the atoms that are bound to the user-specified selection. Note that this new selection does not necessarily include the indices of the original user-specified selection.
            
            '''

        if self.__parent_molecule.information.get_bonds() is None:
            print "You need to define the bonds to use select_all_atoms_bound_to_selection()."
            return
        bonds_to_consider = self.__parent_molecule.get_bonds()[selection]
        return numpy.unique(numpy.nonzero(bonds_to_consider)[1])

        

    def select_branch(self, root_atom_index, directionality_atom_index):
        '''Identify an isolated "branch" of a molecular model. Assumes the atoms with indices root_atom_index and directionality_atom_index are bound to one another and that the branch starts at root_atom_index one and "points" in the direction of directionality_atom_index.
            
        Arguments:
        root_atom_index -- An int, the index of the first atom in the branch (the "root").
        directionality_atom_index -- An int, the index of the second atom in the branch, used to establish directionality
        
        Returns:
        A numpy array containing the indices of the atoms of the branch.
        
        '''
        
        # note that this function is mostly retained for legacy reasons. the old version of pymolecule
        # had a branch-identification function. 

        if self.__parent_molecule.get_bonds() is None:
            print "To identify atoms in the same molecule as the atoms of a selection, you need to define the bonds."
            return

        #Make sure atoms are neighboring
        if not directionality_atom_index in self.select_all_atoms_bound_to_selection(numpy.array([root_atom_index])):
            print "The root and directionality atoms, with indices " + str(root_atom_index) + " and " + str(directionality_atom_index) + ", respectively, are not neighboring atoms."
            return

        # first, set up the two indices need to manage the growing list of
        # connected atoms.
        # current_index is the index in the list that you're currently considering
        current_index = 1
            
        # create an "empty" array to store the indices of the connected atoms
        indices_of_this_branch = [root_atom_index, directionality_atom_index] # can't know ahead of time what size, so let's use a python list # -99999 * numpy.ones(len(self.__parent_molecule.information.get_coordinates()), dtype=int) # assume initially that all the atoms belong to this molecule. this list will be shortened, possibly, later if that assumption is incorrect.
            
        while True:
            # get all the neighbors of the current atom
            try: current_atom_index = indices_of_this_branch[current_index]
            except: break # this error because you've reached the end of the larger molecule
            
            neighbors_indices = self.__parent_molecule.select_all_atoms_bound_to_selection(numpy.array([current_atom_index]))
            
            # get the ones in neighbors_indices that are not in indices_of_this_molecule
            new_ones = numpy.setdiff1d(neighbors_indices, indices_of_this_branch)
            indices_of_this_branch.extend(new_ones)
            
            # prepare to look at the next atom in the list
            current_index = current_index + 1

        return numpy.array(indices_of_this_branch)

    def select_atoms_from_same_molecule(self, selection):
        '''Selects all the atoms that belong to the same molecule as a user-defined selection, assuming that the pymolecule.Molecule object actually contains multiple physically distinct molecules that are not bound to each other via covalent bonds.
        
            Arguments:
            selection -- A numpy.array containing the indices of the user-defined selection.
            
            Returns:
            A numpy.array containing the indices of the atoms belonging to the same molecules as the atoms of the user-defined selection.
            
            '''
            
        # If your "Molecule" object actually contains several molecules, this one
        # selects all the atoms from any molecule containing any atom in the selection
        # note that bonds must be defined

        if self.__parent_molecule.get_bonds() is None:
            print "To identify atoms in the same molecule as the atoms of a selection, you need to define the bonds."
            return

        indices = []
        for index in selection:

            # first, set up the two indices need to manage the growing list of
            # connected atoms.
            # current_index is the index in the list that you're currently considering
            current_index = 0
                
            # create an "empty" array to store the indices of the connected atoms
            indices_of_this_molecule = [index] # can't know ahead of time what size, so let's use a python list # -99999 * numpy.ones(len(self.__parent_molecule.information.get_coordinates()), dtype=int) # assume initially that all the atoms belong to this molecule. this list will be shortened, possibly, later if that assumption is incorrect.
                
            while True:
                # get all the neighbors of the current atom
                try: current_atom_index = indices_of_this_molecule[current_index]
                except: break # this error because you've reached the end of the larger molecule
                
                neighbors_indices = self.__parent_molecule.select_all_atoms_bound_to_selection(numpy.array([current_atom_index]))
                
                # get the ones in neighbors_indices that are not in indices_of_this_molecule
                new_ones = numpy.setdiff1d(neighbors_indices, indices_of_this_molecule)
                indices_of_this_molecule.extend(new_ones)
                
                # prepare to look at the next atom in the list
                current_index = current_index + 1
                
            #indices_of_this_molecule = indices_of_this_molecule[:current_index-1] # so the list is prunes down.
            indices.append(indices_of_this_molecule)

        # now merge and remove redundancies
        return numpy.unique(numpy.hstack(indices))
    
    def selections_of_constituent_molecules(self):
        '''Identifies the indices of atoms belonging to separate molecules, assuming that the pymolecule.Molecule object actually contains multiple physically distinct molecules that are not bound to each other via covalent bonds.
        
            Returns:
            A python list of numpy.array objects containing the indices of the atoms belonging to each molecule of the composite pymolecule.Molecule object.
            
            '''

        # If your pymolecule.Molecule object contains multiple molecules (e.g., several
        # chains), this will return a list of selections corresponding to the atoms of
        # each molecule
        
        atoms_not_yet_considered = self.select_all()
        selections = []
        
        while len(atoms_not_yet_considered) > 0:
            # add the atoms in the same molecule as the first atom in atoms_not_yet_considered
            this_molecule_atoms = self.select_atoms_from_same_molecule(numpy.array([atoms_not_yet_considered[0]]))
            
            # now remove these from the atoms_not_yet_considered list
            atoms_not_yet_considered = numpy.setxor1d(this_molecule_atoms, atoms_not_yet_considered, True)
            
            # save the atoms of this molecule
            selections.append(this_molecule_atoms)

        return selections

    def select_atoms_near_other_selection(self, selection, cutoff):
        '''Selects all atoms that are near the atoms of a user-defined selection.
        
            Arguments:
            selection -- A numpy.array containing the indices of the user-defined selection.
            cutoff -- A float, the distance cutoff (in Angstroms).
            
            Returns:
            A numpy.array containing the indices of all atoms near the user-defined selection, not including the atoms of the user-defined selection themselves.
            
            '''

        # note that this does not return a selection that includes the input selection.
        # merge selections as required to get a selection that also includes the input.
        
        invert_selection = self.invert_selection(selection)
        
        selection_coors = self.__parent_molecule.get_coordinates()[selection]
        inversion_coors = self.__parent_molecule.get_coordinates()[invert_selection]
        
        indices_of_nearby = invert_selection[numpy.unique(numpy.nonzero(cdist(inversion_coors, selection_coors) < cutoff)[0])]
        return indices_of_nearby

    def select_atoms_in_same_residue(self, selection):
        '''Selects all atoms that are in the same residue as any of the atoms of a user-defined seleciton. Residues are considered unique if they have a unique combination of resname, resseq, and chainid fields.
        
            Arguments:
            selection -- A numpy.array containing the indices of the user-defined selection.
            
            Returns:
            A numpy.array containing the indices of all atoms in the same residue as any of the atoms of the user-defined selection.
            
            '''

        # get string ids representing the residues of all atoms
        keys = numpy.core.defchararray.add(self.__parent_molecule.get_atom_information()['resname_stripped'], '-')
        keys = numpy.core.defchararray.add(keys, numpy.array([str(t) for t in self.__parent_molecule.get_atom_information()['resseq']]))
        keys = numpy.core.defchararray.add(keys, '-')
        keys = numpy.core.defchararray.add(keys, self.__parent_molecule.get_atom_information()['chainid_stripped'])

        # get the unique keys of the selection
        unique_keys_of_selection = numpy.unique(keys[selection])
        
        # now get all the atoms of these selection keys
        
        # the below works, but is slow for large systems
        #residues = self.__parent_molecule.selections_of_residues()
        #new_selection = numpy.array([], dtype=int)
        #for key in unique_keys_of_selection:
        #    print key
        #    new_selection = numpy.append(new_selection, residues[key])
        
        # let's use this instead, faster for large systems.
        new_selection = numpy.array([], dtype=int)
        for key in unique_keys_of_selection:
            new_selection = numpy.append(new_selection, numpy.nonzero(keys == key)[0])
        
        return new_selection
    
    def invert_selection(self, selection):
        '''Inverts a user-defined selection (i.e., identifies all atoms that are not in the seleciton).
        
            Arguments:
            selection -- A numpy.array containing the indices of the user-defined selection.
            
            Returns:
            A numpy.array containing the indices of all atoms that are not in the user-defined seleciton.
            
            '''

        # selection is a list of atom indices
        all_atoms = numpy.arange(0,len(self.__parent_molecule.get_atom_information()), 1, dtype=int)
        remaining_indicies = numpy.delete(all_atoms, selection)
        return remaining_indicies
    
    def select_all(self):
        '''Selects all the atoms in a pymolecule.Molecule object.
        
            Returns:
            A numpy.array containing the indices of all atoms in the pymolecule.Molecule object.
            
            '''

        return self.select_atoms({})

    def select_close_atoms_from_different_molecules(self, other_mol, cutoff, pairwise_comparison=True, terminate_early=False):
        '''Effectively detects steric clashes between self and another pymolecule.Molecule.
        
            Arguments
            other_mol -- A pymolecule.Molecule object of the other molecule.
            cutoff -- A float, the user-defined distance cutoff in Angstroms.
            pairwise_comparison -- An optional boolean, whether or not to perform a simple pairwise distance comparison (if True) or to use a more sophisitcated method (if False). True by default.
            terminate_early = An optional boolean, whether or not to stop looking for steric clashes once one is found. False by default.
        
            Returns:
            A tuple containing two elements. The first is a numpy.array containing the indices of all nearby atoms from this pymolecule.Molecule object (self). The second is a numpy.array containing the indices of all nearby atoms from the other molecule.
        
            '''
        
        if pairwise_comparison == True:
            
            dists = cdist(self.__parent_molecule.get_coordinates(), other_mol.get_coordinates())
            close_ones = numpy.nonzero(dists < cutoff)
            close_ones_from_mol_parent_molecule = numpy.unique(close_ones[0])
            close_ones_from_mol_other_mol = numpy.unique(close_ones[1])
            
            return (close_ones_from_mol_parent_molecule, close_ones_from_mol_other_mol)
        else: # so do the more complex hierarchical comparison
            # first, do some quick and easy checks
            self_coordinates = self.__parent_molecule.get_coordinates()
            other_coordinates = other_mol.get_coordinates()
            self_hierarchy = self.__parent_molecule.get_hierarchy()
            other_hierarchy = other_mol.get_hierarchy()
            
            margin = numpy.array([cutoff, cutoff, cutoff])
            self_min = numpy.min(self_coordinates,0) - margin
            other_mol_max = numpy.max(other_coordinates,0) + margin
            
            if self_min[0] > other_mol_max[0]: return (numpy.array([]),numpy.array([]))
            if self_min[1] > other_mol_max[1]: return (numpy.array([]),numpy.array([]))
            if self_min[2] > other_mol_max[2]: return (numpy.array([]),numpy.array([]))
            
            self_max = numpy.max(self_coordinates,0) + margin
            other_mol_min = numpy.min(other_coordinates,0) - margin
        
            if other_mol_min[0] > self_max[0]: return (numpy.array([]),numpy.array([]))
            if other_mol_min[1] > self_max[1]: return (numpy.array([]),numpy.array([]))
            if other_mol_min[2] > self_max[2]: return (numpy.array([]),numpy.array([]))
            
            # now assign spheres to the whole molecule, the chains, the residues
            # note that this won't recalculate the data if it's already been calculated
            self.__parent_molecule.define_molecule_chain_residue_spherical_boundaries()
            other_mol.define_molecule_chain_residue_spherical_boundaries()
            
            # if the whole molecules are too far away, give up
            self_cent = self_hierarchy['spheres']['molecule']['center']
            self_rad = self_hierarchy['spheres']['molecule']['radius']
            other_cent = other_hierarchy['spheres']['molecule']['center']
            other_rad = other_hierarchy['spheres']['molecule']['radius']
            mol_dist = numpy.linalg.norm(self_cent - other_cent)
            
            if mol_dist > self_rad + other_rad + cutoff: return (numpy.array([]),numpy.array([])) # the molecules are too far away to clash
    
            # check the chains
            chain_distances = cdist(self_hierarchy['spheres']['chains']['centers'], other_mol.get_hierarchy()['spheres']['chains']['centers'])
            sum1_matrix = numpy.hstack([numpy.array([self_hierarchy['spheres']['chains']['radii']]).T for t in range(len(other_hierarchy['spheres']['chains']['radii']))])
            sum2_matrix = numpy.vstack([numpy.array([other_hierarchy['spheres']['chains']['radii']]) for t in range(len(self_hierarchy['spheres']['chains']['radii']))])
            sum_matrix = sum1_matrix + sum2_matrix + cutoff
            indicies_of_clashing_chains = numpy.nonzero(chain_distances < sum_matrix)
    
            if len(indicies_of_clashing_chains[0]) == 0: return (numpy.array([]),numpy.array([])) # the chains don't clash, so no atoms can either
            
            # check the residues
            residue_distances = cdist(self_hierarchy['spheres']['residues']['centers'], other_mol.get_hierarchy()['spheres']['residues']['centers'])
            sum1_matrix = numpy.hstack([numpy.array([self_hierarchy['spheres']['residues']['radii']]).T for t in range(len(other_hierarchy['spheres']['residues']['radii']))])
            sum2_matrix = numpy.vstack([numpy.array([other_hierarchy['spheres']['residues']['radii']]) for t in range(len(self_hierarchy['spheres']['residues']['radii']))])
            sum_matrix = sum1_matrix + sum2_matrix + cutoff
            
            indicies_of_clashing_residues = numpy.nonzero(residue_distances < sum_matrix)
            
            if len(indicies_of_clashing_residues[0]) == 0: return (numpy.array([]),numpy.array([])) # the residues don't clash, so no atoms can either
            
            # now time to check the atoms
            self_close_atom_indices = numpy.array([], dtype=int)
            other_close_atom_indices = numpy.array([], dtype=int)
            
            for i in range(len(indicies_of_clashing_residues[0])):
                self_res_index = indicies_of_clashing_residues[0][i]
                other_res_index = indicies_of_clashing_residues[1][i]
                
                self_res_name = self_hierarchy['spheres']['residues']['keys'][self_res_index]
                other_res_name = other_hierarchy['spheres']['residues']['keys'][other_res_index]
                
                self_res_indicies = self_hierarchy['residues']['indices'][self_res_name]
                other_res_indicies = other_hierarchy['residues']['indices'][other_res_name]
                
                self_coors = self_coordinates[self_res_indicies]
                other_coors = other_coordinates[other_res_indicies]
                
                some_self_indices, some_other_indices = numpy.nonzero(cdist(self_coors, other_coors) < cutoff)
                if len(some_self_indices) != 0 or len(some_other_indices) != 0: # so there are some
                    self_close_atom_indices = numpy.append(self_close_atom_indices, self_res_indicies[some_self_indices])
                    other_close_atom_indices = numpy.append(other_close_atom_indices, other_res_indicies[some_other_indices])
                    
                    if terminate_early == True: # so don't keep looking once you've found something
                        return (self_close_atom_indices, other_close_atom_indices)

            return (numpy.unique(self_close_atom_indices), numpy.unique(other_close_atom_indices)) # so nothing was found in the end

    def get_molecule_from_selection(self, selection, serial_reindex=True, resseq_reindex=False):
        '''Creates a pymolecule.Molecule from a user-defined atom selection.
        
            Arguments
            selection -- A numpy.array containing the indices of the atoms in the user-defined selection.
            serial_reindex -- An optional boolean, whether or not to reindex the atom serial fields. Default is True.
            resseq_reindex -- An optional boolean, whether or not to reindex the atom resseq fields. Default is False.
        
            Returns:
            A pymolecule.Molecule object containing the atoms of the user-defined selection.
        
            '''

        new_mol = Molecule()
        new_mol.set_coordinates(self.__parent_molecule.get_coordinates()[selection])
        
        # try to get the undo coordinates as well, though they may not have been set
        try: new_mol.information.set_coordinates_undo_point(self.__parent_molecule.get_coordinates_undo_point()[selection])
        except: new_mol.information.set_coordinates_undo_point(None)
        
        new_mol.set_atom_information(self.__parent_molecule.get_atom_information()[selection])
        
        if not self.__parent_molecule.get_bonds() is None:
            new_mol.set_bonds(self.__parent_molecule.get_bonds()[selection])
            new_mol.set_bonds(new_mol.get_bonds()[:,selection])
        else: new_mol.set_bonds(None)
        
        # note that hierarchy will have to be recalculated
        
        if serial_reindex == True: new_mol.information.serial_reindex()
        if resseq_reindex == True: new_mol.information.resseq_reindex()
        return new_mol
    
    def selections_of_chains(self):
        '''Identifies the atom selections of each chain.
        
            Returns:
            A dictionary. The keys of the dictionary correspond to the chainids, and the values are numpy.array objects containing the indices of the associated chain atoms.
        
            '''

        if not 'chains' in self.__parent_molecule.get_hierarchy().keys() : # so it hasn't already been calculated
            unique_chainids = numpy.unique(self.__parent_molecule.get_atom_information()['chainid_stripped'])

            self.__parent_molecule.get_hierarchy()['chains'] = {}
            self.__parent_molecule.get_hierarchy()['chains']['indices'] = {}
            for chainid in unique_chainids:
                self.__parent_molecule.get_hierarchy()['chains']['indices'][chainid] = self.__parent_molecule.select_atoms({'chainid_stripped':chainid})
                
        return self.__parent_molecule.get_hierarchy()['chains']['indices']
    
    def selections_of_residues(self):
        '''Identifies the atom selections of each residue.
        
            Returns:
            A dictionary. The keys of this dictionary correspond to the unique resname-resseq-chainid residue identifiers, and the values are numpy.array objects containing the indices of the associated residue atoms.
        
            '''

        if not 'residues' in self.__parent_molecule.get_hierarchy().keys() : # so it hasn't already been calculated

            keys = numpy.core.defchararray.add(self.__parent_molecule.get_atom_information()['resname_stripped'], '-')
            keys = numpy.core.defchararray.add(keys, numpy.array([str(t) for t in self.__parent_molecule.get_atom_information()['resseq']]))
            keys = numpy.core.defchararray.add(keys, '-')
            keys = numpy.core.defchararray.add(keys, self.__parent_molecule.get_atom_information()['chainid_stripped'])

            unique_resnames = numpy.unique(keys)
            
            self.__parent_molecule.get_hierarchy()['residues'] = {}
            self.__parent_molecule.get_hierarchy()['residues']['indices'] = {}
            for key in unique_resnames:
                resname, resseq, chainid = key.split('-')
                resseq = int(resseq)
                
                self.__parent_molecule.get_hierarchy()['residues']['indices'][key] = self.__parent_molecule.select_atoms({'chainid_stripped':chainid, 'resname_stripped':resname, 'resseq':resseq})
                
        return self.__parent_molecule.get_hierarchy()['residues']['indices']

    def in_same_ring(self, index1, index2):
        '''Determines if two atoms in a Molecule are in the same ring
            
        Arguments:
        index1 -- index of the first atom to see if it is in the ring
        index2 -- index of the second atom to see if it is in the ring
        
        Returns: A set of vertices of a path that has been traversed
        
        '''
        
        if index1 == index2: return True
        
        paths = []
        
        paths = self.__ring_recursive_walk(index1, index2, [], 0)
        
        if len(paths) == 0:
            print "No paths found between two indices"
            return False
        
        #Remove paths that do not start or end at the correct location
        for path in paths:
            if not (index2 in path and index1 in path): paths.remove(path)
            
        #Now need to find the intersection between each combination of paths
        for path1, path2 in itertools.combinations(paths, 2):
            intersection = set(path1).intersection(set(path2))
            
            #Only two elements are in index1 and index2
            if len(intersection) == 2 and index1 in intersection and index2 in intersection:
                return True
            
        return False
    
    def __ring_recursive_walk(self, start, end, already_crossed, ringsize):
        '''Recursive helper method to traverse a graph to search for a circular subgraph
        
        Arguments:
        start -- graph vertex that is currently being traversed
        end -- vertex to be reached
        alreadyCrossed -- list of vertices already processed
        ringsize - counts number of vertices traversed
        
        Returns: A set of vertices of a path that has been traversed
        
        '''
        
        paths = []
        
        ring_size += 1
        already_crossed.append(start)
        
        #Base case 1: Second point is reached
        if start == end:
            paths.append(already_crossed)
            return paths
        
        #Base case 2: Max ring size is reached
        if ring_size >= self.__max_ring_size():
            paths.append(already_crossed)
            return paths
        
        #Base case 3: No new neighbors
        #Get a list of all the atoms that atom:index is connected to that haven't been previously evaluated
        neighbors = self.self.__parent_molecule.select_all_atoms_bound_to_selection(numpy.array([start]))[:]
        
        for neighbor_index in neighbors:
            if neighbor_index in already_crossed:
                neighbors.remove(neighbor_index)
                
        if len(neighbors) == 0:
            paths.append(already_crossed)
            return paths
        
        for neighbor in neighbors:
            paths.extend(self.__ring_recursive_walk(neighbor, end, already_crossed[:], ring_size))
            
        return paths
        

class Manipulation():
    '''A class for translating and rotating the atomic coordinates of a pymolecule.Molecule object'''
    
    def __init__(self, parent_molecule_object):
        '''Initializes the pymolecule.Manipulation class.
                
            Arguments:
            parent_molecule_object -- The pymolecule.Molecule object associated with this class.
            
            '''
        
        self.__parent_molecule = parent_molecule_object

    def set_coordinate_undo_point(self):
        '''Sets ("saves") the undo point of the atom coordinates. Any subsequent manipulations of atomic coordinates can be "undone" by reseting to this configuration via the coordinate_undo function.'''

        self.__parent_molecule.set_coordinates_undo_point(self.__parent_molecule.get_coordinates().copy())

    def coordinate_undo(self):
        '''Resets the coordinates of all atoms to those saved using the set_coordinate_undo_point function.'''
        
        self.__parent_molecule.set_coordinates(self.__parent_molecule.get_coordinates_undo_point().copy())
    
    def set_atom_location(self, atom_index, new_location):
        '''Translates the entire molecular model (without rotating) so that the atom with the specified index is located at the specified coordinate.
                
            Arguments:
            atom_index -- An int, the index of the target atom.
            new_location -- A numpy.array specifying the new (x, y, z) coordinate of the specified atom.
            
            Returns: A numpy.array specifying the (delta_x, delta_y, delta_z) vector by which the pmolecule.Molecule was translated.
            
            '''
        
        if new_location.shape == (3,): new_location = numpy.array([new_location])
        
        currentloc = self.__parent_molecule.get_coordinates()[atom_index]
        delta = new_location - currentloc
        
        self.translate_molecule(delta)

        return delta
    
    def translate_molecule(self, delta):
        '''Translate all the atoms of the molecular model by a specified vector.
            
        Arguments:
        delta -- A numpy.array (delta_x, delta_y, delta_z) specifying the amount to move each atom along the x, y, and z coordinates.
        
        '''

        if delta.shape == (3,): delta = numpy.array([delta])
        
        self.__parent_molecule.set_coordinates(self.__parent_molecule.get_coordinates() + delta)

        if 'spheres' in self.__parent_molecule.get_hierarchy().keys():
            # so update location of hierarchical elements
            self.__parent_molecule.get_hierarchy()['spheres']['molecule']['center'] = self.__parent_molecule.get_hierarchy()['spheres']['molecule']['center'] + delta
            self.__parent_molecule.get_hierarchy()['spheres']['chains']['centers'] = self.__parent_molecule.get_hierarchy()['spheres']['chains']['centers'] + delta
            self.__parent_molecule.get_hierarchy()['spheres']['residues']['centers'] = self.__parent_molecule.get_hierarchy()['spheres']['residues']['centers'] + delta

    def rotate_molecule_around_a_line_between_points(self, line_point1, line_point2, rotate):
        '''rotate the molecular model about a line segment. The end points of the line segment are explicitly specified coordinates.
            
            Arguments:
            line_point1 -- A numpy.array (x, y, z) corresponding to one end of the line segment.
            line_point2 -- A numpy.array (x, y, z) corresponding to the other end of the line segment.
            rotate -- A float, the angle of rotation, in radians.
            
            '''

        if line_point1.shape == (1,3): line_point1 = line_point1[0]
        if line_point2.shape == (1,3): line_point2 = line_point2[0]

        a = line_point1[0]
        b = line_point1[1]
        c = line_point1[2]
        #d = line_point2[0]
        #e = line_point2[1]
        #f = line_point2[2]
        
        delta = line_point2 - line_point1
        
        u = delta[0] #d-a
        v = delta[1] #e-b
        w = delta[2] #f-c
        
        v_2_plus_w_2 = numpy.power(v, 2) + numpy.power(w, 2)
        u_2_plus_w_2 = numpy.power(u, 2) + numpy.power(w, 2)
        u_2_plus_v_2 = numpy.power(u, 2) + numpy.power(v, 2)
        u_2_plus_v_2_plus_w_2 = u_2_plus_v_2 + numpy.power(w, 2)

        cos = numpy.cos(rotate)
        sin = numpy.sin(rotate)
        
        ux_plus_vy_plus_wz = numpy.sum(self.__parent_molecule.get_coordinates() * delta,1)

        # Now rotate molecule. In a perform world, I'd have an awesome
        # better numpified version of this, perhaps with tensor or matrix multiplication
        
        for t in range(len(self.__parent_molecule.get_coordinates())): # so t is an atom index
            x_not, y_not, z_not = self.__parent_molecule.get_coordinates()[t]
            ux_plus_vy_plus_wz = u*x_not + v*y_not + w*z_not
            
            self.__parent_molecule.get_coordinates()[t][0] = (a*v_2_plus_w_2 + u*(-b*v-c*w+ux_plus_vy_plus_wz)+(-a*v_2_plus_w_2+u*(b*v+c*w-v*y_not-w*z_not)+v_2_plus_w_2*x_not)*cos+numpy.sqrt(u_2_plus_v_2_plus_w_2)*(-c*v+b*w-w*y_not+v*z_not)*sin)#/u_2_plus_v_2_plus_w_2
            self.__parent_molecule.get_coordinates()[t][1] = (b*u_2_plus_w_2 + v*(-a*u-c*w+ux_plus_vy_plus_wz)+(-b*u_2_plus_w_2+v*(a*u+c*w-u*x_not-w*z_not)+u_2_plus_w_2*y_not)*cos+numpy.sqrt(u_2_plus_v_2_plus_w_2)*(c*u-a*w+w*x_not-u*z_not)*sin)#/u_2_plus_v_2_plus_w_2
            self.__parent_molecule.get_coordinates()[t][2] = (c*u_2_plus_v_2 + w*(-a*u-b*v+ux_plus_vy_plus_wz)+(-c*u_2_plus_v_2+w*(a*u+b*v-u*x_not-v*y_not)+u_2_plus_v_2*z_not)*cos+numpy.sqrt(u_2_plus_v_2_plus_w_2)*(-b*u+a*v-v*x_not+u*y_not)*sin)#/u_2_plus_v_2_plus_w_2
            
        self.__parent_molecule.set_coordinates(self.__parent_molecule.get_coordinates() * (1.0/u_2_plus_v_2_plus_w_2))

        # here I'm going to just delete the hierarchical info because I'm lazy.
        try: del self.__parent_molecule.get_hierarchy()['spheres'] # calculated bounding spheres, if any, are no longer valid.
        except: pass

    def rotate_molecule_around_a_line_between_atoms(self, line_point1_index, line_point2_index, rotate):
        '''Rotate the molecular model about a line segment. The end points of the line segment are atoms of specified indices.
            
            Arguments:
            line_point1_index -- An int, the index of the first atom at one end of the line segment.
            line_point2_index -- An int, the index of the second atom at the other end of the line segment.
            rotate -- A float, the angle of rotation, in radians.
            
            '''
        
        pt1 = self.__parent_molecule.get_coordinates()[line_point1_index]
        pt2 = self.__parent_molecule.get_coordinates()[line_point2_index]
        self.rotate_molecule_around_a_line_between_points(pt1, pt2, rotate)

        try: del self.__parent_molecule.get_hierarchy()['spheres'] # calculated bounding spheres, if any, are no longer valid.
        except: pass
      
    def rotate_molecule_around_pivot_point(self, pivot, thetax, thetay, thetaz):
        '''Rotate the molecular model around a specified atom.
            
            Arguments:
            pivot -- A numpy.array, the (x, y, z) coordinate about which the molecular model will be rotated.
            thetax -- A float, the angle to rotate relative to the x axis, in radians.
            thetay -- A float, the angle to rotate relative to the y axis, in radians.
            thetaz -- A float, the angle to rotate relative to the z axis, in radians.
            
            '''
        
        if pivot.shape == (3,): pivot = numpy.array([pivot])
        
        # First, move the Molecule so the pivot is at the origin
        self.__parent_molecule.set_coordinates(self.__parent_molecule.get_coordinates() - pivot)
        
        # do the rotation
        sinx = numpy.sin(thetax)
        siny = numpy.sin(thetay)
        sinz = numpy.sin(thetaz)
        cosx = numpy.cos(thetax)
        cosy = numpy.cos(thetay)
        cosz = numpy.cos(thetaz)
        
        rot_matrix = numpy.array([[(cosy * cosz), (sinx * siny * cosz + cosx * sinz), (sinx * sinz - cosx * siny * cosz)], [-(cosy * sinz), (cosx * cosz - sinx * siny * sinz), (cosx * siny * sinz + sinx * cosz)], [siny, -(sinx * cosy), (cosx * cosy)]])
        self.__parent_molecule.set_coordinates(numpy.dot(rot_matrix, self.__parent_molecule.get_coordinates().T).T)
        
        # now move the pivot point back to it's old location
        self.__parent_molecule.set_coordinates(self.__parent_molecule.get_coordinates() + pivot)

        try: del self.__parent_molecule.information.hierarchy['spheres'] # calculated bounding spheres, if any, are no longer valid.
        except: pass

    def rotate_molecule_around_pivot_atom(self, pivot_index, thetax, thetay, thetaz):
        '''Rotate the molecular model around a specified atom.
            
            Arguments:
            pivot_index -- An int, the index of the atom about which the molecular model will be rotated.
            thetax -- A float, the angle to rotate relative to the x axis, in radians.
            thetay -- A float, the angle to rotate relative to the y axis, in radians.
            thetaz -- A float, the angle to rotate relative to the z axis, in radians.
            
            '''
        
        pivot = self.__parent_molecule.get_coordinates()[pivot_index]
        self.rotate_molecule_around_pivot_point(pivot, thetax, thetay, thetaz)

        try: del self.__parent_molecule.get_hierarchy()['spheres'] # calculated bounding spheres, if any, are no longer valid.
        except: pass

class Geometry():
    '''A class containing a few gemoetry functions. Note that numpy should be used for most geometry functions.'''

    def __init__(self, parent_molecule_object):
        '''Initializes the pymolecule.Geometry class.
                
            Arguments:
            parent_molecule_object -- The pymolecule.Molecule object associated with this class.
            
            '''
        
        self.__parent_molecule = parent_molecule_object

    def get_angle_between_three_points(self, pt1, pt2, pt3):
        '''Computes the angle (in radians) formed by three points (numpy.array objects).
            
            Arguments
            pt1 -- A numpy.array (x, y, z) representing the first of the three 3D points.
            pt2 -- A numpy.array (x, y, z) representing the second of the three 3D points.
            pt3 -- A numpy.array (x, y, z) representing the third of the three 3D points.
            
            Returns:
            A float containing the angle between the three points, in radians
            
            '''
        
        vector1 = pt1 - pt2
        vector2 = pt3 - pt2
        
        vector1_mag = numpy.linalg.norm(vector1)
        vector2_mag = numpy.linalg.norm(vector2)
        
        #Make sure vectors aren't <0,0,0>
        if vector1_mag < 1e-10 or vector2_mag < 1e-10:
            print "One of vectors to determine angle is < 0, 0, 0 >...returning 0."
            return 0
        
        vector1 = vector1 / vector1_mag
        vector2 = vector2 / vector2_mag
        dot_prod = numpy.dot(vector1, vector2)
        
        #Prevent errors that can rarely occur
        if dot_prod > 1.0: dot_prod = 1.0
        if dot_prod < -1.0: dot_prod = -1.0
        
        return numpy.arccos(dot_prod)

    def get_dihedral_angle(self, pt1, pt2, pt3, pt4):
        '''Calculates the dihedral angle formed by four points (numpy.array objects).
                
            Arguments:
            pt1 -- A numpy.array (x, y, z) representing the first 3D point.
            pt2 -- A numpy.array (x, y, z) representing the second 3D point.
            pt3 -- A numpy.array (x, y, z) representing the third 3D point.
            pt4 -- A numpy.array (x, y, z) representing the fourth 3D point.
            
            Returns:
            A float containing the dihedral angle between the four points, in radians.
            
            '''
        
        b1 = pt2 - pt1
        b2 = pt3 - pt2
        b3 = pt4 - pt3
        
        b2Xb3 = numpy.cross(b2, b3)
        b1Xb2 = numpy.cross(b1, b2)
        
        b1XMagb2 = numpy.linalg.norm(b2) * b1
        
        return numpy.arctan2(numpy.dot(b1XMagb2, b2Xb3), numpy.dot(b1Xb2, b2Xb3))

    def is_planar(self, pt1, pt2, pt3, pt4, planarity_cutoff=0.3):
        '''Checks whether four points (numpy.array) lie in a common plane.
            
            Arguments:
            pt1 -- A numpy.array (x, y, z) representing a 3D point.
            pt2 -- A numpy.array (x, y, z) representing a 3D point.
            pt3 -- A numpy.array (x, y, z) representing a 3D point.
            pt4 -- A numpy.array (x, y, z) representing a 3D point.
            planarity_cutoff -- An optional float. How much the points can deviate (in Angstroms) and still be considered planar. The default is 0.2.
                
            Returns:
            A boolean, whether the 4 points can be considered planar.
                
            '''
        planarity_deviation = self.get_planarity_deviation(pt1, pt2, pt3, pt4)
        return (planarity_deviation < planarity_cutoff)
            
    def get_planarity_deviation(self, pt1, pt2, pt3, pt4):
        '''Determines how close four points (numpy.array objects) come to lying in a common plane.
                
            Arguments:
            pt1 -- A numpy.array (x, y, z) representing a 3D point.
            pt2 -- A numpy.array (x, y, z) representing a 3D point.
            pt3 -- A numpy.array (x, y, z) representing a 3D point.
            pt4 -- A numpy.array (x, y, z) representing a 3D point.
                
            Returns:
            A float, the minimum distance between one point and the plane formed by the other three.
                
            '''
        
        # note that minimal efforts were made to "numpify" this section. It's mostly legacy code.
        
        x1 = pt1[0]
        y1 = pt1[1]
        z1 = pt1[2]
        x2 = pt2[0]
        y2 = pt2[1]
        z2 = pt2[2]
        x3 = pt3[0]
        y3 = pt3[1]
        z3 = pt3[2]
        x4 = pt4[0]
        y4 = pt4[1]
        z4 = pt4[2]
        
        A = (y1*(z2-z3))+(y2*(z3-z1))+(y3*(z1-z2))
        B = (z1*(x2-x3))+(z2*(x3-x1))+(z3*(x1-x2))
        C = (x1*(y2-y3))+(x2*(y3-y1))+(x3*(y1-y2))
        D = ((-x1)*((y2*z3)-(y3*z2)))+((-x2)*((y3*z1)-(y1*z3)))+((-x3)*((y1*z2)-(y2*z1)))
        denom = numpy.sqrt(numpy.power(A,2) + numpy.power(B,2) + numpy.power(C,2))
        if denom == 0: return 0 # implies straight line
        distance1=numpy.fabs((A*x4)+(B*y4)+(C*z4)+D)/denom
        
        A1 = (y1*(z2-z4))+(y2*(z4-z1))+(y4*(z1-z2))
        B1 = (z1*(x2-x4))+(z2*(x4-x1))+(z4*(x1-x2))
        C1 = (x1*(y2-y4))+(x2*(y4-y1))+(x4*(y1-y2))
        D1 = ((-x1)*((y2*z4)-(y4*z2)))+((-x2)*((y4*z1)-(y1*z4)))+((-x4)*((y1*z2)-(y2*z1)))
        distance2=(numpy.fabs((A1*x3)+(B1*y3)+(C1*z3)+D1))/(numpy.sqrt(numpy.power(A1,2) + numpy.power(B1,2) + numpy.power(C1,2)))
        
        A2 = (y1*(z4-z3))+(y4*(z3-z1))+(y3*(z1-z4))
        B2 = (z1*(x4-x3))+(z4*(x3-x1))+(z3*(x1-x4))
        C2 = (x1*(y4-y3))+(x4*(y3-y1))+(x3*(y1-y4))
        D2 = ((-x1)*((y4*z3)-(y3*z4)))+((-x4)*((y3*z1)-(y1*z3)))+((-x3)*((y1*z4)-(y4*z1)))
        distance3=(numpy.fabs((A2*x2)+(B2*y2)+(C2*z2)+D2))/(numpy.sqrt(numpy.power(A2,2) + numpy.power(B2,2) + numpy.power(C2,2)))
        
        A3 = (y4*(z2-z3))+(y2*(z3-z4))+(y3*(z4-z2))
        B3 = (z4*(x2-x3))+(z2*(x3-x4))+(z3*(x4-x2))
        C3 = (x4*(y2-y3))+(x2*(y3-y4))+(x3*(y4-y2))
        D3 = ((-x4)*((y2*z3)-(y3*z2)))+((-x2)*((y3*z4)-(y4*z3)))+((-x3)*((y4*z2)-(y2*z4)))
        distance4=(numpy.fabs((A3*x1)+(B3*y1)+(C3*z1)+D3))/(numpy.sqrt(numpy.power(A3,2) + numpy.power(B3,2) + numpy.power(C3,2)))
        
        return numpy.min(numpy.array([distance1, distance2, distance3, distance4]))

class OtherMolecules():
    '''A class for characterizing the relationships between multiple pymolecule.'''

    def __init__(self, parent_molecule_object):
        '''Initializes the pymolecule.OtherMolecules class.
                
            Arguments:
            parent_molecule_object -- The pymolecule.Molecule object associated with this class.
            
            '''
        
        self.__parent_molecule = parent_molecule_object

    def get_other_molecule_aligned_to_this(self, other_mol, tethers, weight_mat = None):
        '''Aligns a molecule to self (this pymolecule.Molecule object) using a quaternion RMSD alignment.
                
            Arguments:
            other_mol -- A pymolecule.Molecule that is to be aligned to this one.
            tethers -- A tuple of two numpy.array objects, where each array contains the indices of self
                and other_mol, respectively, such that equivalent atoms are listed in the same order.
                So, for example, if (atom 1, self = atom 3, other) and (atom2, self = atom6, other)
                than the tethers would be (numpy.array([1,2]), numpy.array([3,6])).
            
            '''

        # Adapted from Itzhack Y. Bar-Itzhack. New Method for Extracting the Quaternion from a Rotation Matrix. Journal of Guidance, Control, and Dynamics 2000
        if tethers is None: raise Exception('No tethers specified for RMSD alignment')
        elif tethers.shape[0] != 2: raise Exception('Tethers should have only 2 rows')
        
        #If weight_matrix isn't specified, then treat all atoms equally
        if weight_mat is None: weight_mat = numpy.identity(tethers.shape[1])

        # get the atoms corresponding to the tethers, in tether order
        self_static_atom_coordinates = self.__parent_molecule.get_coordinates()[tethers[0]]
        other_dynamic_atom_coordinates = other_mol.get_coordinates()[tethers[1]]
        
        # translate the tether atoms to the origin
        center_self = numpy.mean(self_static_atom_coordinates, 0)
        center_other = numpy.mean(other_dynamic_atom_coordinates, 0)
        
        self_static_atom_coordinates = self_static_atom_coordinates - center_self
        other_dynamic_atom_coordinates = other_dynamic_atom_coordinates - center_other
        
        # get optimal rotation
        M = numpy.dot(numpy.dot(numpy.transpose(self_static_atom_coordinates), weight_mat), other_dynamic_atom_coordinates)

        #Create symmetric 4x4 matrix K from M
        K = numpy.array([[M[0,0] + M[1,1] + M[2,2], M[1,2] - M[2,1], M[2,0] - M[0,2], M[0,1] - M[1,0]],
                         [M[1,2] - M[2,1], M[0,0] - M[1,1] - M[2,2], M[1,0] + M[0,1], M[2,0] + M[0,2]],
                         [M[2,0] - M[0,2], M[1,0] + M[0,1], M[1,1] - M[0,0] - M[2,2], M[1,2] + M[2,1]],
                         [M[0,1] - M[1,0], M[2,0] + M[0,2], M[1,2] + M[2,1], M[2,2] - M [0,0] - M[1,1]]])
    
        #Find eigenvector associated with the most positive eigenvalue of K.  Multiple quaternions can
        E,V = numpy.linalg.eig(K)
        index = numpy.argmax(E)
        eigenvector = V[:,index]
        rot_quat = Quaternion(eigenvector[0], eigenvector[1], eigenvector[2], eigenvector[3])
        
        rot_mat = rot_quat.to_matrix()
    
        #Apply translation and rotation to the other molecule
        new_mol = other_mol.copy()
        
        new_mol.set_coordinates(new_mol.information.get_coordinates() - center_other)
        new_mol.set_coordinates(numpy.dot(new_mol.information.get_coordinates(),rot_mat))
        new_mol.set_coordinates(new_mol.information.get_coordinates() + center_self)

        return new_mol

    def steric_clash_with_another_molecule(self, other_mol, cutoff, pairwise_comparison=True):
        '''Detects steric clashes between the pymolecule.Molecule (self) and another pymolecule.Molecule.
        
            Arguments
            other_mol -- The pymolecule.Molecule object that will be evaluated for steric clashes.
            cutoff -- A float, the user-defined distance cutoff in Angstroms.
            pairwise_comparison -- An optional boolean, whether or not to perform a simple pairwise distance comparison (if True) or to use a more sophisitcated method (if False). True by default.
        
            Returns:
            A boolean.  True if steric clashes are present, False if they are not
        
            '''
        
        if pairwise_comparison == True: # so use a simple pairwise comparison to find close atoms
            indices1, indices2 = self.__parent_molecule.select_close_atoms_from_different_molecules(other_mol, cutoff, True)
        else: # so the more sophisticated heirarchical method
            indices1, indices2 = self.__parent_molecule.select_close_atoms_from_different_molecules(other_mol, cutoff, False, True) # terminate early is true because you don't want all close ones

        if len(indices1) == 0 and len(indices2) == 0: return False
        else: return True

    def merge_with_another_molecule(self, other_molecule):
        '''Merges two molecular models into a single model.
                
            Arguments:
            other_molecule -- A molecular model (pymolecule.Molecule object).
            
            Returns:
            A single pymolecule.Molecule object containing the atoms of this model combined with the atoms of other_molecule.
            
            '''

        merged = self.__parent_molecule.copy()
        
        # if masses have been assigned to either molecule, they must be assigned to both
        if 'mass' in merged.information.get_atom_information().dtype.names or 'mass' in self.__parent_molecule.get_atom_information().dtype.names:
            self.__parent_molecule.assign_masses()
            merged.information.assign_masses()
            
        merged.filename = ""
        merged.get_remarks().extend(other_molecule.get_remarks())
        merged.set_atom_information(numpy.lib.recfunctions.stack_arrays((merged.get_atom_information(), other_molecule.get_atom_information()),usemask=False))
        
        merged.set_coordinates(numpy.vstack((merged.get_coordinates(), other_molecule.get_coordinates())))
        
        merged.set_coordinates_undo_point(None) 
        
        # merge the bonds, though note that bonds between the two molecules will not be set
        if not merged.get_bonds() is None and not other_molecule.get_bonds() is None:
            bonds1 = merged.get_bonds().copy()
            bonds2 = other_molecule.get_bonds().copy()
            
            bonds1_v2 = numpy.hstack((bonds1.todense(), numpy.zeros((bonds1.shape[0], bonds2.shape[0]))))
            bonds2_v2 = numpy.hstack((numpy.zeros((bonds2.shape[0], bonds1.shape[0])), bonds2.todense()))
    
            merged.set_bonds(numpy.vstack((bonds1_v2, bonds2_v2)))
        else: merged.set_bonds(None)
        
        # the molecule center will be redefined, so you might as well start the hierarchy all over
        try: del merged.information.hierarchy['spheres']
        except: pass

        return merged
    
    def get_distance_to_another_molecule(self, other_molecule, pairwise_comparison=True):
        '''Computes the minimum distance between any of the atoms of this molecular model and any of the atoms of a second specified model.
            
            Arguments:
            other_molecule -- a pymolecule.Molecule, the other molecular model.
            pairwise_comparison -- An optional boolean, whether or not to perform a simple pairwise distance comparison (if True) or to use a more sophisitcated method (if False). True by default.

            Returns:
            A float, the minimum distance between any two atoms of the two specified molecular models (self and other_molecule).
            
            '''
        
        if pairwise_comparison == True: 
            return numpy.amin(cdist(self.__parent_molecule.get_coordinates(), other_molecule.get_coordinates()))
        else: # so use the more sophisticated methods for comparison
            # note that this is not the fastest way to do this, but it uses existing functions
            # and is still pretty fast, so I'm going to stick with it.
            
            # first, get a cutoff distance. Let's just do a quick survey of the two molecules to pick a good one.
            self_tmp = self.__parent_molecule.get_coordinates()[numpy.arange(0,len(self.__parent_molecule.get_coordinates()), len(self.__parent_molecule.get_coordinates())/10.0, dtype=int)]
            other_tmp = other_molecule.get_coordinates()[numpy.arange(0,len(other_molecule.get_coordinates()), len(other_molecule.get_coordinates())/10.0, dtype=int)]
            cutoff = numpy.amin(cdist(self_tmp, other_tmp))
            
            # now get all the indices that come within that cutoff
            self_indices, other_indices = self.__parent_molecule.select_close_atoms_from_different_molecules(other_molecule, cutoff, False)
            
            self_coors = self.__parent_molecule.get_coordinates()[self_indices]
            self_other = other_molecule.get_coordinates()[other_indices]
            
            return numpy.amin(cdist(self_coors, self_other))

    def get_rmsd_equivalent_atoms_specified(self, other_mol, tethers):
        '''Calculates the RMSD between this pymolecule.Molecle object and another, where equivalent atoms are explicitly specified.
            
            Arguments:
            other_mol -- The other pymolecule.Molecule object.
            tethers -- A tuple of two numpy.array objects, where each array contains the indices of self
                and other_mol, respectively, such that equivalent atoms are listed in the same order.
                So, for example, if (atom 1, self = atom 3, other) and (atom2, self = atom6, other)
                than the tethers would be (numpy.array([1,2]), numpy.array([3,6])).
                
            Returns:
            A float, the RMSD between self and other_mol.
                
            '''

        if len(self.__parent_molecule.get_coordinates()) != len(other_mol.get_coordinates()):
            print "Cannot calculate RMSD: number of atoms are not equal."
            print "\t" + str(len(self.__parent_molecule.get_coordinates())) + " vs. " + str(len(other_mol.get_coordinates())) + " atoms."
            return 99999999.0
        
        self_coors_in_order = self.__parent_molecule.get_coordinates()[tethers[0]]
        other_coors_in_order = other_mol.get_coordinates()[tethers[1]]
        
        delta = self_coors_in_order - other_coors_in_order
        norm_squared = numpy.sum(delta**2,axis=-1)
        rmsd = numpy.power(numpy.sum(norm_squared) / len(norm_squared), 0.5)
        return rmsd
    
    def get_rmsd_order_dependent(self, other_mol):
        '''Calculates the RMSD between two structures, where equivalent atoms are listed in the same order.
            
            Arguments:
            other_mol -- The other pymolecule.Molecule object.
                
            Returns:
            A float, the RMSD between self and other_mol.
                
            '''
        
        self_index_in_order = numpy.arange(0,len(self.__parent_molecule.get_coordinates()),1,dtype=int)
        other_index_in_order = numpy.arange(0,len(other_mol.get_coordinates()),1,dtype=int)
        
        return self.get_rmsd_equivalent_atoms_specified(other_mol, (self_index_in_order, other_index_in_order))

class Quaternion:
    '''A class supporting quaternion arithmetic'''
    
    def __init__(self, s, x, y, z):
        '''Initializes the pymolecule.Quaternion class.
                
            Arguments:
            s -- ????
            x -- ????
            y -- ????
            z -- ????
            
            '''
        
        self.v = numpy.array([s,x,y,z])
    
    def __str__(self):
        '''String containing quaternion information in the form of s x y z
            
            Returns:
            A string, containing all information about this quaternion
            
            '''
        
        return "" + str(self.v[0]) + "\t" + str(self.v[1]) + "\t" + str(self.v[2]) + "\t" + str(self.v[3])
    
    def copy(self):
        '''Returns a copy of self'''
        
        return Quaternion(self.v[0], self.v[1], self.v[2], self.v[3])
    
    def load_from_mat(self, m):
        '''Converts a rotation matrix that is pure orthogonal (det(matrix)=1) into a Quaternion. Adapted from http://www.euclideanspace.com/maths/geometry/rotations/conversions/matrixToQuaternion/index.htm 
                
            Arguments:
            m -- A 2D numpy.array representing a pure orthogonal matrix
                
            '''
        #Make sure m is a 3x3 array
        if m.shape[0] != 3 or m.shape[1] != 3:
            print "Could not load quaternion from matrix...size is not (3x3)"
            return
        
        #Check that matrix is orthogonal. m_T = m_inv
        if not numpy.array_equal(numpy.transpose(m),numpy.linalg.inv(m)):
            print "Load Quaternion error. Matrix is not orthogonal"
            return
        
        #Need to make sure that the matrix is special orthogonal
        if numpy.fabs(1-numpy.linalg.det(m)) > 0.000001: #Done for rounding errors
            print "Load Quaternion error.  Determinant is not 1"
            return
        
        #First calculate the sum of the diagonal elements
        t = m.trace()
        
        if t > 0:
            S = numpy.sqrt(t + 1.0) * 2
            self.v[0] = .25 * S
            self.v[1] = (m[2,1] - m[1,2]) / S
            self.v[2] = (m[0,2] - m[2,0]) / S
            self.v[3] = (m[1,0] - m[0,1]) / S
        elif m[0,0] > m[1,1] and m[0,0] > m[2,2]:
            S = numpy.sqrt(1.0 + m[0,0] - m[1,1] - m[2,2]) * 2
            self.v[0] = (m[2,1] - m[1,2]) / S
            self.v[1] = .25 * S
            self.v[2] = (m[0,1] + m[1,0]) / S
            self.v[3] = (m[0,2] + m[2,0]) / S
        elif m[1,1] > m[2,2]:
            S = numpy.sqrt(1.0 + m[1,1] - m[0,0] - m[2,2]) * 2
            self.v[0] = (m[0,2] - m[2,0]) / S
            self.v[1] = (m[0,1] + m[1,0]) / S
            self.v[2] = .25 * S
            self.v[3] = (m[2,1] + m[1,2]) / S
        else:
            S = numpy.sqrt(1.0) * 2
            self.v[0] = (m[1,0] - m[0,1]) / S
            self.v[1] = (m[0,2] + m[2,0]) / S
            self.v[2] = (m[2,1] + m[1,2]) / S
            self.v[3] = .25 * S
    
    def rep_as_44_matrix(self):
        '''Creates a 4x4 matrix representation of the Quaternion.
            
            Returns:
            A 4x4 numpy array
            
            '''
        
        n = self.normalize()
        qw = n.v[0]
        qx = n.v[1]
        qy = n.v[2]
        qz = n.v[3]
        
        return numpy.array([[qw,qx,qy,qz],[-qx,qw,-qz,qy],[-qy,qz,qw,-qx],[-qz,-qy,qx,qw]])

    def to_matrix(self):
        '''Converts to a normalized 3x3 matrix
            
            Returns:
            A 3x3 numpy.array, corresponding to the quaternion
            
            '''
        
        #First normalize
        n = self.normalize()
        qw = n.v[0]
        qx = n.v[1]
        qy = n.v[2]
        qz = n.v[3]
        return numpy.array(
                           [[1.0 - 2.0 * qy * qy - 2.0 * qz * qz, 2.0 * qx * qy - 2.0 * qz * qw, 2.0 * qx * qz + 2.0 * qy * qw],
                            [2.0 * qx * qy + 2.0 * qz * qw, 1.0 - 2.0 * qx * qx - 2.0 * qz * qz, 2.0 * qy * qz - 2.0 * qx * qw],
                            [2.0 * qx * qz - 2.0 * qy * qw,2.0 * qy * qz + 2.0 * qx * qw, 1.0 - 2.0 * qy * qy - 2.0 * qx * qx]])
                
    def add(self, q2):
        '''Adds two quaternions
            
        Arguments:
        q2 -- A quaternion, to be added to self
        
        Returns:
        A Quaternion, with the values corresponding to self + q2
        
        '''
        
        return Quaternion(self.v[0] + q2.v[0], self.v[1] + q2.v[1], self.v[2] + q2.v[2], self.v[3] + q2.v[3])
    
    def invert(self):
        '''Takes the inverse of the quaternion for "division"
            
            Returns:
            A Quaternion, with the values corresponding to self^-1
            
            '''
        
        return Quaternion(self.v[0], -1 * self.v[1], -1 * self.v[2], -1 * self.v[3])

    def minus(self, q2):
        '''Multiplies two quaternions
        
        Arguments:
        q2 -- A quaternion, to be subtracted from self
        
        Returns:
        A Quaternion, with the values corresponding to self - q2
        
        '''
        
        return Quaternion(self.v[0] - q2.v[0], self.v[1] - q2.v[1], self.v[2] - q2.v[2], self.v[3] - q2.v[3])
    
    def multiply(self, q2):
        '''Multiplies two quaternions
            
            Arguments:
            q2 -- A quaternion, to be multiplied with self
            
            Returns:
            A Quaternion, with the values corresponding to self * q2
            
            '''
        
        return Quaternion(self.v[0] * q2.v[0] - self.v[1] * q2.v[1] - self.v[2] * q2.v[2] - self.v[3] * q2.v[3],
                          self.v[1] * q2.v[0] + self.v[0] * q2.v[1] + self.v[2] * q2.v[3] - self.v[3] * q2.v[2],
                          self.v[0] * q2.v[2] - self.v[1] * q2.v[3] + self.v[2] * q2.v[0] + self.v[3] * q2.v[1],
                          self.v[0] * q2.v[3] + self.v[1] * q2.v[2] - self.v[2] * q2.v[1] + self.v[3] * q2.v[0])
    
    def normalize(self):
        '''Normalizes the quaternion
            
            Returns:
            A normalized Quaternion
            
            '''
        
        #First normalize
        n = numpy.sqrt(numpy.power(self.v[0],2) + numpy.power(self.v[1],2) + numpy.power(self.v[2],2) + numpy.power(self.v[3],2))
        
        return Quaternion(self.v[0] / n, self.v[1] / n, self.v[2] / n, self.v[3] / n)
                
    def scale(self, scalar):
        '''Scales a quaternion
                
            Arguments:
            scalar -- the value to scale the quaternion by
            
            Returns:
            A Quaternion, with the values corresponding to self * scalar
                
            '''
    
        return Quaternion(self.v[0] * scalar, self.v[1] * scalar, self.v[2] * scalar, self.v[3] * scalar)

# here's the actual Molecule class
class Molecule:
    '''Loads, saves, and manupulates molecuar models. The main pymolecule class.'''
    
    def __init__ (self):
        '''Initializes the variables of the Molecule class.'''
        
        self.fileio = FileIO(self)
        self.atoms_and_bonds = AtomsAndBonds(self)
        self.selections = Selections(self)
        self.manipulation = Manipulation(self)
        self.information = Information(self)
        self.other_molecule = OtherMolecules(self)
        self.geometry = Geometry(self)
  
    # Information methods
    def get_coordinates(self): return self.information.get_coordinates()
    def get_filename(self): return self.information.get_filename()
    def get_remarks(self): return self.information.get_remarks()
    def get_atom_information(self): return self.information.get_atom_information()
    def get_coordinates_undo_point(self): return self.information.get_coordinates_undo_point()
    def get_bonds(self): return self.information.get_bonds()
    def get_hierarchy(self): return self.information.get_hierarchy()
    def get_constants(self): return self.information.get_constants()
    def get_center_of_mass(self, selection=None): return self.information.get_center_of_mass(selection)
    def get_geometric_center(self, selection=None): return self.information.get_geometric_center(selection)
    def get_total_mass(self, selection=None): return self.information.get_total_mass(selection)
    def get_total_number_of_atoms(self, selection=None): return self.information.get_total_number_of_atoms(selection)
    def get_total_number_of_heavy_atoms(self, selection=None): return self.information.get_total_number_of_heavy_atoms(selection)
    def get_bounding_box(self, selection=None, padding=0.0): return self.information.get_bounding_box(selection, padding)
    def get_bounding_sphere(self, selection=None, padding=0.0): return self.information.get_bounding_sphere(selection, padding)
    
    def set_filename(self,filename): self.information.set_filename(filename)
    def set_remarks(self,remarks): self.information.set_remarks(remarks)
    def set_atom_information(self,atom_information): self.information.set_atom_information(atom_information)
    def set_coordinates(self,coordinates): self.information.set_coordinates(coordinates)
    def set_coordinates_undo_point(self,coordinates_undo_point): self.information.set_coordinates_undo_point(coordinates_undo_point)
    def set_bonds(self,bonds): self.information.set_bonds(bonds)
    def set_hierarchy(self,hierarchy): self.information.set_hierarchy(hierarchy)
    
    def assign_masses(self): self.information.assign_masses()
    def assign_elements_from_atom_names(self, selection=None): self.information.assign_elements_from_atom_names(selection)
    def define_molecule_chain_residue_spherical_boundaries(self): self.information.define_molecule_chain_residue_spherical_boundaries()
    def serial_reindex(self): self.information.serial_reindex()
    def resseq_reindex(self): self.information.resseq_reindex()
    
    # File I/O class methods
    def load_pym_into(filename): self.fileio.load_pym_into(filename)
    def load_pdb_into(self, filename, bonds_by_distance=True, serial_reindex=True, resseq_reindex=False): self.fileio.load_pdb_into(filename, bonds_by_distance, serial_reindex, resseq_reindex)
    def load_pdb_into_using_file_object(self, file_obj, bonds_by_distance=True, serial_reindex=True, resseq_reindex=False): self.fileio.load_pdb_into_using_file_object(file_obj, bonds_by_distance, serial_reindex, resseq_reindex)
    def save_pym(self, filename, save_bonds=False, save_filename=False, save_remarks=False, save_hierarchy=False, save_coordinates_undo_point=False): self.fileio.save_pym(filenmae, save_bonds, save_filename, save_remarks, save_hierarchy, save_coordinates_undo_point)
    def save_pdb(self, filename="", serial_reindex=True, resseq_reindex=False, return_text = False): self.fileio.save_pdb(filename, serial_reindex, resseq_reindex, return_text)
    
    # Atoms and Bonds class methods
    def get_number_of_bond_partners_of_element(self, atom_index, the_element): return self.atoms_and_bonds.get_number_of_bond_partners_of_element(atom_index, the_element)
    def get_index_of_first_bond_partner_of_element(self, atom_index, the_element): return self.atoms_and_bonds.get_index_of_first_bond_partner_of_element(atom_index, the_element)
    
    def create_bonds_by_distance(self, remove_old_bond_data=True, delete_excessive_bonds=True): self.atoms_and_bonds.create_bonds_by_distance(remove_old_bond_data, delete_excessive_bonds)
    def delete_bond(self, index1, index2): self.atoms_and_bonds.delete_bond(index1, index2)
    def add_bond(self, index1, index2, order = 1): self.atoms_and_bonds.add_bond(index1, index2, order)
    def delete_atom(self, index): self.atoms_and_bonds.delete_atom(index)
    def add_atom(self, record_name="ATOM", serial=1, name="X", resname="XXX", chainid="X", resseq=1, occupancy=0.0, tempfactor=0.0, charge='', element="X", coordinates=numpy.array([0.0, 0.0, 0.0]), autoindex = True): self.atoms_and_bonds.add_atom(record_name, serial, name, resname, chainid, resseq, occupancy, tempfactor, charge, element, coordinates, autoindex)
    
    # Selections class
    def get_molecule_from_selection(self, selection, serial_reindex=True, resseq_reindex=False): return self.selections.get_molecule_from_selection(selection, serial_reindex, resseq_reindex)

    def select_atoms(self, selection_criteria): return self.selections.select_atoms(selection_criteria)
    def select_atoms_in_bounding_box(self, bounding_box): return self.selections.select_atoms_in_bounding_box(bounding_box)
    def select_branch(self, root_atom_index, directionality_atom_index): return self.selections.select_branch(root_atom_index, directionality_atom_index)
    def select_all_atoms_bound_to_selection(self, selections): return self.selections.select_all_atoms_bound_to_selection(selections)
    def select_atoms_from_same_molecule(self, selection): return self.selections.select_atoms_from_same_molecule(root_atom_index)
    def selections_of_constituent_molecules(self): return self.selections.selections_of_constituent_molecules()
    def select_atoms_near_other_selection(self, selection, cutoff): return self.selections.select_atoms_near_other_selection(selection, cutoff)
    def select_atoms_in_same_residue(self, selection): return self.selections.select_atoms_in_same_residue(selection)
    def invert_selection(self, selection): return self.selections.invert_selection(selection)
    def select_all(self): return self.selections.select_all()
    def select_close_atoms_from_different_molecules(self, other_mol, cutoff, pairwise_comparison=True, terminate_early=False): return self.selections.select_close_atoms_from_different_molecules(other_mol, cutoff, pairwise_comparison, terminate_early)
    def selections_of_chains(self): return self.selections.selections_of_chains()
    def selections_of_residues(self): return self.selections.selections_of_residues()
    
    # Manipulation class
    def set_atom_location(self, atom_index, new_location): return self.manipulation.set_atom_location(atom_index, new_location)

    def coordinate_undo(self): self.manipulation.coordinate_undo()
    def translate_molecule(self, delta): self.manipulation.translate_molecule(delta)
    def rotate_molecule_around_a_line_between_points(self, line_point1, line_point2, rotate): self.manipulation.rotate_molecule_around_a_line_between_points(line_point1, line_point2, rotate)
    def rotate_molecule_around_a_line_between_atoms(self, line_point1_index, line_point2_index, rotate): self.manipulation.rotate_molecule_around_a_line_between_atoms(line_point1_index, line_point2_index, rotate)
    def rotate_molecule_around_pivot_point(self, pivot, thetax, thetay, thetaz): self.manipulation.rotate_molecule_around_pivot_point(pivot, thetax, thetay, thetaz)
    def rotate_molecule_around_pivot_atom(self, pivot_index, thetax, thetay, thetaz): self.manipulation.rotate_molecule_around_pivot_atom(pivot_index, thetax, thetay, thetaz)
        
    # Geometry class
    def get_angle_between_three_points(self, pt1, pt2, pt3): return self.geometry.get_angle_between_three_points(pt1, pt2, pt3)
    def get_dihedral_angle(self, pt1, pt2, pt3, pt4): return self.geometry.get_dihedral_angle(pt1, pt2, pt3, pt4)
    def get_planarity_deviation(self, pt1, pt2, pt3, pt4): return self.geometry.get_planarity_deviation(pt1, pt2, pt3, pt4)
    def is_planar(self, pt1, pt2, pt3, pt4, planarity_cutoff=0.2): return self.geometry.is_planar(pt1, pt2, pt3, pt4, planarity_cutoff)
    
    # Other molecule class
    def get_other_molecule_aligned_to_this(self, other_mol, tethers): return self.other_molecule.get_other_molecule_aligned_to_this(other_mol, tethers) #Add Weight Matrix
    def get_distance_to_another_molecule(self, other_molecule, pairwise_comparison=True): return self.other_molecule.get_distance_to_another_molecule(other_moelcule, pairwise_comparison)
    def get_rmsd_equivalent_atoms_specified(self, other_mol, tethers): return self.other_molecule.get_rmsd_equivalent_atoms_specified(other_mol, tethers)
    def get_rmsd_order_dependent(self, other_mol): return self.other_molecule.get_rmsd_order_dependent(other_mol)
    def steric_clash_with_another_molecule(self, other_mol, cutoff, pairwise_comparison=True): return self.other_molecule.steric_clash_with_another_molecule(other_mol, cutoff, pairwise_comparison)
    def merge_with_another_molecule(self, other_molecule): return self.other_molecule.merge_with_another_molecule(other_molecule)
        
    ######## supporting functions ########
    
    def numpy_structured_array_remove_field(self, narray, field_names): # surprised this doesn't come with numpy
        '''Removes a specific field name from a structured numpy array.
                
            Arguments:
            narray -- A structured numpy array.
            field_names -- A list of strings, where each string is one of the field names of narray.
            
            Returns:
            A structured numpy array identical to narray, but with the field names in field_names removed.
            
            '''
        
        names = list(narray.dtype.names) # now remove the coordinates from the atom_information object to save memory
        for f in field_names: names.remove(f)
        return narray[names]
    
    def delete_row_and_col_csr(mat, i):
        ### This code is necessary for deleting atoms now that we keep bond info in a csr_matrix
        ### CURRENTLY UNTESTED AND NOT USED IN CODE!!!!!
        if not isinstance(mat, scipy.sparse.csr_matrix):
            raise ValueError("works only for CSR format -- use .tocsr() first")
        n = mat.indptr[i+1] - mat.indptr[i]
        if n > 0:
            mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
            mat.data = mat.data[:-n]
            mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
            mat.indices = mat.indices[:-n]
        mat.indptr[i:-1] = mat.indptr[i+1:]
        mat.indptr[i:] -= n
        mat.indptr = mat.indptr[:-1]
        mat._shape = (mat._shape[0]-1, mat._shape[1])
        
        
        mat = mat.transpose()
        
        n = mat.indptr[i+1] - mat.indptr[i]
        if n > 0:
            mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
            mat.data = mat.data[:-n]
            mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
            mat.indices = mat.indices[:-n]
        mat.indptr[i:-1] = mat.indptr[i+1:]
        mat.indptr[i:] -= n
        mat.indptr = mat.indptr[:-1]
        mat._shape = (mat._shape[0]-1, mat._shape[1])
        mat = mat.transpose()
        return mat
        
    
    
    def __is_number(self, s):
        
        '''Determines whether or not a string represents a number.
                
            Arguments:
            s -- A string (e.g., "5.4").
            
            Returns:
            A boolean, whether or not the string can be represented by a float.
            
            '''

        try:
            float(s)
            return True
        except ValueError:
            return False
        
    def copy(self):
        '''Returns an exact copy (pymolecule.Molecule) of this Molecule object.  Undo points are NOT copied.
            
            Returns:
            A pymolecule.Molecule, containing to the same atomic information as this pymolecule.Molecule object.
                
            '''
        
        new_molecule = Molecule()
        new_molecule.set_filename(self.get_filename()[:])
        new_molecule.set_remarks(self.get_remarks()[:])
        new_molecule.set_atom_information(self.get_atom_information().copy())
        new_molecule.set_coordinates(self.get_coordinates().copy())
        
        if not self.get_bonds() is None: new_molecule.set_bonds(self.get_bonds().copy())
        else: new_molecule.set_bonds(None)
        
        new_molecule.set_hierarchy(copy.deepcopy(self.get_hierarchy()))
        
        return new_molecule
