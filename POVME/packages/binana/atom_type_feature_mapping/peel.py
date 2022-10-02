# BINANA is released under the GNU General Public License (see http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't hesitate to contact me, Jacob
# Durrant, at jdurrant [at] ucsd [dot] edu. If you use BINANA in your work, please cite
# Durrant, J. D. and J. A. McCammon (2011). "BINANA: A novel algorithm for ligand-binding
# characterization." J Mol Graph Model 29(6): 888-893.

import math
import os
import sys
import re
import textwrap
#import toroidPoints
import POVME.packages.pymolecule.pymolecule as pymolecule
import numpy
import gzip
import copy
import time
import scipy.ndimage.interpolation as sni
from scipy import append

defaultParams={}
defaultParams['default_gaussian_variance'] = 2.0
#defaultParams['close_contacts_dist1_cutoff'] = 2.5
#defaultParams['close_contacts_dist2_cutoff'] = 4.0
#defaultParams['electrostatic_dist_cutoff'] = 4.0
#defaultParams['active_site_flexibility_dist_cutoff'] = 4.0
#defaultParams['hydrophobic_dist_cutoff'] = 4.0
#defaultParams['hydrophilic_dist_cutoff'] = 4.0
defaultParams['hydroph_gaussian_variance'] = 2.0
defaultParams['hydroph_gaussian_cutoff'] = 7.0
defaultParams['hydrogen_bond_dist_cutoff'] = 4.0
defaultParams['hydrogen_bond_acceptor_radius'] = 2.0
defaultParams['hydrogen_bond_acceptor_gaussian_variance'] = 1.0
defaultParams['hydrogen_bond_donor_gaussian_mean'] = 2.85
defaultParams['hydrogen_bond_donor_gaussian_variance'] = 0.3
defaultParams['hydrogen_bond_angle_cutoff'] = 40.0
#defaultParams['pi_padding_dist'] = 0.75
#defaultParams['pi_pi_interacting_dist_cutoff'] = 4.0
defaultParams['aromatic_height_cutoff'] = 6.0
defaultParams['aromatic_gaussian_mean'] = 2.0
defaultParams['aromatic_gaussian_variance'] = 0.7 #This is made up
#defaultParams['pi_stacking_angle_tolerance'] = 30.0
#defaultParams['T_stacking_angle_tolerance'] = 30.0
#defaultParams['T_stacking_closest_dist_cutoff'] = 5.0
#defaultParams['cation_pi_dist_cutoff'] = 6.0
#defaultParams['salt_bridge_dist_cutoff'] = 5.5
defaultParams['receptor'] = ''
#defaultParams['ligand'] = ''
defaultParams['output_dir'] = ''
defaultParams['output_file'] = ''



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




class point:
    #x=99999.0
    #y=99999.0
    #z=99999.0

    #def __init__ (self, x, y ,z):
    #    self.x = float(x)
    #    self.y = float(y)
    #    self.z = float(z)

    def __repr__(self):
        return 'point(%f, %f, %f)' %(self.x, self.y, self.z)

    def __init__(self, coords):
        self.x = coords[0]
        self.y = coords[1]
        self.z = coords[2]

    def __getitem__(self, index):
        if index == 0:
            return self.x
        elif index == 1:
            return self.y
        elif index == 2:
            return self.z
        else:
            raise Exception('Index %i out of range - Points only have 3 members corresponding to their coordinates' %(index))

    def __setitem__(self, index, value):
        if index == 0:
            self.x = value
        elif index == 1:
            self.y = value
        elif index == 2:
            self.z = value
        else:
            raise Exception('Index %i out of range - Points only have 3 members corresponding to their coordinates' %(index))

    def __tuple__(self):
        return (self.x, self.y, self.z)

    def __list__(self):
        return [self.x, self.y, self.z]

    def copy_of(self):
        return point([self.x, self.y, self.z])

    def print_coords(self):
        print(str(self.x)+"\t"+str(self.y)+"\t"+str(self.z))

    def coords(self):
        return [self.x,self.y,self.z]

    def scalar_mult_new(self,c):
        return point([self.x*c, self.y*c, self.z*c])

    def scalar_mult_inplace(self,c):
        self.x*=c
        self.y*=c
        self.z*=c

    def point_sum_new(self, otherPoint):
        return point([self.x+otherPoint.x, self.y+otherPoint.y, self.z+otherPoint.z])

    def snap(self,reso): # snap the point to a grid
        self.x = round(self.x / reso) * reso
        self.y = round(self.y / reso) * reso
        self.z = round(self.z / reso) * reso

    def dist_to(self,apoint):
        return math.sqrt(math.pow(self.x - apoint.x,2) + math.pow(self.y - apoint.y,2) + math.pow(self.z - apoint.z,2))

    def vector_to_new(self, otherPoint):
        delX = otherPoint[0]-self.x
        delY = otherPoint[1]-self.y
        delZ = otherPoint[2]-self.z
        return point([delX, delY, delZ])

    def vector_sub_new(self, otherPoint):
        delX = self.x-otherPoint[0]
        delY = self.y-otherPoint[1]
        delZ = self.z-otherPoint[2]
        return point([delX, delY, delZ])


    def vector_to_inplace(self, otherPoint):
        self.x = otherPoint[0]-self.x
        self.y = otherPoint[1]-self.y
        self.z = otherPoint[2]-self.z

    def vector_sub_inplace(self, otherPoint):
        self.x = self.x-otherPoint[0]
        self.y = self.y-otherPoint[1]
        self.z = self.z-otherPoint[2]

    def description(self):
        return str(self.x) + " " + str(self.y) + " " + str(self.z)

    def magnitude(self):
        return self.dist_to(point([0,0,0]))

    def CreatePDBLine(self, index):

        output = "ATOM "
        output = output + str(index).rjust(6) + "X".rjust(5) + "XXX".rjust(4)
        output = output + ("%.3f" % self.x).rjust(18)
        output = output + ("%.3f" % self.y).rjust(8)
        output = output + ("%.3f" % self.z).rjust(8)
        output = output + "X".rjust(24)
        return output

class atom:

    def __init__ (self):
        self.atomname = ""
        self.residue = ""
        self.coordinates = point([99999, 99999, 99999])
        self.element = ""
        self.PDBIndex = ""
        self.line=""
        self.atomtype=""
        self.IndeciesOfAtomsConnecting=[]
        self.charge = 0
        self.resid = 0
        self.chain = ""
        self.structure = ""
        self.comment = ""

    def copy_of(self):
        theatom = atom()
        theatom.atomname = self.atomname
        theatom.residue = self.residue
        theatom.coordinates = self.coordinates.copy_of()
        theatom.element = self.element
        theatom.PDBIndex = self.PDBIndex
        theatom.line= self.line
        theatom.atomtype= self.atomtype
        theatom.IndeciesOfAtomsConnecting = self.IndeciesOfAtomsConnecting[:]
        theatom.charge = self.charge
        theatom.resid = self.resid
        theatom.chain = self.chain
        theatom.structure = self.structure
        theatom.comment = self.comment

        return theatom

    def string_id(self):
        toreturn = ""
        if self.chain.strip() != '': toreturn = toreturn + self.chain.strip() + ":"
        toreturn = toreturn + self.residue.strip() + '(' + str(self.resid) + "):" + self.atomname.strip() + '(' + str(self.PDBIndex) + ')'
        return toreturn

    def CreatePDBLine(self, index):

        output = "ATOM "
        output = output + str(index).rjust(6) + self.atomname.rjust(5) + self.residue.rjust(4)
        output = output + ("%.3f" % self.coordinates.x).rjust(18)
        output = output + ("%.3f" % self.coordinates.y).rjust(8)
        output = output + ("%.3f" % self.coordinates.z).rjust(8)
        output = output + self.element.rjust(24)
        return output

    def NumberOfNeighbors(self):
        return len(self.IndeciesOfAtomsConnecting)

    def AddNeighborAtomIndex(self, index):
        if not (index in self.IndeciesOfAtomsConnecting):
            self.IndeciesOfAtomsConnecting.append(index)

    def SideChainOrBackBone(self): # only really applies to proteins, assuming standard atom names
        if self.atomname.strip() == "CA" or self.atomname.strip() == "C" or self.atomname.strip() == "O" or self.atomname.strip() == "N":
            return "BACKBONE"
        else:
            return "SIDECHAIN"

    def ReadPDBLine(self, Line):
        self.line = Line
        self.atomname = Line[11:16].strip()

        if len(self.atomname)==1:
            self.atomname = self.atomname + "  "
        elif len(self.atomname)==2:
            self.atomname = self.atomname + " "
        elif len(self.atomname)==3:
            self.atomname = self.atomname + " " # This line is necessary for babel to work, though many PDBs in the PDB would have this line commented out

        self.coordinates = point([float(Line[30:38]), float(Line[38:46]), float(Line[46:54])])

        # now atom type (for pdbqt)
        self.atomtype = Line[77:79].strip().upper()

        if Line[69:76].strip() != "":
            self.charge = float(Line[69:76])
        else:
            self.charge = 0.0

        if self.element == "": # try to guess at element from name
            two_letters = self.atomname[0:2].strip().upper()
            if two_letters=='BR':
                self.element='BR'
            elif two_letters=='CL':
                self.element='CL'
            elif two_letters=='BI':
                self.element='BI'
            elif two_letters=='AS':
                self.element='AS'
            elif two_letters=='AG':
                self.element='AG'
            elif two_letters=='LI':
                self.element='LI'
            #elif two_letters=='HG':
            #    self.element='HG'
            elif two_letters=='MG':
                self.element='MG'
            elif two_letters=='MN':
                self.element='MN'
            elif two_letters=='RH':
                self.element='RH'
            elif two_letters=='ZN':
                self.element='ZN'
            elif two_letters=='FE':
                self.element='FE'
            else: #So, just assume it's the first letter.
                # Any number needs to be removed from the element name
                self.element = self.atomname
                self.element = self.element.replace('0','')
                self.element = self.element.replace('1','')
                self.element = self.element.replace('2','')
                self.element = self.element.replace('3','')
                self.element = self.element.replace('4','')
                self.element = self.element.replace('5','')
                self.element = self.element.replace('6','')
                self.element = self.element.replace('7','')
                self.element = self.element.replace('8','')
                self.element = self.element.replace('9','')
                self.element = self.element.replace('@','')

                self.element = self.element[0:1].strip().upper()

        self.PDBIndex = Line[6:12].strip()
        self.residue = Line[16:20]
        self.residue = " " + self.residue[-3:] # this only uses the rightmost three characters, essentially removing unique rotamer identification

        try:
            self.resid = int(Line[23:26])
        except ValueError:
            pass

        self.chain = Line[21:22]
        if self.residue.strip() == "": self.residue = " MOL"

                

class PDB:
    
    def __init__ (self):
        self.AllAtoms={}
        self.NonProteinAtoms = {}
        self.max_x = -9999.99
        self.min_x = 9999.99
        self.max_y = -9999.99
        self.min_y = 9999.99
        self.max_z = -9999.99
        self.min_z = 9999.99
        self.rotateable_bonds_count = 0
        self.functions = MathFunctions()
        self.protein_resnames = ["ALA", "ARG", "ASN", "ASP", "ASH", "ASX", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH", "GLX", "GLY", "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
        self.aromatic_rings = []
        self.charges = [] # a list of points

    def LoadPDB(self, FileName, min_x=-9999.99, max_x=9999.99, min_y=-9999.99, max_y=9999.99, min_z=-9999.99, max_z=9999.99):

        autoindex = 1

        self.__init__()

        # Now load the file into a list
        myfile = open(FileName,"r")
        lines = myfile.readlines()
        myfile.close()

        atom_already_loaded = [] # going to keep track of atomname_resid_chain pairs, to make sure redundants aren't loaded. This basically
                                 # gets rid of rotomers, I think.

        for t in range(0,len(lines)):
            line=lines[t]

            if "between atoms" in line and " A " in line:
                    self.rotateable_bonds_count = self.rotateable_bonds_count + 1

            if len(line) >= 7:
                if line[0:4]=="ATOM" or line[0:6]=="HETATM": # Load atom data (coordinates, etc.)
                    TempAtom = atom()
                    TempAtom.ReadPDBLine(line)

                    if TempAtom.coordinates.x > min_x and TempAtom.coordinates.x < max_x and TempAtom.coordinates.y > min_y and TempAtom.coordinates.y < max_y and TempAtom.coordinates.z > min_z and TempAtom.coordinates.z < max_z:

                        if self.max_x < TempAtom.coordinates.x: self.max_x = TempAtom.coordinates.x
                        if self.max_y < TempAtom.coordinates.y: self.max_y = TempAtom.coordinates.y
                        if self.max_z < TempAtom.coordinates.z: self.max_z = TempAtom.coordinates.z

                        if self.min_x > TempAtom.coordinates.x: self.min_x = TempAtom.coordinates.x
                        if self.min_y > TempAtom.coordinates.y: self.min_y = TempAtom.coordinates.y
                        if self.min_z > TempAtom.coordinates.z: self.min_z = TempAtom.coordinates.z

                        key = TempAtom.atomname.strip() + "_" + str(TempAtom.resid) + "_" + TempAtom.residue.strip() + "_" + TempAtom.chain.strip() # this string unique identifies each atom

                        if key in atom_already_loaded and TempAtom.residue.strip() in self.protein_resnames: # so this is a protein atom that has already been loaded once
                            self.printout("Warning: Duplicate protein atom detected: \"" + TempAtom.line.strip() + "\". Not loading this duplicate.")
                            print("")

                        if not key in atom_already_loaded or not TempAtom.residue.strip() in self.protein_resnames: # so either the atom hasn't been loaded, or else it's a non-protein atom
                                                                                                            # so note that non-protein atoms can have redundant names, but protein atoms cannot.
                                                                                                            # This is because protein residues often contain rotamers
                            atom_already_loaded.append(key) # so each atom can only be loaded once. No rotamers.
                            self.AllAtoms[autoindex] = TempAtom # So you're actually reindexing everything here.
                            if not TempAtom.residue[-3:] in self.protein_resnames: self.NonProteinAtoms[autoindex] = TempAtom

                            autoindex = autoindex + 1

        self.CheckProteinFormat()

        self.CreateBondsByDistance() # only for the ligand, because bonds can be inferred based on atomnames from PDB
        self.assign_aromatic_rings()
        self.assign_charges()

    def printout(self, thestring):
        lines = textwrap.wrap(thestring, 80)
        for line in lines:
            print(line)

    def SavePDB(self, filename):
        f = open(filename, 'w')
        towrite = self.SavePDBString()
        if towrite.strip() == "": towrite = "ATOM      1  X   XXX             0.000   0.000   0.000                       X" # just so no PDB is empty, VMD will load them all
        f.write(towrite)
        f.close()

    def SavePDBString(self):

        ToOutput = ""

        # write coordinates
        for atomindex in self.AllAtoms:
            ToOutput = ToOutput + self.AllAtoms[atomindex].CreatePDBLine(atomindex) + "\n"

        return ToOutput

    def AddNewAtom(self, atom):

        # first get available index
        t = 1
        while t in list(self.AllAtoms.keys()):
            t = t + 1

        # now add atom
        self.AllAtoms[t] = atom

    def SetResname(self, resname):
        for atomindex in self.AllAtoms:
            self.AllAtoms[atomindex].residue = resname

    def connected_atoms_of_given_element(self, index, connected_atom_element):
        atom = self.AllAtoms[index]
        connected_atoms = []
        for index2 in atom.IndeciesOfAtomsConnecting:
            atom2 = self.AllAtoms[index2]
            if atom2.element == connected_atom_element:
                connected_atoms.append(index2)
        return connected_atoms

    def connected_heavy_atoms(self, index):
        atom = self.AllAtoms[index]
        connected_atoms = []
        for index2 in atom.IndeciesOfAtomsConnecting:
            atom2 = self.AllAtoms[index2]
            if atom2.element != "H": connected_atoms.append(index2)
        return connected_atoms

    def CheckProteinFormat(self):
        curr_res = ""
        first = True
        residue = []

        for atom_index in self.AllAtoms:
            atom = self.AllAtoms[atom_index]

            key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain

            if first == True:
                curr_res = key
                first = False

            if key != curr_res:

                self.CheckProteinFormat_process_residue(residue, last_key)

                residue = []
                curr_res = key

            residue.append(atom.atomname.strip())
            last_key = key

        self.CheckProteinFormat_process_residue(residue, last_key)


    def CheckProteinFormat_process_residue(self, residue, last_key):
        temp = last_key.strip().split("_")
        resname = temp[0]
        real_resname = resname[-3:]
        resid = temp[1]
        chain = temp[2]

        if real_resname in self.protein_resnames: # so it's a protein residue

            if not "N" in residue:
                self.printout('Warning: There is no atom named "N" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine secondary structure. If this residue is far from the active site, this warning may not affect the NNScore.')
                print("")
            if not "C" in residue:
                self.printout('Warning: There is no atom named "C" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine secondary structure. If this residue is far from the active site, this warning may not affect the NNScore.')
                print("")
            if not "CA" in residue:
                self.printout('Warning: There is no atom named "CA" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine secondary structure. If this residue is far from the active site, this warning may not affect the NNScore.')
                print("")

            if real_resname == "GLU" or real_resname == "GLH" or real_resname == "GLX":
                if not "OE1" in residue:
                    self.printout('Warning: There is no atom named "OE1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "OE2" in residue:
                    self.printout('Warning: There is no atom named "OE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")

            if real_resname == "ASP" or real_resname == "ASH" or real_resname == "ASX":
                if not "OD1" in residue:
                    self.printout('Warning: There is no atom named "OD1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "OD2" in residue:
                    self.printout('Warning: There is no atom named "OD2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")

            if real_resname == "LYS" or real_resname == "LYN":
                if not "NZ" in residue:
                    self.printout('Warning: There is no atom named "NZ" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")

            if real_resname == "ARG":
                if not "NH1" in residue:
                    self.printout('Warning: There is no atom named "NH1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "NH2" in residue:
                    self.printout('Warning: There is no atom named "NH2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")

            if real_resname == "HIS" or real_resname == "HID" or real_resname == "HIE" or real_resname == "HIP":
                if not "NE2" in residue:
                    self.printout('Warning: There is no atom named "NE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "ND1" in residue:
                    self.printout('Warning: There is no atom named "ND1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-cation and salt-bridge interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")

            if real_resname == "PHE":
                if not "CG" in residue:
                    self.printout('Warning: There is no atom named "CG" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CD1" in residue:
                    self.printout('Warning: There is no atom named "CD1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CD2" in residue:
                    self.printout('Warning: There is no atom named "CD2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CE1" in residue:
                    self.printout('Warning: There is no atom named "CE1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CE2" in residue:
                    self.printout('Warning: There is no atom named "CE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CZ" in residue:
                    self.printout('Warning: There is no atom named "CZ" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")

            if real_resname == "TYR":
                if not "CG" in residue:
                    self.printout('Warning: There is no atom named "CG" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CD1" in residue:
                    self.printout('Warning: There is no atom named "CD1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CD2" in residue:
                    self.printout('Warning: There is no atom named "CD2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CE1" in residue:
                    self.printout('Warning: There is no atom named "CE1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CE2" in residue:
                    self.printout('Warning: There is no atom named "CE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CZ" in residue:
                    self.printout('Warning: There is no atom named "CZ" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")

            if real_resname == "TRP":
                if not "CG" in residue:
                    self.printout('Warning: There is no atom named "CG" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CD1" in residue:
                    self.printout('Warning: There is no atom named "CD1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CD2" in residue:
                    self.printout('Warning: There is no atom named "CD2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "NE1" in residue:
                    self.printout('Warning: There is no atom named "NE1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CE2" in residue:
                    self.printout('Warning: There is no atom named "CE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CE3" in residue:
                    self.printout('Warning: There is no atom named "CE3" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CZ2" in residue:
                    self.printout('Warning: There is no atom named "CZ2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CZ3" in residue:
                    self.printout('Warning: There is no atom named "CZ3" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CH2" in residue:
                    self.printout('Warning: There is no atom named "CH2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")

            if real_resname == "HIS" or real_resname == "HID" or real_resname == "HIE" or real_resname == "HIP":
                if not "CG" in residue:
                    self.printout('Warning: There is no atom named "CG" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "ND1" in residue:
                    self.printout('Warning: There is no atom named "ND1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CD2" in residue:
                    self.printout('Warning: There is no atom named "CD2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "CE1" in residue:
                    self.printout('Warning: There is no atom named "CE1" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")
                if not "NE2" in residue:
                    self.printout('Warning: There is no atom named "NE2" in the protein residue ' + last_key + '. Please use standard naming conventions for all protein residues. This atom is needed to determine pi-pi and pi-cation interactions. If this residue is far from the active site, this warning may not affect the NNScore.')
                    print("")


    # Functions to determine the bond connectivity based on distance
    # ==============================================================

    def CreateBondsByDistance(self):
        for AtomIndex1 in self.NonProteinAtoms:
            atom1 = self.NonProteinAtoms[AtomIndex1]
            if not atom1.residue[-3:] in self.protein_resnames: # so it's not a protein residue
                for AtomIndex2 in self.NonProteinAtoms:
                    if AtomIndex1 != AtomIndex2:
                        atom2 = self.NonProteinAtoms[AtomIndex2]
                        if not atom2.residue[-3:] in self.protein_resnames: # so it's not a protein residue
                            dist = self.functions.distance(atom1.coordinates, atom2.coordinates)

                            if (dist < self.BondLength(atom1.element, atom2.element) * 1.2):
                                atom1.AddNeighborAtomIndex(AtomIndex2)
                                atom2.AddNeighborAtomIndex(AtomIndex1)

    def BondLength(self, element1, element2):

        '''Bond lengths taken from Handbook of Chemistry and Physics. The information provided there was very specific,
        so I tried to pick representative examples and used the bond lengths from those. Sitautions could arise where these
        lengths would be incorrect, probably slight errors (<0.06) in the hundreds.'''

        distance = 0.0
        if element1 == "C" and element2 == "C": distance = 1.53
        if element1 == "N" and element2 == "N": distance = 1.425
        if element1 == "O" and element2 == "O": distance = 1.469
        if element1 == "S" and element2 == "S": distance = 2.048
        if (element1 == "C" and element2 == "H") or (element1 == "H" and element2 == "C"): distance = 1.059
        if (element1 == "C" and element2 == "N") or (element1 == "N" and element2 == "C"): distance = 1.469
        if (element1 == "C" and element2 == "O") or (element1 == "O" and element2 == "C"): distance = 1.413
        if (element1 == "C" and element2 == "S") or (element1 == "S" and element2 == "C"): distance = 1.819
        if (element1 == "N" and element2 == "H") or (element1 == "H" and element2 == "N"): distance = 1.009
        if (element1 == "N" and element2 == "O") or (element1 == "O" and element2 == "N"): distance = 1.463
        if (element1 == "O" and element2 == "S") or (element1 == "S" and element2 == "O"): distance = 1.577
        if (element1 == "O" and element2 == "H") or (element1 == "H" and element2 == "O"): distance = 0.967
        if (element1 == "S" and element2 == "H") or (element1 == "H" and element2 == "S"): distance = 2.025/1.5 # This one not from source cited above. Not sure where it's from, but it wouldn't ever be used in the current context ("AutoGrow")
        if (element1 == "S" and element2 == "N") or (element1 == "H" and element2 == "N"): distance = 1.633

        if (element1 == "C" and element2 == "F") or (element1 == "F" and element2 == "C"): distance = 1.399
        if (element1 == "C" and element2 == "CL") or (element1 == "CL" and element2 == "C"): distance = 1.790
        if (element1 == "C" and element2 == "BR") or (element1 == "BR" and element2 == "C"): distance = 1.910
        if (element1 == "C" and element2 == "I") or (element1 == "I" and element2 == "C"): distance=2.162

        if (element1 == "S" and element2 == "BR") or (element1 == "BR" and element2 == "S"): distance = 2.321
        if (element1 == "S" and element2 == "CL") or (element1 == "CL" and element2 == "S"): distance = 2.283
        if (element1 == "S" and element2 == "F") or (element1 == "F" and element2 == "S"): distance = 1.640
        if (element1 == "S" and element2 == "I") or (element1 == "I" and element2 == "S"): distance= 2.687

        if (element1 == "P" and element2 == "BR") or (element1 == "BR" and element2 == "P"): distance = 2.366
        if (element1 == "P" and element2 == "CL") or (element1 == "CL" and element2 == "P"): distance = 2.008
        if (element1 == "P" and element2 == "F") or (element1 == "F" and element2 == "P"): distance = 1.495
        if (element1 == "P" and element2 == "I") or (element1 == "I" and element2 == "P"): distance= 2.490
        if (element1 == "P" and element2 == "O") or (element1 == "O" and element2 == "P"): distance= 1.6 # estimate based on eye balling Handbook of Chemistry and Physics

        if (element1 == "N" and element2 == "BR") or (element1 == "BR" and element2 == "N"): distance = 1.843
        if (element1 == "N" and element2 == "CL") or (element1 == "CL" and element2 == "N"): distance = 1.743
        if (element1 == "N" and element2 == "F") or (element1 == "F" and element2 == "N"): distance = 1.406
        if (element1 == "N" and element2 == "I") or (element1 == "I" and element2 == "N"): distance= 2.2

        if (element1 == "SI" and element2 == "BR") or (element1 == "BR" and element2 == "SI"): distance = 2.284
        if (element1 == "SI" and element2 == "CL") or (element1 == "CL" and element2 == "SI"): distance = 2.072
        if (element1 == "SI" and element2 == "F") or (element1 == "F" and element2 == "SI"): distance = 1.636
        if (element1 == "SI" and element2 == "P") or (element1 == "P" and element2 == "SI"): distance= 2.264
        if (element1 == "SI" and element2 == "S") or (element1 == "S" and element2 == "SI"): distance= 2.145
        if (element1 == "SI" and element2 == "SI") or (element1 == "SI" and element2 == "SI"): distance= 2.359
        if (element1 == "SI" and element2 == "C") or (element1 == "C" and element2 == "SI"): distance= 1.888
        if (element1 == "SI" and element2 == "N") or (element1 == "N" and element2 == "SI"): distance= 1.743
        if (element1 == "SI" and element2 == "O") or (element1 == "O" and element2 == "SI"): distance= 1.631

        return distance;

    # Functions to identify positive charges
    # ======================================

    def assign_charges(self):
        # Get all the quaternary amines on non-protein residues (these are the only non-protein groups that will be identified as positively charged)
        AllCharged = []
        for atom_index in self.NonProteinAtoms:
            atom = self.NonProteinAtoms[atom_index]
            if atom.element == "MG" or atom.element == "MN" or atom.element == "RH" or atom.element == "ZN" or atom.element == "FE" or atom.element == "BI" or atom.element == "AS" or atom.element == "AG":
                    chrg = self.charged(atom.coordinates, [atom_index], True)
                    self.charges.append(chrg)

            if atom.element == "N":
                if atom.NumberOfNeighbors() == 4: # a quartinary amine, so it's easy
                    indexes = [atom_index]
                    indexes.extend(atom.IndeciesOfAtomsConnecting)
                    chrg = self.charged(atom.coordinates, indexes, True) # so the indicies stored is just the index of the nitrogen and any attached atoms
                    self.charges.append(chrg)
                elif atom.NumberOfNeighbors() == 3: # maybe you only have two hydrogen's added, by they're sp3 hybridized. Just count this as a quartinary amine, since I think the positive charge would be stabalized.
                    nitrogen = atom
                    atom1 = self.AllAtoms[atom.IndeciesOfAtomsConnecting[0]]
                    atom2 = self.AllAtoms[atom.IndeciesOfAtomsConnecting[1]]
                    atom3 = self.AllAtoms[atom.IndeciesOfAtomsConnecting[2]]
                    angle1 = self.functions.angle_between_three_points(atom1.coordinates, nitrogen.coordinates, atom2.coordinates) * 180.0 / math.pi
                    angle2 = self.functions.angle_between_three_points(atom1.coordinates, nitrogen.coordinates, atom3.coordinates) * 180.0 / math.pi
                    angle3 = self.functions.angle_between_three_points(atom2.coordinates, nitrogen.coordinates, atom3.coordinates) * 180.0 / math.pi
                    average_angle = (angle1 + angle2 + angle3) / 3
                    if math.fabs(average_angle - 109.0) < 5.0:
                        indexes = [atom_index]
                        indexes.extend(atom.IndeciesOfAtomsConnecting)
                        chrg = self.charged(nitrogen.coordinates, indexes, True) # so indexes added are the nitrogen and any attached atoms.
                        self.charges.append(chrg)

            if atom.element == "C": # let's check for guanidino-like groups (actually H2N-C-NH2, where not CN3.)
                if atom.NumberOfNeighbors() == 3: # the carbon has only three atoms connected to it
                    nitrogens = self.connected_atoms_of_given_element(atom_index,"N")
                    if len(nitrogens) >= 2: # so carbon is connected to at least two nitrogens
                        # now we need to count the number of nitrogens that are only connected to one heavy atom (the carbon)
                        nitrogens_to_use = []
                        all_connected = atom.IndeciesOfAtomsConnecting[:]
                        not_isolated = -1

                        for atmindex in nitrogens:
                            if len(self.connected_heavy_atoms(atmindex)) == 1:
                                nitrogens_to_use.append(atmindex)
                                all_connected.remove(atmindex)

                        if len(all_connected) > 0: not_isolated = all_connected[0] # get the atom that connects this charged group to the rest of the molecule, ultimately to make sure it's sp3 hybridized

                        if len(nitrogens_to_use) == 2 and not_isolated != -1: # so there are at two nitrogens that are only connected to the carbon (and probably some hydrogens)

                            # now you need to make sure not_isolated atom is sp3 hybridized
                            not_isolated_atom = self.AllAtoms[not_isolated]
                            if (not_isolated_atom.element == "C" and not_isolated_atom.NumberOfNeighbors()==4) or (not_isolated_atom.element == "O" and not_isolated_atom.NumberOfNeighbors()==2) or not_isolated_atom.element == "N" or not_isolated_atom.element == "S" or not_isolated_atom.element == "P":

                                pt = self.AllAtoms[nitrogens_to_use[0]].coordinates.copy_of()
                                pt.x = pt.x + self.AllAtoms[nitrogens_to_use[1]].coordinates.x
                                pt.y = pt.y + self.AllAtoms[nitrogens_to_use[1]].coordinates.y
                                pt.z = pt.z + self.AllAtoms[nitrogens_to_use[1]].coordinates.z
                                pt.x = pt.x / 2.0
                                pt.y = pt.y / 2.0
                                pt.z = pt.z / 2.0

                                indexes = [atom_index]
                                indexes.extend(nitrogens_to_use)
                                indexes.extend(self.connected_atoms_of_given_element(nitrogens_to_use[0],"H"))
                                indexes.extend(self.connected_atoms_of_given_element(nitrogens_to_use[1],"H"))

                                chrg = self.charged(pt, indexes, True) # True because it's positive
                                self.charges.append(chrg)

            if atom.element == "C": # let's check for a carboxylate
                if atom.NumberOfNeighbors() == 3: # a carboxylate carbon will have three items connected to it.
                    oxygens = self.connected_atoms_of_given_element(atom_index,"O")
                    if len(oxygens) == 2: # a carboxylate will have two oxygens connected to it.
                        # now, each of the oxygens should be connected to only one heavy atom (so if it's connected to a hydrogen, that's okay)
                        if len(self.connected_heavy_atoms(oxygens[0])) == 1 and len(self.connected_heavy_atoms(oxygens[1])) == 1:
                            # so it's a carboxylate! Add a negative charge.
                            pt = self.AllAtoms[oxygens[0]].coordinates.copy_of()
                            pt.x = pt.x + self.AllAtoms[oxygens[1]].coordinates.x
                            pt.y = pt.y + self.AllAtoms[oxygens[1]].coordinates.y
                            pt.z = pt.z + self.AllAtoms[oxygens[1]].coordinates.z
                            pt.x = pt.x / 2.0
                            pt.y = pt.y / 2.0
                            pt.z = pt.z / 2.0
                            chrg = self.charged(pt, [oxygens[0], atom_index, oxygens[1]], False)
                            self.charges.append(chrg)

            if atom.element == "P": # let's check for a phosphate or anything where a phosphorus is bound to two oxygens where both oxygens are bound to only one heavy atom (the phosphorus). I think this will get several phosphorus substances.
                oxygens = self.connected_atoms_of_given_element(atom_index,"O")
                if len(oxygens) >=2: # the phosphorus is bound to at least two oxygens
                    # now count the number of oxygens that are only bound to the phosphorus
                    count = 0
                    for oxygen_index in oxygens:
                        if len(self.connected_heavy_atoms(oxygen_index)) == 1: count = count + 1
                    if count >=2: # so there are at least two oxygens that are only bound to the phosphorus
                        indexes = [atom_index]
                        indexes.extend(oxygens)
                        chrg = self.charged(atom.coordinates, indexes, False)
                        self.charges.append(chrg)

            if atom.element == "S": # let's check for a sulfonate or anything where a sulfur is bound to at least three oxygens and at least three are bound to only the sulfur (or the sulfur and a hydrogen).
                oxygens = self.connected_atoms_of_given_element(atom_index,"O")
                if len(oxygens) >=3: # the sulfur is bound to at least three oxygens
                    # now count the number of oxygens that are only bound to the sulfur
                    count = 0
                    for oxygen_index in oxygens:
                        if len(self.connected_heavy_atoms(oxygen_index)) == 1: count = count + 1
                    if count >=3: # so there are at least three oxygens that are only bound to the sulfur
                        indexes = [atom_index]
                        indexes.extend(oxygens)
                        chrg = self.charged(atom.coordinates, indexes, False)
                        self.charges.append(chrg)

        # Now that you've found all the positive charges in non-protein residues, it's time to look for aromatic rings in protein residues
        curr_res = ""
        first = True
        residue = []

        for atom_index in self.AllAtoms:
            atom = self.AllAtoms[atom_index]

            key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain

            if first == True:
                curr_res = key
                first = False

            if key != curr_res:

                self.assign_charged_from_protein_process_residue(residue, last_key)

                residue = []
                curr_res = key

            residue.append(atom_index)
            last_key = key

        self.assign_charged_from_protein_process_residue(residue, last_key)

    def assign_charged_from_protein_process_residue(self, residue, last_key):
        temp = last_key.strip().split("_")
        resname = temp[0]
        real_resname = resname[-3:]
        resid = temp[1]
        chain = temp[2]

        if real_resname == "LYS" or real_resname == "LYN": # regardless of protonation state, assume it's charged.
            for index in residue:
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "NZ":

                    # quickly go through the residue and get the hydrogens attached to this nitrogen to include in the index list
                    indexes = [index]
                    for index2 in residue:
                        atom2 = self.AllAtoms[index2]
                        if atom2.atomname.strip() == "HZ1": indexes.append(index2)
                        if atom2.atomname.strip() == "HZ2": indexes.append(index2)
                        if atom2.atomname.strip() == "HZ3": indexes.append(index2)

                    chrg = self.charged(atom.coordinates, indexes, True)
                    self.charges.append(chrg)
                break

        if real_resname == "ARG":
            charge_pt = point([0.0,0.0,0.0])
            count = 0.0
            indices = []
            for index in residue:
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "NH1":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atomname.strip() == "NH2":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atomname.strip() == "2HH2": indices.append(index)
                if atom.atomname.strip() == "1HH2": indices.append(index)
                if atom.atomname.strip() == "CZ": indices.append(index)
                if atom.atomname.strip() == "2HH1": indices.append(index)
                if atom.atomname.strip() == "1HH1": indices.append(index)

            if count != 0.0:

                charge_pt.x = charge_pt.x / count
                charge_pt.y = charge_pt.y / count
                charge_pt.z = charge_pt.z / count

                if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
                    chrg = self.charged(charge_pt, indices, True)
                    self.charges.append(chrg)

        if real_resname == "HIS" or real_resname == "HID" or real_resname == "HIE" or real_resname == "HIP": # regardless of protonation state, assume it's charged. This based on "The Cation-Pi Interaction," which suggests protonated state would be stabalized. But let's not consider HIS when doing salt bridges.
            charge_pt = point([0.0,0.0,0.0])
            count = 0.0
            indices = []
            for index in residue:
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "NE2":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atomname.strip() == "ND1":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atomname.strip() == "HE2": indices.append(index)
                if atom.atomname.strip() == "HD1": indices.append(index)
                if atom.atomname.strip() == "CE1": indices.append(index)
                if atom.atomname.strip() == "CD2": indices.append(index)
                if atom.atomname.strip() == "CG": indices.append(index)

            if count != 0.0:
                charge_pt.x = charge_pt.x / count
                charge_pt.y = charge_pt.y / count
                charge_pt.z = charge_pt.z / count
                if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
                    chrg = self.charged(charge_pt, indices, True)
                    self.charges.append(chrg)

        if real_resname == "GLU" or real_resname == "GLH" or real_resname == "GLX": # regardless of protonation state, assume it's charged. This based on "The Cation-Pi Interaction," which suggests protonated state would be stabalized.
            charge_pt = point([0.0,0.0,0.0])
            count = 0.0
            indices = []
            for index in residue:
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "OE1":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atomname.strip() == "OE2":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atomname.strip() == "CD": indices.append(index)

            if count != 0.0:
                charge_pt.x = charge_pt.x / count
                charge_pt.y = charge_pt.y / count
                charge_pt.z = charge_pt.z / count
                if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
                    chrg = self.charged(charge_pt, indices, False) # False because it's a negative charge
                    self.charges.append(chrg)

        if real_resname == "ASP" or real_resname == "ASH" or real_resname == "ASX": # regardless of protonation state, assume it's charged. This based on "The Cation-Pi Interaction," which suggests protonated state would be stabalized.
            charge_pt = point([0.0,0.0,0.0])
            count = 0.0
            indices = []
            for index in residue:
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "OD1":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atomname.strip() == "OD2":
                    charge_pt.x = charge_pt.x + atom.coordinates.x
                    charge_pt.y = charge_pt.y + atom.coordinates.y
                    charge_pt.z = charge_pt.z + atom.coordinates.z
                    indices.append(index)
                    count = count + 1.0
                if atom.atomname.strip() == "CG": indices.append(index)

            if count != 0.0:
                charge_pt.x = charge_pt.x / count
                charge_pt.y = charge_pt.y / count
                charge_pt.z = charge_pt.z / count
                if charge_pt.x != 0.0 or charge_pt.y != 0.0 or charge_pt.z != 0.0:
                    chrg = self.charged(charge_pt, indices, False) # False because it's a negative charge
                    self.charges.append(chrg)

    class charged():
        def __init__(self, coordinates, indices, positive):
            self.coordinates = coordinates
            self.indices = indices
            self.positive = positive # true or false to specifiy if positive or negative charge

    # Functions to identify aromatic rings
    # ====================================

    def add_aromatic_marker(self, indicies_of_ring):
        # first identify the center point
        points_list = []
        total = len(indicies_of_ring)
        x_sum = 0.0
        y_sum = 0.0
        z_sum = 0.0

        for index in indicies_of_ring:
            atom = self.AllAtoms[index]
            points_list.append(atom.coordinates)
            x_sum = x_sum + atom.coordinates.x
            y_sum = y_sum + atom.coordinates.y
            z_sum = z_sum + atom.coordinates.z

        if total == 0: return # to prevent errors in some cases

        center = point([x_sum / total, y_sum / total, z_sum / total])

        # now get the radius of the aromatic ring
        radius = 0.0
        for index in indicies_of_ring:
            atom = self.AllAtoms[index]
            dist = center.dist_to(atom.coordinates)
            if dist > radius: radius = dist

        # now get the plane that defines this ring
        if len(indicies_of_ring) < 3:
            return # to prevent an error in some cases. If there aren't three point, you can't define a plane
        elif len(indicies_of_ring) == 3:
            A = self.AllAtoms[indicies_of_ring[0]].coordinates
            B = self.AllAtoms[indicies_of_ring[1]].coordinates
            C = self.AllAtoms[indicies_of_ring[2]].coordinates
        elif len(indicies_of_ring) == 4:
            A = self.AllAtoms[indicies_of_ring[0]].coordinates
            B = self.AllAtoms[indicies_of_ring[1]].coordinates
            C = self.AllAtoms[indicies_of_ring[3]].coordinates
        else: # best, for 5 and 6 member rings
            A = self.AllAtoms[indicies_of_ring[0]].coordinates
            B = self.AllAtoms[indicies_of_ring[2]].coordinates
            C = self.AllAtoms[indicies_of_ring[4]].coordinates

        AB = self.functions.vector_subtraction(B,A)
        AC = self.functions.vector_subtraction(C,A)
        ABXAC = self.functions.CrossProduct(AB,AC)

        # formula for plane will be ax + by + cz = d
        x1 = self.AllAtoms[indicies_of_ring[0]].coordinates.x
        y1 = self.AllAtoms[indicies_of_ring[0]].coordinates.y
        z1 = self.AllAtoms[indicies_of_ring[0]].coordinates.z

        a = ABXAC.x
        b = ABXAC.y
        c = ABXAC.z
        d = a*x1 + b*y1 + c*z1

        ar_ring = self.aromatic_ring(center, indicies_of_ring, [a,b,c,d], radius)
        self.aromatic_rings.append(ar_ring)

    class aromatic_ring():
        def __init__(self, center, indices, plane_coeff, radius):
            self.center = center
            self.indices = indices
            self.plane_coeff = plane_coeff # a*x + b*y + c*z = dI think that
            self.radius = radius

    def assign_aromatic_rings(self):
        # Get all the rings containing each of the atoms in the ligand
        AllRings = []
        for atom_index in self.NonProteinAtoms:
            AllRings.extend(self.all_rings_containing_atom(atom_index))

        for ring_index_1 in range(len(AllRings)):
            ring1 = AllRings[ring_index_1]
            if len(ring1) != 0:
                for ring_index_2 in range(len(AllRings)):
                    if ring_index_1 != ring_index_2:
                        ring2 = AllRings[ring_index_2]
                        if len(ring2) != 0:
                            if self.set1_is_subset_of_set2(ring1, ring2) == True: AllRings[ring_index_2] = []

        while [] in AllRings: AllRings.remove([])

        # Now we need to figure out which of these ligands are aromatic (planar)

        for ring_index in range(len(AllRings)):
            ring = AllRings[ring_index]
            is_flat = True
            for t in range(-3, len(ring)-3):
                pt1 = self.NonProteinAtoms[ring[t]].coordinates
                pt2 = self.NonProteinAtoms[ring[t+1]].coordinates
                pt3 = self.NonProteinAtoms[ring[t+2]].coordinates
                pt4 = self.NonProteinAtoms[ring[t+3]].coordinates

                # first, let's see if the last atom in this ring is a carbon connected to four atoms. That would be a quick way of telling this is not an aromatic ring
                cur_atom = self.NonProteinAtoms[ring[t+3]]
                if cur_atom.element == "C" and cur_atom.NumberOfNeighbors() == 4:
                    is_flat = False
                    break

                # now check the dihedral between the ring atoms to see if it's flat
                angle = self.functions.dihedral(pt1, pt2, pt3, pt4) * 180 / math.pi
                if (angle > -165 and angle < -15) or (angle > 15 and angle < 165): # 15 degress is the cutoff #, ring[t], ring[t+1], ring[t+2], ring[t+3] # range of this function is -pi to pi
                    is_flat = False
                    break

                # now check the dihedral between the ring atoms and an atom connected to the current atom to see if that's flat too.
                for substituent_atom_index in cur_atom.IndeciesOfAtomsConnecting:
                    pt_sub = self.NonProteinAtoms[substituent_atom_index].coordinates
                    angle = self.functions.dihedral(pt2, pt3, pt4, pt_sub) * 180 / math.pi
                    if (angle > -165 and angle < -15) or (angle > 15 and angle < 165): # 15 degress is the cutoff #, ring[t], ring[t+1], ring[t+2], ring[t+3] # range of this function is -pi to pi
                        is_flat = False
                        break

            if is_flat == False: AllRings[ring_index] = []
            if len(ring) < 5: AllRings[ring_index] = [] # While I'm at it, three and four member rings are not aromatic
            if len(ring) > 6: AllRings[ring_index] = [] # While I'm at it, if the ring has more than 6, also throw it out. So only 5 and 6 member rings are allowed.



        while [] in AllRings: AllRings.remove([])

        for ring in AllRings:
            self.add_aromatic_marker(ring)

        # Now that you've found all the rings in non-protein residues, it's time to look for aromatic rings in protein residues
        curr_res = ""
        first = True
        residue = []

        for atom_index in self.AllAtoms:
            atom = self.AllAtoms[atom_index]

            key = atom.residue + "_" + str(atom.resid) + "_" + atom.chain

            if first == True:
                curr_res = key
                first = False

            if key != curr_res:

                self.assign_aromatic_rings_from_protein_process_residue(residue, last_key)

                residue = []
                curr_res = key

            residue.append(atom_index)
            last_key = key

        self.assign_aromatic_rings_from_protein_process_residue(residue, last_key)

    def assign_aromatic_rings_from_protein_process_residue(self, residue, last_key):
        temp = last_key.strip().split("_")
        resname = temp[0]
        real_resname = resname[-3:]
        resid = temp[1]
        chain = temp[2]

        if real_resname == "PHE":
            indicies_of_ring = []

            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CG": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CD1": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CE1": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CZ": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CE2": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CD2": indicies_of_ring.append(index)

            self.add_aromatic_marker(indicies_of_ring)

        if real_resname == "TYR":
            indicies_of_ring = []

            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CG": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CD1": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CE1": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CZ": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CE2": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CD2": indicies_of_ring.append(index)

            self.add_aromatic_marker(indicies_of_ring)

        if real_resname == "HIS" or real_resname == "HID" or real_resname == "HIE" or real_resname == "HIP":
            indicies_of_ring = []

            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CG": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "ND1": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CE1": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "NE2": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CD2": indicies_of_ring.append(index)

            self.add_aromatic_marker(indicies_of_ring)

        if real_resname == "TRP":
            indicies_of_ring = []

            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CG": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CD1": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "NE1": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CE2": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CD2": indicies_of_ring.append(index)

            self.add_aromatic_marker(indicies_of_ring)

            indicies_of_ring = []

            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CE2": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CD2": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CE3": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CZ3": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CH2": indicies_of_ring.append(index)
            for index in residue: # written this way because order is important
                atom = self.AllAtoms[index]
                if atom.atomname.strip() == "CZ2": indicies_of_ring.append(index)

            self.add_aromatic_marker(indicies_of_ring)

    def set1_is_subset_of_set2(self, set1, set2):
        is_subset = True
        for item in set1:
            if not item in set2:
                is_subset = False
                break
        return is_subset

    def all_rings_containing_atom(self, index):

        AllRings = []

        atom = self.AllAtoms[index]
        for conneceted_atom in atom.IndeciesOfAtomsConnecting:
            self.ring_recursive(conneceted_atom, [index], index, AllRings)

        return AllRings

    def ring_recursive(self, index, AlreadyCrossed, orig_atom, AllRings):

        if len(AlreadyCrossed) > 6: return # since you're only considering aromatic rings containing 5 or 6 members anyway, save yourself some time.

        atom = self.AllAtoms[index]

        temp = AlreadyCrossed[:]
        temp.append(index)

        for conneceted_atom in atom.IndeciesOfAtomsConnecting:
            if not conneceted_atom in AlreadyCrossed:
                self.ring_recursive(conneceted_atom, temp, orig_atom, AllRings)
            if conneceted_atom == orig_atom and orig_atom != AlreadyCrossed[-1]:
                AllRings.append(temp)

    # Functions to assign secondary structure to protein residues
    # ===========================================================

    def assign_secondary_structure(self):
        # first, we need to know what resid's are available
        resids = []
        last_key = "-99999_Z"
        for atom_index in self.AllAtoms:
            atom = self.AllAtoms[atom_index]
            key = str(atom.resid) + "_" + atom.chain
            if key != last_key:
                last_key = key
                resids.append(last_key)

        structure = {}
        for resid in resids:
            structure[resid] = "OTHER"

        atoms = []

        for atom_index in self.AllAtoms:
            atom = self.AllAtoms[atom_index]
            if atom.SideChainOrBackBone() == "BACKBONE":
                if len(atoms) < 8:
                    atoms.append(atom)
                else:
                    atoms.pop(0)
                    atoms.append(atom)

                    # now make sure the first four all have the same resid and the last four all have the same resid
                    if atoms[0].resid == atoms[1].resid and atoms[0].resid == atoms[2].resid and atoms[0].resid == atoms[3].resid and atoms[0] != atoms[4].resid and atoms[4].resid == atoms[5].resid and atoms[4].resid == atoms[6].resid and atoms[4].resid == atoms[7].resid and atoms[0].resid + 1 == atoms[7].resid and atoms[0].chain == atoms[7].chain:
                        resid1 = atoms[0].resid
                        resid2 = atoms[7].resid

                        # Now give easier to use names to the atoms
                        for atom in atoms:
                            if atom.resid == resid1 and atom.atomname.strip() == "N": first_N = atom
                            if atom.resid == resid1 and atom.atomname.strip() == "C": first_C = atom
                            if atom.resid == resid1 and atom.atomname.strip() == "CA": first_CA = atom

                            if atom.resid == resid2 and atom.atomname.strip() == "N": second_N = atom
                            if atom.resid == resid2 and atom.atomname.strip() == "C": second_C = atom
                            if atom.resid == resid2 and atom.atomname.strip() == "CA": second_CA = atom

                        # Now compute the phi and psi dihedral angles
                        phi = self.functions.dihedral(first_C.coordinates, second_N.coordinates, second_CA.coordinates, second_C.coordinates) * 180.0 / math.pi
                        psi = self.functions.dihedral(first_N.coordinates, first_CA.coordinates, first_C.coordinates, second_N.coordinates) * 180.0 / math.pi

                        # Now use those angles to determine if it's alpha or beta
                        if phi > -145 and phi < -35 and psi > -70 and psi < 50:
                            key1 = str(first_C.resid) + "_" + first_C.chain
                            key2 = str(second_C.resid) + "_" + second_C.chain
                            structure[key1] = "ALPHA"
                            structure[key2] = "ALPHA"
                        if (phi >= -180 and phi < -40 and psi <= 180 and psi > 90) or (phi >= -180 and phi < -70 and psi <= -165): # beta. This gets some loops (by my eye), but it's the best I could do.
                            key1 = str(first_C.resid) + "_" + first_C.chain
                            key2 = str(second_C.resid) + "_" + second_C.chain
                            structure[key1] = "BETA"
                            structure[key2] = "BETA"

        # Now update each of the atoms with this structural information
        for atom_index in self.AllAtoms:
            atom = self.AllAtoms[atom_index]
            key = str(atom.resid) + "_" + atom.chain
            atom.structure = structure[key]

        # Some more post processing.
        CA_list = [] # first build a list of the indices of all the alpha carbons
        for atom_index in self.AllAtoms:
            atom = self.AllAtoms[atom_index]
            if atom.residue.strip() in self.protein_resnames and atom.atomname.strip() == "CA": CA_list.append(atom_index)

        # some more post processing.
        change = True
        while change == True:
            change = False

            # A residue of index i is only going to be in an alpha helix its CA is within 6 A of the CA of the residue i + 3
            for CA_atom_index in CA_list:
                CA_atom = self.AllAtoms[CA_atom_index]
                if CA_atom.structure == "ALPHA": # so it's in an alpha helix
                    another_alpha_is_close = False
                    for other_CA_atom_index in CA_list: # so now compare that CA to all the other CA's
                        other_CA_atom = self.AllAtoms[other_CA_atom_index]
                        if other_CA_atom.structure == "ALPHA": # so it's also in an alpha helix
                            if other_CA_atom.resid - 3 == CA_atom.resid or other_CA_atom.resid + 3 == CA_atom.resid: # so this CA atom is one of the ones the first atom might hydrogen bond with
                                if other_CA_atom.coordinates.dist_to(CA_atom.coordinates) < 6.0: # so these two CA atoms are close enough together that their residues are probably hydrogen bonded
                                    another_alpha_is_close = True
                                    break
                    if another_alpha_is_close == False:
                        self.set_structure_of_residue(CA_atom.chain, CA_atom.resid, "OTHER")
                        change = True

            # Alpha helices are only alpha helices if they span at least 4 residues (to wrap around and hydrogen bond). I'm going to require them to span at least 5 residues, based on examination of many structures.
            for index_in_list in range(len(CA_list)-5):

                index_in_pdb1 = CA_list[index_in_list]
                index_in_pdb2 = CA_list[index_in_list+1]
                index_in_pdb3 = CA_list[index_in_list+2]
                index_in_pdb4 = CA_list[index_in_list+3]
                index_in_pdb5 = CA_list[index_in_list+4]
                index_in_pdb6 = CA_list[index_in_list+5]

                atom1 = self.AllAtoms[index_in_pdb1]
                atom2 = self.AllAtoms[index_in_pdb2]
                atom3 = self.AllAtoms[index_in_pdb3]
                atom4 = self.AllAtoms[index_in_pdb4]
                atom5 = self.AllAtoms[index_in_pdb5]
                atom6 = self.AllAtoms[index_in_pdb6]

                if atom1.resid + 1 == atom2.resid and atom2.resid + 1 == atom3.resid and atom3.resid + 1 == atom4.resid and atom4.resid + 1 == atom5.resid and atom5.resid + 1 == atom6.resid: # so they are sequential

                    if atom1.structure != "ALPHA" and atom2.structure == "ALPHA" and atom3.structure != "ALPHA":
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        change = True
                    if atom2.structure != "ALPHA" and atom3.structure == "ALPHA" and atom4.structure != "ALPHA":
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        change = True
                    if atom3.structure != "ALPHA" and atom4.structure == "ALPHA" and atom5.structure != "ALPHA":
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        change = True
                    if atom4.structure != "ALPHA" and atom5.structure == "ALPHA" and atom6.structure != "ALPHA":
                        self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
                        change = True

                    if atom1.structure != "ALPHA" and atom2.structure == "ALPHA" and atom3.structure == "ALPHA" and atom4.structure != "ALPHA":
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        change = True
                    if atom2.structure != "ALPHA" and atom3.structure == "ALPHA" and atom4.structure == "ALPHA" and atom5.structure != "ALPHA":
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        change = True
                    if atom3.structure != "ALPHA" and atom4.structure == "ALPHA" and atom5.structure == "ALPHA" and atom6.structure != "ALPHA":
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
                        change = True

                    if atom1.structure != "ALPHA" and atom2.structure == "ALPHA" and atom3.structure == "ALPHA" and atom4.structure == "ALPHA" and atom5.structure != "ALPHA":
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        change = True
                    if atom2.structure != "ALPHA" and atom3.structure == "ALPHA" and atom4.structure == "ALPHA" and atom5.structure == "ALPHA" and atom6.structure != "ALPHA":
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
                        change = True

                    if atom1.structure != "ALPHA" and atom2.structure == "ALPHA" and atom3.structure == "ALPHA" and atom4.structure == "ALPHA" and atom5.structure == "ALPHA" and atom6.structure != "ALPHA":
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        self.set_structure_of_residue(atom4.chain, atom4.resid, "OTHER")
                        self.set_structure_of_residue(atom5.chain, atom5.resid, "OTHER")
                        change = True

            # now go through each of the BETA CA atoms. A residue is only going to be called a beta sheet if CA atom is within 6.0 A of another CA beta, same chain, but index difference > 2.
            for CA_atom_index in CA_list:
                CA_atom = self.AllAtoms[CA_atom_index]
                if CA_atom.structure == "BETA": # so it's in a beta sheet
                    another_beta_is_close = False
                    for other_CA_atom_index in CA_list:
                        if other_CA_atom_index != CA_atom_index: # so not comparing an atom to itself
                            other_CA_atom = self.AllAtoms[other_CA_atom_index]
                            if other_CA_atom.structure == "BETA": # so you're comparing it only to other BETA-sheet atoms
                                if other_CA_atom.chain == CA_atom.chain: # so require them to be on the same chain. needed to indecies can be fairly compared
                                    if math.fabs(other_CA_atom.resid - CA_atom.resid) > 2: # so the two residues are not simply adjacent to each other on the chain
                                        if CA_atom.coordinates.dist_to(other_CA_atom.coordinates) < 6.0: # so these to atoms are close to each other
                                            another_beta_is_close = True
                                            break
                    if another_beta_is_close == False:
                        self.set_structure_of_residue(CA_atom.chain, CA_atom.resid, "OTHER")
                        change = True

            # Now some more post-processing needs to be done. Do this again to clear up mess that may have just been created (single residue beta strand, for example)
            # Beta sheets are usually at least 3 residues long

            for index_in_list in range(len(CA_list)-3):

                index_in_pdb1 = CA_list[index_in_list]
                index_in_pdb2 = CA_list[index_in_list+1]
                index_in_pdb3 = CA_list[index_in_list+2]
                index_in_pdb4 = CA_list[index_in_list+3]

                atom1 = self.AllAtoms[index_in_pdb1]
                atom2 = self.AllAtoms[index_in_pdb2]
                atom3 = self.AllAtoms[index_in_pdb3]
                atom4 = self.AllAtoms[index_in_pdb4]

                if atom1.resid + 1 == atom2.resid and atom2.resid + 1 == atom3.resid and atom3.resid + 1 == atom4.resid: # so they are sequential

                    if atom1.structure != "BETA" and atom2.structure == "BETA" and atom3.structure != "BETA":
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        change = True
                    if atom2.structure != "BETA" and atom3.structure == "BETA" and atom4.structure != "BETA":
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        change = True
                    if atom1.structure != "BETA" and atom2.structure == "BETA" and atom3.structure == "BETA" and atom4.structure != "BETA":
                        self.set_structure_of_residue(atom2.chain, atom2.resid, "OTHER")
                        self.set_structure_of_residue(atom3.chain, atom3.resid, "OTHER")
                        change = True

    def set_structure_of_residue(self, chain, resid, structure):
        for atom_index in self.AllAtoms:
            atom = self.AllAtoms[atom_index]
            if atom.chain == chain and atom.resid == resid:
                atom.structure = structure

class MathFunctions:
    def planarity(self, point1, point2, point3, point4):

            x1 = point1.x
            y1 = point1.y
            z1 = point1.z
            x2 = point2.x
            y2 = point2.y
            z2 = point2.z
            x3 = point3.x
            y3 = point3.y
            z3 = point3.z
            x4 = point4.x
            y4 = point4.y
            z4 = point4.z


            A = (y1*(z2-z3))+(y2*(z3-z1))+(y3*(z1-z2))
            B = (z1*(x2-x3))+(z2*(x3-x1))+(z3*(x1-x2))
            C = (x1*(y2-y3))+(x2*(y3-y1))+(x3*(y1-y2))
            D = ((-x1)*((y2*z3)-(y3*z2)))+((-x2)*((y3*z1)-(y1*z3)))+((-x3)*((y1*z2)-(y2*z1)))
            distance=(math.fabs((A*x4)+(B*y4)+(C*z4)+D))/(math.sqrt(math.pow(A,2) + math.pow(B,2) + math.pow(C,2)))

            A1 = (y1*(z2-z4))+(y2*(z4-z1))+(y4*(z1-z2))
            B1 = (z1*(x2-x4))+(z2*(x4-x1))+(z4*(x1-x2))
            C1 = (x1*(y2-y4))+(x2*(y4-y1))+(x4*(y1-y2))
            D1 = ((-x1)*((y2*z4)-(y4*z2)))+((-x2)*((y4*z1)-(y1*z4)))+((-x4)*((y1*z2)-(y2*z1)))
            distance1=(math.fabs((A1*x3)+(B1*y3)+(C1*z3)+D1))/(math.sqrt(math.pow(A1,2) + math.pow(B1,2) + math.pow(C1,2)))

            A2 = (y1*(z4-z3))+(y4*(z3-z1))+(y3*(z1-z4))
            B2 = (z1*(x4-x3))+(z4*(x3-x1))+(z3*(x1-x4))
            C2 = (x1*(y4-y3))+(x4*(y3-y1))+(x3*(y1-y4))
            D2 = ((-x1)*((y4*z3)-(y3*z4)))+((-x4)*((y3*z1)-(y1*z3)))+((-x3)*((y1*z4)-(y4*z1)))
            distance2=(math.fabs((A2*x2)+(B2*y2)+(C2*z2)+D2))/(math.sqrt(math.pow(A2,2) + math.pow(B2,2) + math.pow(C2,2)))

            A3 = (y4*(z2-z3))+(y2*(z3-z4))+(y3*(z4-z2))
            B3 = (z4*(x2-x3))+(z2*(x3-x4))+(z3*(x4-x2))
            C3 = (x4*(y2-y3))+(x2*(y3-y4))+(x3*(y4-y2))
            D3 = ((-x4)*((y2*z3)-(y3*z2)))+((-x2)*((y3*z4)-(y4*z3)))+((-x3)*((y4*z2)-(y2*z4)))
            distance3=(math.fabs((A3*x1)+(B3*y1)+(C3*z1)+D3))/(math.sqrt(math.pow(A3,2) + math.pow(B3,2) + math.pow(C3,2)))

            final_dist = -1

            if (distance < distance1 and distance < distance2 and distance < distance3):
                    final_dist = distance
            elif (distance1 < distance and distance1 < distance2 and distance1 < distance3):
                    final_dist = distance1
            elif (distance2 < distance and distance2 < distance1 and distance2 < distance3):
                    final_dist = distance2
            elif (distance3 < distance and distance3 < distance1 and distance3 < distance2):
                    final_dist = distance3

            # Now normalize by the length of the longest bond

            return final_dist

    def vector_subtraction(self, vector1, vector2): # vector1 - vector2
        return point([vector1.x - vector2.x, vector1.y - vector2.y, vector1.z - vector2.z])

    def CrossProduct(self, Pt1, Pt2): # never tested
        Response = point([0,0,0])

        Response.x = Pt1.y * Pt2.z - Pt1.z * Pt2.y
        Response.y = Pt1.z * Pt2.x - Pt1.x * Pt2.z
        Response.z = Pt1.x * Pt2.y - Pt1.y * Pt2.x

        return Response;

    def vector_scalar_multiply(self, vector, scalar):
        return point([vector.x * scalar, vector.y * scalar, vector.z * scalar])

    def dot_product(self, point1, point2):
        return point1.x * point2.x + point1.y * point2.y + point1.z * point2.z

    def dihedral(self, point1, point2, point3, point4): # never tested

        b1 = self.vector_subtraction(point2, point1)
        b2 = self.vector_subtraction(point3, point2)
        b3 = self.vector_subtraction(point4, point3)

        b2Xb3 = self.CrossProduct(b2,b3)
        b1Xb2 = self.CrossProduct(b1,b2)

        b1XMagb2 = self.vector_scalar_multiply(b1,b2.magnitude())
        radians = math.atan2(self.dot_product(b1XMagb2,b2Xb3), self.dot_product(b1Xb2,b2Xb3))
        return radians

    def angle_between_three_points(self, point1, point2, point3): # As in three connected atoms
        vector1 = self.vector_subtraction(point1, point2)
        vector2 = self.vector_subtraction(point3, point2)
        return self.angle_between_points(vector1, vector2)

    def angle_between_points(self, point1, point2):
        new_point1 = self.return_normalized_vector(point1)
        new_point2 = self.return_normalized_vector(point2)
        dot_prod = self.dot_product(new_point1, new_point2)
        if dot_prod > 1.0: dot_prod = 1.0 # to prevent errors that can rarely occur
        if dot_prod < -1.0: dot_prod = -1.0
        return math.acos(dot_prod)

    def return_normalized_vector(self, vector):
        dist = self.distance(point([0,0,0]), vector)
        return point([vector.x/dist, vector.y/dist, vector.z/dist])

    def distance(self, point1, point2):
        deltax = point1.x - point2.x
        deltay = point1.y - point2.y
        deltaz = point1.z - point2.z

        return math.sqrt(math.pow(deltax,2) + math.pow(deltay,2) + math.pow(deltaz,2))

    def project_point_onto_plane(self, apoint, plane_coefficients): # essentially finds the point on the plane that is closest to the specified point
        # the plane_coefficients are [a,b,c,d], where the plane is ax + by + cz = d

        # First, define a plane using cooeficients a, b, c, d such that ax + by + cz = d
        a = plane_coefficients[0]
        b = plane_coefficients[1]
        c = plane_coefficients[2]
        d = plane_coefficients[3]

        # Now, define a point in space (s,u,v)
        s = apoint.x
        u = apoint.y
        v = apoint.z

        # the formula of a line perpendicular to the plan passing through (s,u,v) is:
        #x = s + at
        #y = u + bt
        #z = v + ct

        t = (d - a*s - b*u - c*v) / (a*a + b*b + c*c)

        # here's the point closest on the plane
        x = s + a*t
        y = u + b*t
        z = v + c*t

        return point(x,y,z)

    ### Simple quaternion math using tuples
    ## from http://stackoverflow.com/questions/4870393/rotating-coordinate-system-via-a-quaternion
    def normalize(self, v, tolerance=0.00001):
        mag2 = sum(n * n for n in v)
        if abs(mag2 - 1.0) > tolerance:
            mag = numpy.power(mag2, 0.5)
            v = tuple(n / mag for n in v)
        return v

    def q_mult(self, q1, q2):
        w1, x1, y1, z1 = q1
        w2, x2, y2, z2 = q2
        w = w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2
        x = w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2
        y = w1 * y2 + y1 * w2 + z1 * x2 - x1 * z2
        z = w1 * z2 + z1 * w2 + x1 * y2 - y1 * x2
        return w, x, y, z


    def q_conjugate(self, q):
        q = self.normalize(q)
        w, x, y, z = q
        return (w, -x, -y, -z)

    def qv_mult(self, q1, v1):
        v1 = self.normalize(v1)
        q2 = (0.0,) + v1
        return self.q_mult(self.q_mult(q1, q2), self.q_conjugate(q1))[1:]

    def axisangle_to_q(self, v, theta):
        v = self.normalize(v)
        x, y, z = v
        theta /= 2
        w =     numpy.cos(theta)
        x = x * numpy.sin(theta)
        y = y * numpy.sin(theta)
        z =     z * numpy.sin(theta)
        return w, x, y, z

    def q_to_axisangle(self, q):
        w, v = q[0], q[1:]
        theta = numpy.arccos(w) * 2.0
        return self.normalize(v), theta




class featureMap:
    def __init__(self, borders, reso, boolean = False):
        #self.origin = point([float(borders[0])*reso,
        #                     float(borders[2])*reso,
        #                     float(borders[4])*reso])
        self.origin = point([float(borders[0]),
                             float(borders[2]),
                             float(borders[4])])
        #self.borders = borders
        self.reso = reso
        #self.shape = (int(numpy.round(float(math.fabs(borders[1]-borders[0]))/reso))+1,
        #              int(numpy.round(float(math.fabs(borders[3]-borders[2]))/reso))+1,
        #              int(numpy.round(float(math.fabs(borders[5]-borders[4]))/reso))+1)
        # Replacing numpy.round because it's super slow
        self.shape = (int(0.5+(float(math.fabs(borders[1]-borders[0]))/reso))+1,
                      int(0.5+(float(math.fabs(borders[3]-borders[2]))/reso))+1,
                      int(0.5+(float(math.fabs(borders[5]-borders[4]))/reso))+1)

        self.borders = [self.origin[0],
                        self.origin[0]+(self.reso*(self.shape[0]-1)), #We subtract one because the matrix has grid points on top of both the max and min in this dimension
                        self.origin[1],
                        self.origin[1]+(self.reso*(self.shape[1]-1)),
                        self.origin[2],
                        self.origin[2]+(self.reso*(self.shape[2]-1))]

        self.boolean = boolean
        if boolean:
            self.data = numpy.zeros(self.shape, dtype = numpy.bool)
        else:
            self.data = numpy.zeros(self.shape, dtype = numpy.float)

        self.mf = MathFunctions()
        # This is generated when the program wants to add spherical gaussians, but must be initialized here
        self.gaussianLookupTable = None

    @classmethod
    def from3dArray(cls, myArray, gridReso, origin = [0.,0.,0.]):
        shape = myArray.shape
        borders = [origin[0],
                   (origin[0]+shape[0]*gridReso),
                   origin[1],
                   (origin[1]+shape[1]*gridReso),
                   origin[2],
                   (origin[2]+shape[2]*gridReso)]

        thisItem = cls(borders, gridReso)
        thisItem.data = numpy.copy(myArray)
        thisItem.shape = thisItem.data.shape
        return thisItem


    @classmethod
    def fromPovmeList(cls, povmeList, gridReso):
        minX = numpy.min(povmeList[:,0])
        maxX = numpy.max(povmeList[:,0])# + 1e-5
        minY = numpy.min(povmeList[:,1])
        maxY = numpy.max(povmeList[:,1])# + 1e-5
        minZ = numpy.min(povmeList[:,2])
        maxZ = numpy.max(povmeList[:,2])# + 1e-5
        #print [minX, maxX, minY, maxY, minZ, maxZ]
        thisMap = cls([minX, maxX, minY, maxY, minZ, maxZ], gridReso)
        referencePt = povmeList[0,:3]
        referenceInd = thisMap.point_to_nearest_index(referencePt)
        shiftedCoordList = copy.deepcopy(povmeList)
        shiftedCoordList[:,:3] = povmeList[:,:3] - referencePt
        shiftedIndexList = copy.deepcopy(shiftedCoordList)
        shiftedIndexList[:,:3] = (shiftedIndexList[:,:3] / thisMap.getReso()) + referenceInd
        for shiftedIndex in shiftedIndexList:
            thisMap.data[tuple(shiftedIndex[:3])] = shiftedIndex[3]
        #for point in povmeList:
        #    offsetCoord = point[:3]-referencePt
        #    offsetInd = offsetCoord/thisMap.getReso()
        #    index = referenceInd + offsetInd
        #    index = tuple([int(i) for i in index])
            #index = thisMap.point_to_nearest_index(point[:3])
            #print "referencePt, referenceInd, offsetCoord, offsetInd, point, index, thisMap.data.shape"
            #print referencePt, referenceInd, offsetCoord, offsetInd, point, index, thisMap.data.shape
        #    thisMap.data[index] = point[3]
        return thisMap



    def toPovmeList(self):
        nonZeroPts = numpy.transpose(self.data.nonzero())
        #print nonZeroPts.shape
        pointList = numpy.zeros((len(nonZeroPts), 4))
        for index, nonZeroPt in enumerate(nonZeroPts):
            coords = self.index_to_coord(nonZeroPt[0],nonZeroPt[1],nonZeroPt[2])
            pointList[index,0] = coords[0]
            pointList[index,1] = coords[1]
            pointList[index,2] = coords[2]
            #print 'nonZeroPt', nonZeroPt, 'self.data[nonZeroPt]', self.data[nonZeroPt]
            pointList[index,3] = self.data[nonZeroPt[0], nonZeroPt[1], nonZeroPt[2]]
        return pointList

    def getOrigin(self):
        return self.origin

    def getReso(self):
        return self.reso

    def getShape(self):
        return self.shape

    def getData(self):
        return self.data

    def setData(self, newMatrix):
        self.data = numpy.copy(newMatrix)

    def getBoolean(self):
        return self.boolean

    def getBorders(self):
        return self.borders



    def generateTranslations(self,ligandMap,resolution,number):
        #Need to add borders?
        #User should input size of cubes, not number of cubes
        ligandDim = ligandMap.getShape()
        receptorDim = self.getShape()
        receptorSize = numpy.ndarray.size(self.getData())
        cubeLength = (receptorSize / number) ** .33333333
        x = []
        y = []
        z = []
        for n in numpy.arange(0,receptorDim[0]-ligandDim[0],cubeLength):
            x.append(n)
        for n in numpy.arange(0,receptorDim[1]-ligandDim[1],cubeLength):
            y.append(n)
        for n in numpy.arange(0,receptorDim[2]-ligandDim[2],cubeLength):
            z.append(n)

        import itertools
        translation_list = list(itertools.product(x, y, z))

        return translation_list

    def generateRotations(self,number,nsteps):

        if number in list(self.pointDict.keys()):
            spherePoints = self.pointDict[number]

        else:
            number = str(number)
            nsteps = str(nsteps)
            outputFile = number + "_" + nsteps
            os.system('./4d_points_on_sphere.o ' + number + " " + nsteps + " > " + outputFile + '.out')
            spherePoints = numpy.genfromtxt(outputFile + '.out',skip_header=3, skip_footer=0, usecols=(4,5,6,7), comments='}')
            self.pointDict[number] = spherePoints

        return spherePoints




    def is_overlapping(self, otherMap):
        return self.origin == otherMap.getOrigin() and self.borders == otherMap.getBorders() and self.reso == otherMap.getReso()

    def snap_to_grid(self, otherMap):
        '''Snaps this feature map to its closest equivalent on another grid of equal resolution'''
        nearestOtherIndex = otherMap.point_to_nearest_index(self.getOrigin())
        nearestOtherPoint = otherMap.index_to_coord(nearestOtherIndex[0], nearestOtherIndex[1], nearestOtherIndex[2])
        offset =self.getOrigin().vector_to_new(nearestOtherPoint)
        #Shift the origin of this map to its nearest grid point on the other map
        self.origin.vector_sub_inplace(offset.scalar_mult_new(-1.0))
        #print 'self.origin, otherMap.origin', self.origin, otherMap.origin
        #print 'self.borders, otherMap.borders', self.borders, otherMap.borders
        #print 'nearestOtherIndex', nearestOtherIndex
        #print 'nearestOtherPoint', nearestOtherPoint
        #print 'offset', offset
        #print 'old borders', self.borders
        self.borders = [self.borders[0]+offset[0],
                        self.borders[1]+offset[0],
                        self.borders[2]+offset[1],
                        self.borders[3]+offset[1],
                        self.borders[4]+offset[2],
                        self.borders[5]+offset[2]]

        #print 'newBorders', self.borders


    def translate_inplace(self, vector):
        #for feature in self.features.keys():
        self.borders = [self.borders[0]+vector[0],
                        self.borders[1]+vector[0],
                        self.borders[2]+vector[1],
                        self.borders[3]+vector[1],
                        self.borders[4]+vector[2],
                        self.borders[5]+vector[2]]
        self.origin = point([self.origin[0]+vector[0],
                             self.origin[1]+vector[1],
                             self.origin[2]+vector[2]])


    def interpolation_rotate_inplace(self, quaternion):
        [w,x,y,z] = quaternion
        zphi = 57.29578 * numpy.arctan2((w*y + x*z),-(x*y - w*z))
        xtheta = 57.29578 * numpy.arccos(-w**2 - x**2 + y**2 + z**2)
        zpsi = 57.29578 * numpy.arctan2((w*y - x*z),(x*y + w*z))
        #rotatedArrays = copy.deepcopy(ligandMaps)

        for feature in list(rotatedArrays.keys()):
            self.data = sni.rotate(self.data,zphi,axes=(0,1), reshape=False)
            self.data = sni.rotate(self.data,xtheta,axes=(1,2), reshape=False)
            self.data = sni.rotate(self.data,zpsi,axes=(0,1), reshape=False)

    def to_pdb_string(self, isovalue = 0.):
        relPoints = numpy.transpose(numpy.nonzero(self.data > isovalue))
        absPoints = numpy.array([[(self.reso*x)+self.origin[0],
                                  (self.reso*y)+self.origin[1],
                                  (self.reso*z)+self.origin[2]] for x,y,z in relPoints])



        pdbText = numpy_to_pdb(absPoints,'X')
        return pdbText


    def write_pdb(self, fileName, isovalue = 0.):
        pdbText = self.to_pdb_string(isovalue = isovalue)
        fo = open(fileName,'w')
        fo.write(pdbText)
        fo.close()




    def write_dx_file(self, filename): ## Adapted from POVME
        '''
        Generates a DX file

        Arguments:
        freq_mat -- a Nx4 matrix, where the first 3 columns are the x,y,z coords of the point, and the 4th column is the frequency of emptiness for that point in space

        '''

        header_template = """# Data from POVME 2.0
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
        #N = freq_mat.shape[0] # number of data points

        minx = self.borders[0]
        maxx = self.borders[1]
        miny = self.borders[2]
        maxy = self.borders[3]
        minz = self.borders[4]
        maxz = self.borders[5]

        #minx = min(freq_mat[:,0])
        #miny = min(freq_mat[:,1])
        #minz = min(freq_mat[:,2]) # find the upper and lower corners of the grid
        #maxx = max(freq_mat[:,0])
        #maxy = max(freq_mat[:,1])
        #maxz = max(freq_mat[:,2])

        ##widthx = maxx - minx # find the widths of the grid
        ##widthy = maxy - miny
        ##widthz = maxz - minz

        ##xs = numpy.arange(minx, maxx, self.reso)
        ##ys = numpy.arange(miny, maxy, self.reso)
        ##zs = numpy.arange(minz, maxz, self.reso)

        #xs = numpy.unique(freq_mat[:,0])
        #ys = numpy.unique(freq_mat[:,1])
        #zs = numpy.unique(freq_mat[:,2])

        resx = self.reso
        resy = self.reso
        resz = self.reso
        #resx = xs[1]- xs[0]
        #resy = ys[1]- ys[0]
        #resz = zs[1]- zs[0]

        #resx = freq_mat[(widthz+1)*(widthy+1),0] - freq_mat[0,0]
        #resy = freq_mat[widthz+1,1] - freq_mat[0,1] # find the resolution of the grid
        #resz = freq_mat[1,2] - freq_mat[0,2]

        nx = self.shape[0]
        ny = self.shape[1]
        nz = self.shape[2]
        #nx = int(numpy.round((widthx) / resx + 1)) # number of grid points in each dimension
        #ny = int(numpy.round((widthy) / resy + 1)) # need to add one because the subtraction leaves out an entire row
        #nz = int(numpy.round((widthz) / resz + 1))

        N = int(nx * ny * nz)
        # test to make sure all is well with the size of the grid and its dimensions
        assert (nx * ny * nz) == N, "Something is wrong with the freq_mat array: it is not a prismatic shape"

        # 3. write the header and footer
        #dx_file = gzip.open(filename,'wb')
        dx_file = open(filename,'wb')

        header = header_template % (nx, ny, nz, minx, miny, minz, resx, resy, resz, nx, ny, nz, N) # format the header
        footer = footer_template # the footer needs no formatting
        dx_file.write(header)
        newline_counter = 1
        #for i in range(N): # write the data to the DX file
        for x in range(nx):
            for y in range(ny):
                for z in range(nz):
                    #print self.data[x,y,z]
                    dx_file.write("%8.6e" % self.data[x,y,z])
                    if newline_counter == 3:
                        newline_counter = 0
                        dx_file.write("\n")
                    else:
                        dx_file.write(" ")
                    newline_counter += 1
        dx_file.write(footer)
        dx_file.close
        return

    def index_to_coord(self,indexX, indexY, indexZ):
        x = self.origin[0] + (self.reso * indexX)
        y = self.origin[1] + (self.reso * indexY)
        z = self.origin[2] + (self.reso * indexZ)
        return (x,y,z)

    def point_to_nearest_index(self, thisPoint):
        relCoord = self.origin.vector_to_new(thisPoint)
        #closestX = int(numpy.round(relCoord[0]/self.reso))
        #closestY = int(numpy.round(relCoord[1]/self.reso))
        #closestZ = int(numpy.round(relCoord[2]/self.reso))
        #Avoiding use of numpy.round because it's super slow
        closestX = int(0.5+(relCoord[0]/self.reso))
        closestY = int(0.5+(relCoord[1]/self.reso))
        closestZ = int(0.5+(relCoord[2]/self.reso))
        return (closestX, closestY, closestZ)



    def points_near(self, target, cutoff=5, returnDistance=False, returnVector = False):
        #relTarget = self.origin.vector_to_new(target) #This is already corrected for in point_to_nearest_index

        #print 'ZZZZ', target.coords()
        centerIndex = self.point_to_nearest_index(target)
        #print 'ZZZZA', centerIndex
        indexSearchDist = int(numpy.round(float(cutoff) / self.reso) + 1)
        #indexSearchDist = (cutoff / self.reso) + 1

        minXInd = max(0,numpy.round(centerIndex[0] - indexSearchDist))
        maxXInd = min(self.shape[0],numpy.round(centerIndex[0] + indexSearchDist))
        minYInd = max(0,numpy.round(centerIndex[1] - indexSearchDist))
        maxYInd = min(self.shape[1],numpy.round(centerIndex[1] + indexSearchDist))
        minZInd = max(0,numpy.round(centerIndex[2] - indexSearchDist))
        maxZInd = min(self.shape[2],numpy.round(centerIndex[2] + indexSearchDist))
        pointList = []
        for xInd in range(minXInd, maxXInd):
            for yInd in range(minYInd, maxYInd):
                for zInd in range(minZInd, maxZInd):
                    this_coord = self.index_to_coord(xInd, yInd, zInd)
                    this_point = point(this_coord)
                    distance = self.mf.distance(this_point, target)
                    #print 'AAA', this_coord, target
                    if distance <= cutoff and xInd >= minXInd and xInd <= maxXInd and yInd >= minYInd and yInd <= maxYInd and zInd >= minZInd and zInd <= maxZInd:
                        toAppend = [xInd, yInd, zInd]
                        if returnDistance == True:
                            toAppend.append(distance)
                        if returnVector == True:
                            toAppend.append(target.vector_to_new(this_point))
                        pointList.append(toAppend)

        return pointList


    def add_sphere(self, sphereOrigin, sphereRadius, weight = 1.0):
        #print sphereOrigin
        pointList = self.points_near(sphereOrigin, sphereRadius)
        for x, y, z in pointList:
            self.data[x, y, z] = self.data[x, y, z] + weight
        #print 'AS',numpy.nonzero(self.data)

    def add_sphere_dist_function(self, sphereOrigin, sphereRadius, weightFunction,normalize = True):
        #print sphereOrigin
        pointList = self.points_near(sphereOrigin, sphereRadius, returnDistance=True)
        pointList = numpy.array(pointList)
        #print pointList
        #print weightFunction,
        #print pointList
        #print pointList[:,3]
        if len(pointList) > 0:
            weights = numpy.vectorize(weightFunction)(pointList[:,3])
            if normalize:
                weights = numpy.array(weights) / sum(weights)
        else:
            weights = []
        for i, (x, y, z, d) in enumerate(pointList):
            self.data[x, y, z] = self.data[x, y, z] + weights[i]
            #self.data[x, y, z] = self.data[x, y, z] + weightFunction(d)
        #print 'AS',numpy.nonzero(self.data)


    def prepareGaussianLookupTable(self, n=301):
        if self.gaussianLookupTable == None:
            self.gaussianLookupTable = []
            for i in range(n):
                self.gaussianLookupTable.append(numpy.exp(-((float(i)/100)**2)))

    def add_spherical_gaussian(self, origin, sigmaSq, mu = 0., weight = 1.0,normalize = True):
        #print sphereOrigin
        pointList = self.points_near(origin, sigmaSq, returnDistance=True)
        pointList = numpy.array(pointList)
        self.prepareGaussianLookupTable()


        if len(pointList) > 0:

            weights = [weight * self.gaussianLookupTable[int(100*(i-mu)/(2*sigmaSq))] for i in pointList[:,3]]
            if normalize:
                weights = numpy.array(weights) / sum(weights)
        else:
            weights = []
        for i, (x, y, z, d) in enumerate(pointList):
            self.data[x, y, z] = self.data[x, y, z] + weights[i]
            #self.data[x, y, z] = self.data[x, y, z] + weightFunction(d)
        #print 'AS',numpy.nonzero(self.data)


    def add_cone(self, coneOrigin, coneHalfAngle, coneVector, coneHeight, weight = 1.0):
        #Cone half angle should be in degrees - Here we convert to radians
        coneHalfAngle = numpy.deg2rad(coneHalfAngle)
        coneRadius = coneHeight * numpy.arctan(coneHalfAngle)
        coneVectorNorm = coneVector.scalar_mult_new(1./coneVector.magnitude())


        nearbyPointList = self.points_near(coneOrigin, coneHeight/numpy.cos(coneHalfAngle))
        for x, y, z in nearbyPointList:
            thisCoordAbs = point(self.index_to_coord(x, y, z))
            thisCoordRel = coneOrigin.vector_to_new(thisCoordAbs)
            #height check
            pointDepth = self.mf.dot_product(coneVectorNorm, thisCoordRel)
            if pointDepth < 0 or pointDepth > coneHeight:
                continue
            #Get the position of the point's perpendicular shadow on the axis of the cone
            pointShadowRel = coneVectorNorm.scalar_mult_new(pointDepth)
            radialDistance = pointShadowRel.dist_to(thisCoordRel)
            localConeRadius = coneRadius * (pointDepth/coneHeight)
            if radialDistance <= localConeRadius:
                self.data[x,y,z] += weight


    def add_cone_dist_function(self, coneOrigin, coneHalfAngle, coneVector, coneHeight, weightFunction, normalize = True):
        #Cone half angle should be in degrees - Here we convert to radians
        coneHalfAngle = numpy.deg2rad(coneHalfAngle)
        coneRadius = coneHeight * numpy.arctan(coneHalfAngle)
        coneVectorNorm = coneVector.scalar_mult_new(1./coneVector.magnitude())

        nearbyPointList = self.points_near(coneOrigin, coneHeight/numpy.cos(coneHalfAngle))
        weights = []
        pointsToColor = []
        for x, y, z in nearbyPointList:
            thisCoordAbs = point(self.index_to_coord(x, y, z))
            thisCoordRel = coneOrigin.vector_to_new(thisCoordAbs)
            #height check
            pointDepth = self.mf.dot_product(coneVectorNorm, thisCoordRel)
            if pointDepth < 0 or pointDepth > coneHeight:
                continue
            #Get the position of the point's perpendicular shadow on the axis of the cone
            pointShadowRel = coneVectorNorm.scalar_mult_new(pointDepth)
            radialDistance = pointShadowRel.dist_to(thisCoordRel)
            localConeRadius = coneRadius * (pointDepth/coneHeight)
            if radialDistance <= localConeRadius:
                weights.append(weightFunction(thisCoordRel.magnitude()))
                pointsToColor.append([x,y,z])
        if normalize:
            weights = numpy.array(weights) / sum(weights)
        for (x, y, z), weight in zip(pointsToColor,weights):
            self.data[x,y,z] += weight


    def add_disc(self, discOrigin, discHeight, discNormal, discInnerRadius, discOuterRadius, bidirectional=True, weight = 1):
        nearbyPointList = self.points_near(discOrigin, point([discOuterRadius,discHeight,0.]).magnitude())
        discNormal.scalar_mult_inplace(1.0/discNormal.magnitude())
        discVector = discNormal.scalar_mult_new(discHeight)
        if bidirectional:
            self.add_disc(discOrigin, discHeight, discNormal.scalar_mult_new(-1.0), discInnerRadius,
                          discOuterRadius, bidirectional=False, weight = weight)
            #discOrigin = discOrigin.vector_sub_new(discVector)
            #discVector.scalar_mult_inplace(2.0)
            #discHeight *= 2.

        for x,y,z in nearbyPointList:
            thisCoordAbs = point(self.index_to_coord(x,y,z))
            thisCoordRel = discOrigin.vector_to_new(thisCoordAbs)
            #Find projection on normal vector and ensure that point falls inside height
            pointDepth = self.mf.dot_product(discNormal, thisCoordRel)
            if pointDepth > discHeight or pointDepth < 0:
                continue

            pointShadowRel = discNormal.scalar_mult_new(pointDepth)

            #Get distance from point to its perpendicular shadow on the normal and
            #ensure that it's between inner and outer radius
            radialDistance = pointShadowRel.dist_to(thisCoordRel)
            if radialDistance < discOuterRadius and radialDistance > discInnerRadius:
                self.data[x,y,z] += weight

    def add_disc_dist_function(self, discOrigin, discHeight, discNormal, discInnerRadius, discOuterRadius, weightFunction, bidirectional=True):
        nearbyPointList = self.points_near(discOrigin, point([discOuterRadius,discHeight,0.]).magnitude())
        discNormal.scalar_mult_inplace(1.0/discNormal.magnitude())
        #discVector = discNormal.scalar_mult_new(discHeight)
        if bidirectional:
            self.add_disc_dist_function(discOrigin, discHeight, discNormal.scalar_mult_new(-1.0),
                                        discInnerRadius, discOuterRadius, weightFunction, bidirectional = False)
            #discOrigin = discOrigin.vector_sub_new(discVector)
            #discVector.scalar_mult_inplace(2.0)
            #discHeight *= 2.

        pointsToColor = []
        weights = []
        for x,y,z in nearbyPointList:
            thisCoordAbs = point(self.index_to_coord(x,y,z))
            thisCoordRel = discOrigin.vector_to_new(thisCoordAbs)
            #Find projection on normal vector and ensure that point falls inside height
            pointDepth = self.mf.dot_product(discNormal, thisCoordRel)
            if pointDepth > discHeight or pointDepth < 0:
                continue

            pointShadowRel = discNormal.scalar_mult_new(pointDepth)

            #Get distance from point to its perpendicular shadow on the normal and
            #ensure that it's between inner and outer radius
            radialDistance = pointShadowRel.dist_to(thisCoordRel)
            if radialDistance < discOuterRadius and radialDistance > discInnerRadius:
                pointsToColor.append([x,y,z])
                weights.append(weightFunction(pointDepth))

        weights = numpy.array(weights) / sum(weights)
        for (x,y,z), weight in zip(pointsToColor,weights):
            self.data[x,y,z] += weight





class algebra:
    def __init__(self):
        self.stringScoreFuncs = ['hydrophobic_A * hydrophobic_B',
                                 'hbondAcceptor_A * hbondDonor_B',
                                 'hbondDonor_A * hbondAcceptor_B',
                                 'aromatic_A * aromatic_B',
                                 'hydrophobicity_A * hydrophobicity_B']
        self.scoreFuncs = []
        self.vecScoreFuncs = []
        self.prepareScoreFuncs()
        self.lastScoreMaps = None
        self.lastScores = None
        self.lastScoreFuncs = None

    def getScoreFuncs(self):
        return self.stringScoreFuncs

    def setScoreFuncs(self, scoreFuncList):
        self.stringScoreFuncs = scoreFuncList
        self.prepareScoreFuncs()


    def prepareScoreFuncs(self):
        self.scoreFuncs = []
        self.vecScoreFuncs = []
        for stringScoreFunc in self.stringScoreFuncs:
            standardizedStringScoreFunc = stringScoreFunc
            features = []
            variables = re.findall('([a-zA-Z0-9]+)_([AB])',stringScoreFunc)
            #print variables
            c=0
            standardMatNames = []
            for variable in variables:
                standardMatName = 'M%i' %(c)
                standardMatNames.append(standardMatName)
                ### So the thing appended would be like ['M1','hydrophobic','A']
                features.append([standardMatName, variable[0], variable[1]])

                standardizedStringScoreFunc = standardizedStringScoreFunc.replace('_'.join(variable),standardMatName)
                c+=1

            # If it's just two matrices being multiplied then we can use numpy.multiply, which will be very fast
            operators = stringScoreFunc.replace(variables[0][0],'').replace(variables[1][0],'').replace(variables[0][1],'').replace(variables[1][1],'').replace('_','').strip()
            #print operators
            if len(variables) == 2 and operators=='*':
                self.scoreFuncs.append([features, numpy.multiply])
                #self.vecScoreFuncs.append([[standardMatName, variable[0], variable[1]],numpy.multiply])
                self.vecScoreFuncs.append([features,numpy.multiply])

                #print 'SIMPLIFIED SCORE FUNCTION WILL BE USED FOR:', stringScoreFunc
            else:
                scoreFunc =  eval('lambda %s: %s' %(','.join(standardMatNames), standardizedStringScoreFunc))
                #print 'lambda %s: %s' %(','.join(standardMatNames), standardizedStringScoreFunc)
                self.scoreFuncs.append([features, scoreFunc])
                #self.vecScoreFuncs.append([[standardMatName, variable[0], variable[1]],numpy.vectorize(scoreFunc)])
                self.vecScoreFuncs.append([features,numpy.vectorize(scoreFunc)])
        #self.vecScoreFuncs = [[feature, numpy.vectorize(func)] for feature, func in self.scoreFuncs]

    def scoreAll(self,A, B):
        scores = []
        scoreMaps = []

        for features, vecScoreFunc in self.vecScoreFuncs:
            args = []
            ### Each feature is of the format ['M1','hydrophobic','A']
            for feature in features:
                if feature[2] == 'A':
                    args.append(A[feature[1]])
                if feature[2] == 'B':
                    args.append(B[feature[1]])
            scoreMap = self.scoreOne(args, vecScoreFunc)
            #scoreMap = self.scoreOne(A[features[0]], B[features[1]], vecScoreFunc)
            scoreMaps.append(scoreMap)
            volumeNormFactor = pow(A[feature[1]].getReso(), 3)
            termScore = volumeNormFactor * numpy.sum(scoreMap.getData().flat)
            scores.append(termScore)
        self.lastScoreMaps = scoreMaps
        self.lastScores = scores
        self.lastStringScoreFuncs = self.stringScoreFuncs
        return scoreMaps, scores

    def scoreOne(self, args, vecScoreFunc):
        #print A, B
        resultMap = None
        #if not(A.isInSameLocation(B)):
            #A,B = self.commonizeVolume(A,B)
        #for arg in args:
        #    print 'raw origin', arg.getOrigin(), 'shape', arg.getShape(), 'borders',  arg.getBorders()
        commonizedArgs = self.commonizeVolume(args)
        #for arg in commonizedArgs:
        #    print 'commonized', arg.getOrigin(), 'shape', arg.getShape(), 'borders',  arg.getBorders()
        #commonizedArgs = args

        #print commonizedArgs
        resultMatrix = vecScoreFunc(*[i.getData() for i in commonizedArgs])
        resultMap = featureMap(commonizedArgs[0].getBorders(), commonizedArgs[0].getReso(), boolean=False)
        resultMap.setData(resultMatrix)
        return resultMap




    def commonizeVolume(self, fMaps):
        #Look at *smallest* common area
        #bottomLeft = point([max(fMap1.getBorders()[0], fMap2.getBorders()[0]),
        #                    max(fMap1.getBorders()[2], fMap2.getBorders()[2]),
        #                    max(fMap1.getBorders()[4], fMap2.getBorders()[4])])
        #topRight  =  point([min(fMap1.getBorders()[1], fMap2.getBorders()[1]),
        #                    min(fMap1.getBorders()[3], fMap2.getBorders()[3]),
        #                    min(fMap1.getBorders()[5], fMap2.getBorders()[5])])


        #Ensure that maps have the same resolution
        for fMap in fMaps[1:]:
            if fMap.getReso() != fMaps[0].getReso():
                raise Exception("Two grids entered into commonizeVolume don't have the same resolution")
            fMap.snap_to_grid(fMaps[0])

        bottomLeft = point([max([i.getBorders()[0] for i in fMaps]),
                            max([i.getBorders()[2] for i in fMaps]),
                            max([i.getBorders()[4] for i in fMaps])])
        topRight  =  point([min([i.getBorders()[1] for i in fMaps]),
                            min([i.getBorders()[3] for i in fMaps]),
                            min([i.getBorders()[5] for i in fMaps])])
        #print 'In commonizeVolume: [i.getBorders() for i in fMaps]', [i.getBorders() for i in fMaps], ' bottomLeft', bottomLeft, 'topRight',topRight
        #fMap1BottomLeft = fMap1.getOrigin()
        #fMap1TopRight = point([fMap1.getBorders()[1],fMap1.getBorders()[3],fMap1.getBorders()[5]])
        #### I SHOULDN'T BE ROUNDING - I SHOULD TRANSLATE TO PUT THEM ON THE SAME GRID
        #fMap1BLShift = fMap1BottomLeft.vector_to_new(bottomLeft)
        #fMap1TRShift = fMap1TopRight.vector_to_new(topRight)

        fmapBottomLefts = [i.getOrigin() for i in fMaps]
        fmapTopRights = [point([i.getBorders()[1],i.getBorders()[3],i.getBorders()[5]]) for i in fMaps]
        #### I SHOULDN'T BE ROUNDING - I SHOULD TRANSLATE TO PUT THEM ON THE SAME GRID
        BLShiftsAngstroms = [fmapBottomLeft.vector_to_new(bottomLeft) for fmapBottomLeft in fmapBottomLefts]
        TRShiftsAngstroms = [fmapTopRight.vector_to_new(topRight) for fmapTopRight in fmapTopRights]

        BLShiftsReso = [i.scalar_mult_new(1./fMaps[0].getReso()) for i in BLShiftsAngstroms]
        TRShiftsReso = [i.scalar_mult_new(1./fMaps[0].getReso()) for i in TRShiftsAngstroms]

        #fMap2BottomLeft = fMap2.getOrigin()
        #fMap2TopRight = point([fMap2.getBorders()[1],fMap2.getBorders()[3],fMap2.getBorders()[5]])
        #fMap2BLShift = fMap2BottomLeft.vector_to_new(bottomLeft)
        #fMap2TRShift = fMap2TopRight.vector_to_new(topRight)

        #fMap2Shift = bottomLeft.vector_to_new(fMap2.getOrigin())
        #### SCALE SHIFTS TO RESOLUTION AND IMPLEMENT SOLUTION TO OFF-GRID SHIFTS
        # LOOK AT point.scalar_mult_inplace

        newFMaps = [featureMap([bottomLeft[0], topRight[0],
                                bottomLeft[1], topRight[1],
                                bottomLeft[2], topRight[2]],
                               fMaps[i].getReso()) for i in range(len(fMaps))]
        #print 'In commonizeVolume [i.getShape() for i in newFMaps]', [i.getShape() for i in newFMaps]
        #print 'In commonizeVolume [i.getData().shape for i in newFMaps]', [i.getData().shape for i in newFMaps]
        ##newFMap1 = featureMap([bottomLeft[0], topRight[0],
        ##                       bottomLeft[1], topRight[1],
        ##                       bottomLeft[2], topRight[2]],
        ##                      fMap1.getReso())
        ##newFMap2 = featureMap([bottomLeft[0], topRight[0],
        ##                       bottomLeft[1], topRight[1],
        ##                       bottomLeft[2], topRight[2]],
        ##                      fMap2.getReso())
        #newData1 = [newFMap.getData() for newFMap in newFMaps]
        ##print 'SHAPES', fMap1.getData().shape, fMap2.getData().shape
        ##print fMap1BLShift.coords()
        ##print fMap1TRShift.coords()
        ##print fMap1.getData().shape
        ##print [fMap1BLShift[0], fMap1.getShape()[0]+fMap1TRShift[0], fMap1BLShift[1], fMap1.getShape()[1]+fMap1TRShift[1], fMap1BLShift[2], fMap1.getShape()[2]+fMap1TRShift[2]]
        for i in range(len(fMaps)):
            #print 'i, fMaps[i].getShape, BLShifts[i], TRShifts[i]'
            #print i, fMaps[i].getShape(), BLShifts[i], TRShifts[i]

            #newData = fMaps[i].getData()[numpy.round(BLShiftsReso[i][0]):numpy.round(fMaps[i].getShape()[0]+TRShiftsReso[i][0]),
            #                             numpy.round(BLShiftsReso[i][1]):numpy.round(fMaps[i].getShape()[1]+TRShiftsReso[i][1]),
            #                             numpy.round(BLShiftsReso[i][2]):numpy.round(fMaps[i].getShape()[2]+TRShiftsReso[i][2])]
            # Not using numpy.round because it's incredibly slow
            newData = fMaps[i].getData()[int(BLShiftsReso[i][0]+0.5):int(fMaps[i].getShape()[0]+TRShiftsReso[i][0]+0.5),
                                         int(BLShiftsReso[i][1]+0.5):int(fMaps[i].getShape()[1]+TRShiftsReso[i][1]+0.5),
                                         int(BLShiftsReso[i][2]+0.5):int(fMaps[i].getShape()[2]+TRShiftsReso[i][2]+0.5)]


            newFMaps[i].setData(newData)
            #print 'in commonizeVolumes newData.shape:', newData.shape, 'newFMaps[i].getShape()', newFMaps[i].getShape()
            #print newFMaps[i].getData().shape
        #newData2 = newFMap2.getData()
        #newData2= fMap2.getData()[numpy.round(fMap2BLShift[0]):numpy.round(fMap2.getShape()[0]+fMap2TRShift[0]),
        #                          numpy.round(fMap2BLShift[1]):numpy.round(fMap2.getShape()[1]+fMap2TRShift[1]),
        #                          numpy.round(fMap2BLShift[2]):numpy.round(fMap2.getShape()[2]+fMap2TRShift[2])]
        #newFMap2.setData(newData2)
        # newData1 = newFMap1.getData()
        # newData1[fMap1Shift[0]:fMap1Shift[0]+fMap1.getShape()[0],
        #          fMap1Shift[1]:fMap1Shift[1]+fMap1.getShape()[1],
        #          fMap1Shift[2]:fMap1Shift[2]+fMap1.getShape()[2]] = fMap1.getData()
        # newFMap1.setData(newData1)
        # newData2 = newFMap2.getData()
        # newData2[fMap2Shift[0]:fMap2Shift[0]+fMap2.getShape()[0],
        #          fMap2Shift[1]:fMap2Shift[1]+fMap2.getShape()[1],
        #          fMap2Shift[2]:fMap2Shift[2]+fMap2.getShape()[2]] = fMap2.getData()
        # newFMap2.setData(newData2)
        return newFMaps





    def dockOne(self,receptorMaps, ligandMaps,translation_list,spherePoints):

        #successfulDockings = numpy.array(['Translation Steps','Rotation Quaternion','Score'],dtype='string_')
        successfulDockings=[]
        c=0
        for [w,x,y,z] in spherePoints:
            c += 1
            print("Processing rotation %i of %i" %(c, len(spherePoints)))
            #Quaternions to Euler angles:
            #180/3.1415926
            #zphi = 57.29578 * numpy.arctan2((w*y + x*z),-(x*y - w*z))
            #xtheta = 57.29578 * numpy.arccos(-w**2 - x**2 + y**2 + z**2)
            #zpsi = 57.29578 * numpy.arctan2((w*y - x*z),(x*y + w*z))
            rotatedArrays = copy.deepcopy(ligandMaps)
            for feature in list(rotatedArrays.keys()):
                rotatedArrays[feature].interpolation_rotate_inplace([w,x,y,z])
            #print zphi, xtheta, zpsi
            #for feature in rotatedArrays.keys():
            #    rotatedArrays[feature].data = sni.rotate(rotatedArrays[feature].data,zphi,axes=(0,1), reshape=False)
            #    rotatedArrays[feature].data = sni.rotate(rotatedArrays[feature].data,xtheta,axes=(1,2), reshape=False)
            #    rotatedArrays[feature].data = sni.rotate(rotatedArrays[feature].data,zpsi,axes=(0,1), reshape=False)





            for [X,Y,Z] in translation_list:
                #shiftedArray = sni.shift(receptorMap,[X,Y,Z])

                shiftedArrays = copy.deepcopy(rotatedArrays)
                for feature in list(shiftedArrays.keys()):
                    shiftedArrays[feature].translate_inplace([X,Y,Z])
                    #shiftedArrays[feature].borders = [shiftedArrays[feature].borders[0]+X,
                    #                                  shiftedArrays[feature].borders[1]+X,
                    #                                  shiftedArrays[feature].borders[2]+Y,
                    #                                  shiftedArrays[feature].borders[3]+Y,
                    #                                  shiftedArrays[feature].borders[4]+Z,
                    #                                  shiftedArrays[feature].borders[5]+Z]
                    #shiftedArrays[feature].origin = point([shiftedArrays[feature].origin[0]+X,
                    #                                       shiftedArrays[feature].origin[1]+Y,
                    #                                       shiftedArrays[feature].origin[2]+Z])








                scoreMaps, scores = self.scoreAll(receptorMaps, shiftedArrays)
                #totalScore = sum(scores)
                #if rotate3.score.rank in range(50):
                successfulDockings.append([[X,Y,Z],[w,x,y,z],scores])
                #print [[X,Y,Z],[w,x,y,z],scores]
        successfulDockings.sort(key=lambda x: sum(x[2]), reverse = True)
        print("Docking complete!")
        return successfulDockings

    def dockPeel(self,receptorMaps, ligandPeel, translation_list,spherePoints):

        resolution = receptorMaps[list(receptorMaps.keys())[0]].getReso()

        #successfulDockings = numpy.array(['Translation Steps','Rotation Quaternion','Score'],dtype='string_')
        successfulDockings=[]
        c=0
        for [w,x,y,z] in spherePoints:
            c += 1
            print("Processing rotation %i of %i: %s" %(c, len(spherePoints), str([w,x,y,z])))
            #Quaternions to Euler angles:
            #180/3.1415926
            #zphi = 57.29578 * numpy.arctan2((w*y + x*z),-(x*y - w*z))
            #xtheta = 57.29578 * numpy.arccos(-w**2 - x**2 + y**2 + z**2)
            #zpsi = 57.29578 * numpy.arctan2((w*y - x*z),(x*y + w*z))
            rotatedPeel = copy.deepcopy(ligandPeel)

            rotatedPeel.rotate([w,x,y,z])
            #borders = numpy.array([99999.,-99999.,99999.,-99999.,99999.,-99999.])

            #for feature in rotatedPeel.features.keys():
            #    for index in rotatedPeel.features[feature].keys():
            #        oldCoordsList = rotatedPeel.features[feature][index]['coordinates'].coords()
            #        newCoordsList = rotatedPeel.functions.qv_mult((w,x,y,z), oldCoordsList)
            #        #coordsList = [coords.x, coords.y, coords.z]
            #        rotatedPeel.features[feature][index]['coordinates'] = point(newCoordsList)
            #        for key in rotatedPeel.features[feature][index]:
            #            if 'vector' in key:
            #                coordsList = rotatedPeel.features[feature][index][key].coords()
            #                #coordsList = [coords.x, coords.y, coords.z]
            #                rotatedPeel.features[feature][index][key] = point(rotatedPeel.functions.qv_mult((w,x,y,z), coordsList))

                    #if newCoordsList[0] < borders[0]:
                    #    borders[0] = newCoordsList[0]
                    #if newCoordsList[0] > borders[1]:
                    #    borders[1] = newCoordsList[0]
                    #if newCoordsList[1] < borders[2]:
                    #    borders[2] = newCoordsList[1]
                    #if newCoordsList[1] > borders[3]:
                    #    borders[3] = newCoordsList[1]
                    #if newCoordsList[2] < borders[4]:
                    #    borders[4] = newCoordsList[2]
                    #if newCoordsList[2] > borders[5]:
                    #    borders[5] = newCoordsList[2]


            #borders += numpy.array([-8., 8., -8., 8., -8., 8.])\
            borders = rotatedPeel.suggest_borders(padding = 6.0)

            rotatedMaps = rotatedPeel.create_feature_maps(borders, resolution)


            lastTranslationReso = [0.,0.,0.]
            shiftedMaps = copy.deepcopy(rotatedMaps)
            for [absXReso,absYReso,absZReso] in translation_list:
            #for [X,Y,Z] in translation_list:
                #shiftedArray = sni.shift(receptorMap,[X,Y,Z])
                #for feature in shiftedPeel.features.keys():
                    #for index in shiftedPeel.features[feature].keys():
                        #coords = shiftedPeel.features[feature][index]['coordinates'].coords()

                        #shiftedPeel.features[feature][index]['coordinates'] = point(coords + numpy.array([X,Y,Z]))
                dXReso = absXReso - lastTranslationReso[0]
                dYReso = absYReso - lastTranslationReso[1]
                dZReso = absZReso - lastTranslationReso[2]

                absXAngstroms = absXReso * resolution
                absYAngstroms = absYReso * resolution
                absZAngstroms = absZReso * resolution

                lastTranslationReso[0] = absXReso
                lastTranslationReso[1] = absYReso
                lastTranslationReso[2] = absZReso

                dXAngstroms = dXReso * resolution
                dYAngstroms = dYReso * resolution
                dZAngstroms = dZReso * resolution

                #shiftedMaps = copy.deepcopy(rotatedMaps)

                for feature in list(shiftedMaps.keys()):
                    shiftedMaps[feature].translate_inplace([dXAngstroms,dYAngstroms,dZAngstroms])
                    #shiftedMaps[feature].translate_inplace([X,Y,Z])
                    #shiftedMaps[feature].borders = [shiftedMaps[feature].borders[0]+X,
                    #                                shiftedMaps[feature].borders[1]+X,
                    #                                shiftedMaps[feature].borders[2]+Y,
                    #                                shiftedMaps[feature].borders[3]+Y,
                    #                                shiftedMaps[feature].borders[4]+Z,
                    #                                shiftedMaps[feature].borders[5]+Z]
                    #shiftedMaps[feature].origin = point([shiftedMaps[feature].origin[0]+X,
                    #                                     shiftedMaps[feature].origin[1]+Y,
                    #                                     shiftedMaps[feature].origin[2]+Z])

                #if X==0 and Y==0 and Z==0:
                #    shiftedMaps['aromatic'].write_dx_file('aromatic_%.3e_%.3e_%.3e_%.3e_%.3e_%.3e_%.3e.dx' %(w,x,y,z,X,Y,Z))
                #    shiftedMaps['occupancy'].write_dx_file('occ_%.3e_%.3e_%.3e_%.3e_%.3e_%.3e_%.3e.dx' %(w,x,y,z,X,Y,Z))

                scoreMaps, scores = self.scoreAll(receptorMaps, shiftedMaps)

                #totalScore = sum(scores)
                #if rotate3.score.rank in range(50):
                successfulDockings.append([[absXAngstroms,absYAngstroms,absZAngstroms],[w,x,y,z],scores])
                #print [[X,Y,Z],[w,x,y,z],scores]
        successfulDockings.sort(key=lambda x: sum(x[2]), reverse = True)
        print("Docking complete!")
        return successfulDockings

    def printLastScores(self):
        if self.lastScores != None:
            print('\n'.join(['%.5e\t%s'%(score, func) for func, score in zip(self.lastStringScoreFuncs, self.lastScores)]))




class peel:



    # supporting functions
    #def list_alphebetize_and_combine(self, list):
    #    list.sort()
    #    return '_'.join(list)

    def hashtable_entry_add_one(self, hashtable, key, toadd = 1): # note that dictionaries (hashtables) are passed by reference in python
        if key in hashtable:
            hashtable[key] = hashtable[key] + toadd
        else:
            hashtable[key] = toadd

    def center(self, string, length):
        while len(string) < length:
            string = " " + string
            if len(string) < length:
                string = string + " "
        return string

    def bin_atom_dictionary(self):
        pass
        
    def bin_hbond_acceptors(self):
        #bins hbond receptors into neighborhoods so that only nearby ones need to be checked when characterizing grid points
        #self.hba_regions = {}
        #self.hba_region_spacing = 10
        #self.hba_region_skin = 2
        pass

    def bin_hbond_donors(self):
        #bins hbond receptors into neighborhoods so that only nearby ones need to be checked when characterizing grid points
        #self.hba_regions = {}
        #self.hba_region_spacing = 10
        #self.hba_region_skin = 2
        pass

    def bin_hydrophobics(self):
        pass

    def bin_hydrophilics(self):
        pass

    def bin_aromatics(self):
        pass

    def bin_occupancies(self):
        pass


    def color_map_by_aromatic(self, this_feature_map):
        if not(self.aromaticsBinned):
            self.bin_aromatics()
            self.aromaticsBinned = True

        #print 'AROMATICS', self.features['ar_rings']
        for ar_ring in self.features['ar_rings']:

            #this_feature_map.add_disc(self.features['ar_rings'][ar_ring]['coordinates'],
            #                          self.parameters['pi_pi_interacting_dist_cutoff'],
            #                          self.features['ar_rings'][ar_ring]['norm_vector'],
            #                          0.0,
            #                          self.features['ar_rings'][ar_ring]['radius'],
            #                          weight=1.)
            this_feature_map.add_disc_dist_function(self.features['ar_rings'][ar_ring]['coordinates'],
                                                    self.parameters['aromatic_height_cutoff'],
                                                    self.features['ar_rings'][ar_ring]['norm_vector'],
                                                    0.0,
                                                    self.features['ar_rings'][ar_ring]['radius'],
                                                    lambda x: numpy.exp(-numpy.power(x-self.parameters['aromatic_gaussian_mean']/self.parameters['aromatic_gaussian_variance'],2)/2.0))
        #print 'CMBAromatic', numpy.nonzero(this_feature_map.data)

    def color_map_by_atomDictionary(self, this_feature_map, category):
        '''if not(self.atomDictionaryBinned):
            self.bin_atom_dictionary()
            self.atomDictionaryDonorsBinned = True'''
                
        for index in self.atom_dictionary[category]:
            this_feature_map.add_spherical_gaussian(point(self.receptor.information.get_coordinates()[index]),
                                                    self.parameters['default_gaussian_variance'])

    def color_map_by_hbondDonor(self, this_feature_map):
        if not(self.hbondDonorsBinned):
            self.bin_hbond_donors()
            self.hbondDonorsBinned = True

        #print 'HB DONORS', self.features['hb_donors']
        for hbondDonor in self.features['hb_donors']:

            #this_feature_map.add_cone(self.features['hb_donors'][hbondDonor]['coordinates'],
            #                          self.parameters['hydrogen_bond_angle_cutoff'],
            #                          self.features['hb_donors'][hbondDonor]['hbd_vector'],
            #                          self.parameters['hydrogen_bond_dist_cutoff'],
            #                          1)
            this_feature_map.add_cone_dist_function(self.features['hb_donors'][hbondDonor]['coordinates'],
                                                    self.parameters['hydrogen_bond_angle_cutoff'],
                                                    self.features['hb_donors'][hbondDonor]['hbd_vector'],
                                                    self.parameters['hydrogen_bond_dist_cutoff'],
                                                    lambda x: numpy.exp(-numpy.power((x-self.parameters['hydrogen_bond_donor_gaussian_mean'])/self.parameters['hydrogen_bond_donor_gaussian_variance'],2)/2))
        #print 'CMBHBD', numpy.nonzero(this_feature_map.data)

    def color_map_by_hbondAcceptor(self, this_feature_map):
        if not(self.hbondAcceptorsBinned):
            self.bin_hbond_acceptors()
            self.hbondAcceptorsBinned = True


        for hbondAcceptor in self.features['hb_acceptors']:

            #this_feature_map.add_sphere(self.features['hb_acceptors'][hbondAcceptor],
            #                      self.parameters['hydrogen_bond_acceptor_radius'],
            #                      1)
            this_feature_map.add_sphere_dist_function(self.features['hb_acceptors'][hbondAcceptor]['coordinates'],
                                                      self.parameters['hydrogen_bond_acceptor_radius'],
                                                      lambda x: numpy.exp(-numpy.power(x/self.parameters['hydrogen_bond_acceptor_gaussian_variance'],2)/2))
        #print 'CMBHBA', numpy.nonzero(this_feature_map.data)


    def color_map_by_hydrophobic(self, this_feature_map):
        if not(self.hydrophobicsBinned):
            self.bin_hydrophobics()
            self.hydrophobicsBinned = True


        for hydrophobic in self.features['hydrophobics']:

            #this_feature_map.add_sphere(self.features['hydrophobics'][hydrophobic],
            #                      self.parameters['hydroph_dist_cutoff'],
            #                      1)
            #this_feature_map.add_sphere_dist_function(self.features['hydrophobics'][hydrophobic]['coordinates'],
            #                                          self.parameters['hydroph_gaussian_cutoff'],
            #                                          lambda x: numpy.exp(-numpy.power(x/self.parameters['hydroph_gaussian_variance'],2)/2))
            this_feature_map.add_spherical_gaussian(self.features['hydrophobics'][hydrophobic]['coordinates'],
                                          self.parameters['hydroph_gaussian_variance'])
        #print 'CMBHydrophobic', numpy.nonzero(this_feature_map.data)

    def color_map_by_hydrophilic(self, this_feature_map):
        if not(self.hydrophilicsBinned):
            self.bin_hydrophilics()
            self.hydrophilicsBinned = True


        for hydrophilic in self.features['hydrophilics']:

            #this_feature_map.add_sphere(self.features['hydrophilics'][hydrophilic],
            #                      self.parameters['hydroph_dist_cutoff'],
            #                      1)
            #this_feature_map.add_sphere_dist_function(self.features['hydrophilics'][hydrophilic]['coordinates'],
            #                                          self.parameters['hydroph_gaussian_cutoff'],
            #                                          lambda x: numpy.exp(-numpy.power(x/self.parameters['hydroph_gaussian_variance'],2)/2))
            this_feature_map.add_spherical_gaussian(self.features['hydrophilics'][hydrophilic]['coordinates'],
                                          self.parameters['hydroph_gaussian_variance'])
        #print 'CMBHydrophilic', numpy.nonzero(this_feature_map.data)


    def color_map_by_hydrophobicity(self, this_feature_map):

        if not(self.hydrophobicsBinned):
            self.bin_hydrophobics()
            self.hydrophobicsBinned = True

        vectorizedWeightFunc = numpy.vectorize(lambda x: numpy.exp(-x/3), otypes=[numpy.float])

        for hydrophobic in self.features['hydrophobics']:

            this_feature_map.add_sphere_dist_function(self.features['hydrophobics'][hydrophobic]['coordinates'],
                                                      6,
                                                      vectorizedWeightFunc)

        if not(self.hydrophilicsBinned):
            self.bin_hydrophilics()
            self.hydrophilicsBinned = True

        vectorizedWeightFunc = numpy.vectorize(lambda x: -numpy.exp(-x/3), otypes=[numpy.float])


        for hydrophilic in self.features['hydrophilics']:

            this_feature_map.add_sphere_dist_function(self.features['hydrophilics'][hydrophilic]['coordinates'],
                                                      6,
                                                      vectorizedWeightFunc)

        #print 'CMBHydrophobicity', numpy.nonzero(this_feature_map.data)


    def color_map_by_occupancy(self, this_feature_map):
        if not(self.occupanciesBinned):
            self.bin_occupancies()
            self.occupanciesBinned = True
        for this_atom in self.features['occupancies']:
            this_feature_map.add_sphere(self.features['occupancies'][this_atom]['coordinates'],
                                        self.features['occupancies'][this_atom]['radius'],
                                        1)


    def color_map_by_adjacency(self, this_feature_map):
        if not(self.occupanciesBinned):
            self.bin_occupancies()
            self.occupanciesBinned = True
        for this_atom in self.features['occupancies']:
            this_feature_map.add_sphere(self.features['occupancies'][this_atom]['coordinates'],
                                        self.features['occupancies'][this_atom]['radius']*2,
                                        1)





    def create_feature_maps(self, gridBorders, gridReso, features = ['occupancy','adjacency', 'hbondAcceptor','hbondDonor','aromatic','hydrophobic','hydrophilic']):
        # Returns multiple povme maps, each one colored by different ligand features
        # Features to implement: electrostatic, hbondAcceptor, hbondDonor, hbondDonorWDirection, aromatic, aromaticWDirection, hydrophobic
        #featureMaps = [featureMap(gridBorders, gridReso, boolean=False) for i in features]
        featureMaps = {}
        for feature in features:
            if feature == 'occupancy':
                featureMaps[feature] = featureMap(gridBorders, gridReso, boolean=True)
                start = time.time()
                self.color_map_by_occupancy(featureMaps[feature])
                end = time.time()
                print('time to create occupancy map:', end-start)
            if feature == 'adjacency':
                featureMaps[feature] = featureMap(gridBorders, gridReso, boolean=True)
                start = time.time()
                self.color_map_by_adjacency(featureMaps[feature])
                end = time.time()
                print('time to create adjacency map:', end-start)
            if feature == 'hbondAcceptor':
                featureMaps[feature] = featureMap(gridBorders, gridReso, boolean=False)
                start = time.time()
                self.color_map_by_hbondAcceptor(featureMaps[feature])
                end = time.time()
                print('time to create hydrogen bond acceptor map:', end-start)
            if feature == 'hbondDonor':
                featureMaps[feature] = featureMap(gridBorders, gridReso, boolean=False)
                start = time.time()
                self.color_map_by_hbondDonor(featureMaps[feature])
                end = time.time()
                print('time to create hydrogen bond donor map:', end-start)
            if feature == 'aromatic':
                featureMaps[feature] = featureMap(gridBorders, gridReso, boolean=False)
                start = time.time()
                self.color_map_by_aromatic(featureMaps[feature])
                end = time.time()
                print('time to create aromatic map:', end-start)
            if feature == 'hydrophobic':
                featureMaps[feature] = featureMap(gridBorders, gridReso, boolean=False)
                start = time.time()
                self.color_map_by_hydrophobic(featureMaps[feature])
                end = time.time()
                print('time to create hydrophobic map:', end-start)
            if feature == 'hydrophilic':
                featureMaps[feature] = featureMap(gridBorders, gridReso, boolean=False)
                start = time.time()
                self.color_map_by_hydrophilic(featureMaps[feature])
                end = time.time()
                print('time to create hydrophilic map:', end-start)
            if feature == 'hydrophobicity':
                featureMaps[feature] = featureMap(gridBorders, gridReso, boolean=False)
                start = time.time()
                self.color_map_by_hydrophobicity(featureMaps[feature])
                end = time.time()
                print('time to create hydrophobicity map:', end-start)
        #print 'CFM', numpy.nonzero(featureMaps[0].data)
        return featureMaps
    
    def create_atom_map(self, gridBorders, gridReso = 1.0):
        # Returns multiple povme maps by atom categories
        atomMaps = {}
        for category in self.atom_dictionary:
            atomMaps[category] = featureMap(gridBorders, gridReso, boolean=True)
            start = time.time()
            self.color_map_by_atomDictionary(atomMaps[category], category)
            end = time.time()
            print('time to create atom category ' + category + ' map:', end-start)
                
        return atomMaps

    def color_povme_map(self, povmeMatrix, gridReso, features = ['hbondAcceptor','hbondDonor','aromatic','hydrophobic', 'hydrophilic', 'occupancy', 'adjacency'], skin = 0):
        ''' Takes a POVME map (list of points in numpy format) and its resolution.
        Returns a dictionary, with colors as keys and POVME maps for each as values
        NOTE that the order of points in input and output maps are not necessarily the same'''
        #First, get dimensions of the POVME map (remember that it's a list of points)

        minX = numpy.min(povmeMatrix[:,0])
        maxX = numpy.max(povmeMatrix[:,0]) + 1e-5
        minY = numpy.min(povmeMatrix[:,1])
        maxY = numpy.max(povmeMatrix[:,1]) + 1e-5
        minZ = numpy.min(povmeMatrix[:,2])
        maxZ = numpy.max(povmeMatrix[:,2]) + 1e-5
        colorMaps = self.create_feature_maps([minX, maxX, minY, maxY, minZ, maxZ], gridReso, features = features)

        colorPointLists = {}

        povmeList = povmeMatrix.tolist()

        povmeSet = set([ tuple(i) for i in povmeList])
        #print 'BEFORE', len(povmeSet)

        if skin != 0:
            # Create a list of allowable moves which stay within the skin distance of an original point
            moves = []
            skinsq = skin**2
            dx = 0.
            while dx * gridReso <= skin:
                dy = 0.
                while (dx**2 + dy**2) <= skinsq:
                    dz = 0.
                    while (dx**2 + dy**2 + dz**2)  <= skinsq:
                        moves.append((dx, dy, dz))
                        moves.append((dx, dy, -dz))
                        moves.append((dx, -dy, dz))
                        moves.append((dx, -dy, -dz))
                        moves.append((-dx, dy, dz))
                        moves.append((-dx, dy, -dz))
                        moves.append((-dx, -dy, dz))
                        moves.append((-dx, -dy, -dz))
                        dz += gridReso
                    dy += gridReso
                dx += gridReso
            moves = list(set(moves))  #filter out redundant moves (like +- 0). Wasteful? Yes.



            # Then apply all these allowable moves to each original point
            for x,y,z in povmeList:
                for dx, dy, dz in moves:
                    #if not (x+dx, y+dy, z+dz) in povmeSet:
                    povmeSet.add((x+dx, y+dy, z+dz))

        #print 'AFTER', len(povmeSet)


        #Turning into a set so the hash search will be fast
        #povmeSet = set([ tuple(i) for i in povmeList])

        for feature in list(colorMaps.keys()):
            completePointList = colorMaps[feature].toPovmeList()
            #Only return points that were available in the POVME list beforehand
            #THIS COULD BE IMPROVED
            colorPointLists[feature] = numpy.array([row for row in completePointList.tolist() if tuple(row[:3]) in povmeSet])
            #colorPointLists[feature] = numpy.array([row for row in completePointList.tolist() if row[:3] in povmeList.tolist()])

            #thisMap = colorMaps[feature].getData()
            #nonZeroPts = numpy.transpose(thisMap.nonzero())
            #colorPointLists[feature] = numpy.array((len(nonZeroPts), 4))
            #for index, nonZeroPt in enumerate(nonZeroPts):
            #    coords = colorMaps[feature].index_to_coord(nonZeroPt)
            #    colorPointLists[feature][index][0] = thisMap[coords[0]]
            #    colorPointLists[feature][index][1] = thisMap[coords[1]]
            #    colorPointLists[feature][index][2] = thisMap[coords[2]]
            #    colorPointLists[feature][index][3] = thisMap[nonZeroPt]

        return colorPointLists



    def characterize_ligand_features(self):
        '''This is a separate function because pymolecule isn't great at loading ligands.
        Instead, this uses binana's old PDB file loader
        '''
        self.characterize_ligand_occupancies()
        self.characterize_ligand_hydrophilics()
        self.characterize_ligand_hydrophobics()
        self.characterize_ligand_hbond_acceptors()
        self.characterize_ligand_hbond_donors()
        self.characterize_ligand_aromatics()


    def characterize_ligand_hydrophilics(self):

        # Categorize hydrophilic atoms (N+O+S+P) (ligand pdb version)
        for atom in self.receptor.AllAtoms:
            if self.receptor.AllAtoms[atom].element == 'N' or self.receptor.AllAtoms[atom].element == 'O' or self.receptor.AllAtoms[atom].element == 'S' or self.receptor.AllAtoms[atom].element == 'P':
                self.features['hydrophilics'][atom] = {'coordinates':self.receptor.AllAtoms[atom].coordinates}


    def characterize_ligand_hydrophobics(self):

        # Categorize hydrophobic atoms (carbon) (ligand pdb version)
        for atom in self.receptor.AllAtoms:
            if self.receptor.AllAtoms[atom].element == 'C':
                self.features['hydrophobics'][atom] = {'coordinates':self.receptor.AllAtoms[atom].coordinates}


    def characterize_ligand_hbond_acceptors(self):
        # Categorize ligand hbond acceptors(ligand pdb version)
        for atom in self.receptor.AllAtoms:
            if self.receptor.AllAtoms[atom].element == 'O':
                self.features['hb_acceptors'][atom] = {'coordinates':self.receptor.AllAtoms[atom].coordinates}


    def characterize_ligand_hbond_donors(self):


        ####Is it a h-bond donor? (ligand version)
        for receptor_atom_index in self.receptor.AllAtoms:
            #print self.receptor.connected_atoms_of_given_element(receptor_atom_index, "H")
            #adj_hydrogen_indices = self.receptor.connected_atoms_of_given_element(receptor_atom_index, "H")
            this_atom = self.receptor.AllAtoms[receptor_atom_index]
            adjacents = [self.receptor.AllAtoms[i] for i in this_atom.IndeciesOfAtomsConnecting]
            adj_hydrogens = [i for i in adjacents if i.element == 'H']
            if ((this_atom.element == "O") or (this_atom.element == "N")) and len(adj_hydrogens) > 0:
                for adj_hydrogen in adj_hydrogens:
                    #adj_hydrogen = self.receptor.AllAtoms[adj_hydrogen_index]
                    hbd_vector = self.functions.vector_subtraction(adj_hydrogen.coordinates, this_atom.coordinates)
                    self.features['hb_donors'][receptor_atom_index] = {'coordinates':this_atom.coordinates,
                                                           'hbd_vector':hbd_vector}

    def characterize_ligand_charges(self):
        #categorize charged atoms (ligand version)
        for i, charge in enumerate(self.receptor.charges):
            self.chargeds[i] = {'coordinates':charge.coordinates, 'positive':charge.positive}

    def characterize_ligand_aromatics(self):
        #categorize aromatic rings (original version)
        for i, ar_ring in enumerate(self.receptor.aromatic_rings):
            norm_vector = point([ar_ring.plane_coeff[0], ar_ring.plane_coeff[1], ar_ring.plane_coeff[2]])
            self.features['ar_rings'][i] = {'coordinates':ar_ring.center, 'norm_vector':norm_vector, 'radius':ar_ring.radius}


    def characterize_ligand_occupancies(self):
        #Numbers taken from POVME2 / wikipedia/Bondi's compilation (1964)
        radii = {'H':1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'P':1.8, 'S':1.8}
        #Same as above, but setting H to 0.0
        #radii = {'H':0.0, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'P':1.8, 'S':1.8}
        #Half of above
        #radii = {'H':0.0, 'C':0.85, 'N':0.775, 'O':0.76, 'F':0.735, 'P':0.9, 'S':0.9}
        for this_atom in self.receptor.AllAtoms:
            this_element = self.receptor.AllAtoms[this_atom].element
            coordinates_pt = self.receptor.AllAtoms[this_atom].coordinates
            if this_element in list(radii.keys()):
                #print self.receptor.AllAtoms[this_atom].coordinates, radii[this_element]
                self.features['occupancies'][this_atom] = {'coordinates':coordinates_pt, 'radius':radii[this_element]}
            else:
                print("WARNING: Atom type %s is not C, H, O, N P, S, or F. Assuming atomic radius of 2 angstroms" %(this_element))
                self.features['occupancies'][this_atom] = {'coordinates':coordinates_pt, 'radius': 2.0}




















    def characterize_features(self):
        self.characterize_hydrophobics()
        self.characterize_hydrophilics()
        self.characterize_hbond_acceptors()
        self.characterize_hbond_donors()
        self.characterize_aromatics()
        self.characterize_occupancies()


    def characterize_hydrophilics(self):

        # categorize hydrophobic atoms (pymolecule version)
        for receptor_atom_index in self.receptor.selections.select_atoms({'element_stripped':['O','N']}):

            coordinates =  self.receptor.information.get_coordinates()[receptor_atom_index]
            coordinates_pt = point(coordinates)
            self.features['hydrophilics'][receptor_atom_index] = {'coordinates':coordinates_pt}


    def characterize_hydrophobics(self):

        # categorize hydrophobic atoms (pymolecule version)

        for receptor_atom_index in self.receptor.selections.select_atoms({'element_stripped':'C'}):

            coordinates =  self.receptor.information.get_coordinates()[receptor_atom_index]
            coordinates_pt = point(coordinates)
            self.features['hydrophobics'][receptor_atom_index] = {'coordinates':coordinates_pt}



            #Is is hydrophobic? (original version)
            #receptor_atom = receptor.AllAtoms[receptor_atom_index]

            #if receptor_atom.element == "C":
                #self.features['hydrophobics'][receptor_atom_index] = receptor_atom.coordinates

    def characterize_hbond_acceptors(self):
        #Categorize H-bond acceptors (pymolecule version)
        for receptor_atom_index in self.receptor.selections.select_atoms({'element_stripped':'O'}):
            coordinates = self.receptor.information.get_coordinates()[receptor_atom_index]
            coordinates_pt = point(coordinates)
            self.features['hb_acceptors'][receptor_atom_index] = {'coordinates':coordinates_pt}

        #Is it a h-bond acceptor? (original version)
            #if receptor_atom.element == "O":
                #self.features['hb_acceptors'][receptor_atom_index] = receptor_atom.coordinates

    def characterize_hbond_donors(self):
        ####Is it a h-bond donor (pymolecule version)
        #print self.receptor.selections.select_atoms({'element':['O','N']})
        #for receptor_atom_index in self.receptor.selections.select_atoms({'element':['O','N']}):
        #self.receptor.create_bonds_by_distance()\
        print(self.receptor.get_bonds())

        for receptor_atom_index in self.receptor.selections.select_atoms({'element_stripped':['O','N']}):
            #print "BBB", receptor_atom_index, self.receptor.information.get_atom_information()[receptor_atom_index]

            neighbors = self.receptor.selections.select_all_atoms_bound_to_selection([receptor_atom_index])
            neighborSerials = self.receptor.information.get_atom_information()[neighbors]['serial']
            hydrogen_neighbors = self.receptor.selections.select_atoms({'element_stripped':'H', 'serial':list(neighborSerials)})
            for neighbor_index in hydrogen_neighbors:
                #print "CCCC", neighbor
                coordinates = self.receptor.information.get_coordinates()[receptor_atom_index]
                coordinates_pt = point (coordinates)
                bond_vector = self.receptor.information.get_coordinates()[neighbor_index] - self.receptor.information.get_coordinates()[receptor_atom_index]
                bond_vector_pt = point(bond_vector)
                self.features['hb_donors'][neighbor_index] = {'coordinates':coordinates_pt,
                                                              'hbd_vector':bond_vector_pt}
                #self.features['hb_donors'][receptor_atom_index] = self.receptor.information.get_coordinates()[receptor_atom_index]



    def characterize_charges(self):
        #categorize charged atoms (original version)
        #for i, charge in enumerate(self.receptor.charges):
        #    self.chargeds[i] = {'coordinates':charge.coordinates, 'positive':charge.positive}
        pass


    def characterize_aromatics(self):
        #categorize aromatic rings (original version)
        #for i, ar_ring in enumerate(self.receptor.aromatic_rings):
        #    norm_vector = point(ar_ring.plane_coeff[0], ar_ring.plane_coeff[1], ar_ring.plane_coeff[2])
        #    self.features['ar_rings'][i] = {'coordinates':ar_ring.center, 'normal':norm_vector, 'radius':ar_ring.radius}

        #Categorize rings according to residue identity
        #Check phenylalanines
        phe_atoms = self.receptor.selections.select_atoms({"resname_stripped":"PHE"})
        phe_residues = set([])
        for phe_atom in phe_atoms:
            phe_residues.add(self.receptor.information.get_atom_information()[phe_atom][5])
            #print 'atom_information', self.receptor.information.get_atom_information()[phe_atom]
        #print "phe_residues",phe_residues
        for phe_res_num in phe_residues:
            this_phe_CG = self.receptor.selections.select_atoms({"resseq":phe_res_num, "name_stripped":'CG'})
            this_phe_CD1 = self.receptor.selections.select_atoms({"resseq":phe_res_num, "name_stripped":'CD1'})
            this_phe_CD2 = self.receptor.selections.select_atoms({"resseq":phe_res_num, "name_stripped":'CD2'})
            this_phe_CE1 = self.receptor.selections.select_atoms({"resseq":phe_res_num, "name_stripped":'CE1'})
            this_phe_CE2 = self.receptor.selections.select_atoms({"resseq":phe_res_num, "name_stripped":'CE2'})
            this_phe_CZ = self.receptor.selections.select_atoms({"resseq":phe_res_num, "name_stripped":'CZ'})


            coords_CG = self.receptor.information.get_coordinates()[this_phe_CG][0]
            coords_CD1 = self.receptor.information.get_coordinates()[this_phe_CD1][0]
            coords_CD2 = self.receptor.information.get_coordinates()[this_phe_CD2][0]
            coords_CE1 = self.receptor.information.get_coordinates()[this_phe_CE1][0]
            coords_CE2 = self.receptor.information.get_coordinates()[this_phe_CE2][0]
            coords_CZ = self.receptor.information.get_coordinates()[this_phe_CZ][0]
            coords = [coords_CG, coords_CD1, coords_CD2, coords_CE1, coords_CE2, coords_CZ]
            #this_phe_ring_atoms = self.receptor.selections.select_atoms({"resseq":phe_res_num})
            #print 'CD1 coords:', self.receptor.information.get_coordinates()[this_phe_CD1]
            if self.receptor.geometry.is_planar(coords_CG, coords_CD1, coords_CE2, coords_CZ):
                #print "Planar!"
                center_x = sum([float(i[0])/len(coords) for i in coords])
                center_y = sum([float(i[1])/len(coords) for i in coords])
                center_z = sum([float(i[2])/len(coords) for i in coords])
                center = [center_x, center_y, center_z]
                center_pt = point(center)
                norm_vector = numpy.cross(coords_CG-center, coords_CE1-center)
                norm_vector_pt = point(norm_vector)
                radius = 3 # for now
                self.features['ar_rings'][phe_res_num] = {'coordinates':center_pt,
                                         'norm_vector':norm_vector_pt, "radius":radius}
            #print self.receptor.information.get_atom_information()[receptor_atom_index]

        #Check tyrosines

        tyr_atoms = self.receptor.selections.select_atoms({"resname_stripped":"TYR"})
        tyr_residues = set([])
        for tyr_atom in tyr_atoms:
            tyr_residues.add(self.receptor.information.get_atom_information()[tyr_atom][5])

        #print "tyr_residues",tyr_residues
        for tyr_res_num in tyr_residues:
            this_tyr_CG = self.receptor.selections.select_atoms({"resseq":tyr_res_num, "name_stripped":'CG'})
            this_tyr_CD1 = self.receptor.selections.select_atoms({"resseq":tyr_res_num, "name_stripped":'CD1'})
            this_tyr_CD2 = self.receptor.selections.select_atoms({"resseq":tyr_res_num, "name_stripped":'CD2'})
            this_tyr_CE1 = self.receptor.selections.select_atoms({"resseq":tyr_res_num, "name_stripped":'CE1'})
            this_tyr_CE2 = self.receptor.selections.select_atoms({"resseq":tyr_res_num, "name_stripped":'CE2'})
            this_tyr_CZ = self.receptor.selections.select_atoms({"resseq":tyr_res_num, "name_stripped":'CZ'})
            coords_CG = self.receptor.information.get_coordinates()[this_tyr_CG][0]
            coords_CD1 = self.receptor.information.get_coordinates()[this_tyr_CD1][0]
            coords_CD2 = self.receptor.information.get_coordinates()[this_tyr_CD2][0]
            coords_CE1 = self.receptor.information.get_coordinates()[this_tyr_CE1][0]
            coords_CE2 = self.receptor.information.get_coordinates()[this_tyr_CE2][0]
            coords_CZ = self.receptor.information.get_coordinates()[this_tyr_CZ][0]
            coords = [coords_CG, coords_CD1, coords_CD2, coords_CE1, coords_CE2, coords_CZ]
            #this_tyr_ring_atoms = self.receptor.selections.select_atoms({"resseq":tyr_res_num})
            #print 'tyr_res_num', tyr_res_num, 'this_tyr_res_atoms', this_tyr_CG, this_tyr_CD1, this_tyr_CD2,\
            #                                                        this_tyr_CE1, this_tyr_CE2, this_tyr_CZ
            #print 'CD1 coords:', self.receptor.information.get_coordinates()[this_tyr_CD1]
            if self.receptor.geometry.is_planar(coords_CG, coords_CD1, coords_CE2, coords_CZ):
                #print "Planar!"
                center_x = sum([float(i[0])/len(coords) for i in coords])
                center_y = sum([float(i[1])/len(coords) for i in coords])
                center_z = sum([float(i[2])/len(coords) for i in coords])
                center = [center_x, center_y, center_z]
                center_pt = point(center)
                norm_vector = numpy.cross(coords_CG-center, coords_CE1-center)
                norm_vector_pt = point(norm_vector)
                radius = 3 # for now
                self.features['ar_rings'][tyr_res_num] = {'coordinates':center_pt,
                                         'norm_vector':norm_vector_pt, "radius":radius}
            #print self.receptor.information.get_atom_information()[receptor_atom_index]


        #Check histidenes
        his_atoms = self.receptor.selections.select_atoms({"resname_stripped":["HIS","HIE"]})
        his_residues = set([])
        for his_atom in his_atoms:
            his_residues.add(self.receptor.information.get_atom_information()[his_atom][5])

        #print "his_residues",his_residues
        for his_res_num in his_residues:
            this_his_CG = self.receptor.selections.select_atoms({"resseq":his_res_num, "name_stripped":'CG'})
            this_his_ND1 = self.receptor.selections.select_atoms({"resseq":his_res_num, "name_stripped":'ND1'})
            this_his_CD2 = self.receptor.selections.select_atoms({"resseq":his_res_num, "name_stripped":'CD2'})
            this_his_CE1 = self.receptor.selections.select_atoms({"resseq":his_res_num, "name_stripped":'CE1'})
            this_his_NE2 = self.receptor.selections.select_atoms({"resseq":his_res_num, "name_stripped":'NE2'})
            coords_CG = self.receptor.information.get_coordinates()[this_his_CG][0]
            coords_ND1 = self.receptor.information.get_coordinates()[this_his_ND1][0]
            coords_CD2 = self.receptor.information.get_coordinates()[this_his_CD2][0]
            coords_CE1 = self.receptor.information.get_coordinates()[this_his_CE1][0]
            coords_NE2 = self.receptor.information.get_coordinates()[this_his_NE2][0]
            coords = [coords_CG, coords_ND1, coords_CD2, coords_CE1, coords_NE2]
            #this_his_ring_atoms = self.receptor.selections.select_atoms({"resseq":his_res_num})
            #print 'his_res_num', his_res_num, 'this_his_res_atoms', this_his_CG, this_his_CD1, this_his_CD2,\
            #                                                        this_his_CE1, this_his_CE2, this_his_CZ
            #print 'CD1 coords:', self.receptor.information.get_coordinates()[this_his_CD1]
            if self.receptor.geometry.is_planar(coords_CG, coords_ND1, coords_CE1, coords_NE2):
                #print "Planar!"
                center_x = sum([float(i[0])/len(coords) for i in coords])
                center_y = sum([float(i[1])/len(coords) for i in coords])
                center_z = sum([float(i[2])/len(coords) for i in coords])
                center = [center_x, center_y, center_z]
                center_pt = point(center)
                norm_vector = numpy.cross(coords_CG-center, coords_CE1-center)
                norm_vector_pt = point(norm_vector)
                radius = 3 # for now
                self.features['ar_rings'][his_res_num] = {'coordinates':center_pt,
                                         'norm_vector':norm_vector_pt, "radius":radius}
            #print self.receptor.information.get_atom_information()[receptor_atom_index]



        #Check tryptophans

        trp_atoms = self.receptor.selections.select_atoms({"resname_stripped":["TRP","TYP"]})
        trp_residues = set([])
        for trp_atom in trp_atoms:
            trp_residues.add(self.receptor.information.get_atom_information()[trp_atom][5])

        #print "trp_residues",trp_residues
        for trp_res_num in trp_residues:
            this_trp_CG = self.receptor.selections.select_atoms({"resseq":trp_res_num, "name_stripped":'CG'})
            this_trp_CD1 = self.receptor.selections.select_atoms({"resseq":trp_res_num, "name_stripped":'CD1'})
            this_trp_CD2 = self.receptor.selections.select_atoms({"resseq":trp_res_num, "name_stripped":'CD2'})
            this_trp_NE1 = self.receptor.selections.select_atoms({"resseq":trp_res_num, "name_stripped":'NE1'})
            this_trp_CE2 = self.receptor.selections.select_atoms({"resseq":trp_res_num, "name_stripped":'CE2'})
            this_trp_CE3 = self.receptor.selections.select_atoms({"resseq":trp_res_num, "name_stripped":'CE3'})
            this_trp_CZ2 = self.receptor.selections.select_atoms({"resseq":trp_res_num, "name_stripped":'CZ2'})
            this_trp_CZ3 = self.receptor.selections.select_atoms({"resseq":trp_res_num, "name_stripped":'CZ3'})
            this_trp_CH2 = self.receptor.selections.select_atoms({"resseq":trp_res_num, "name_stripped":'CH2'})

            coords_CG = self.receptor.information.get_coordinates()[this_trp_CG][0]
            coords_CD1 = self.receptor.information.get_coordinates()[this_trp_CD1][0]
            coords_CD2 = self.receptor.information.get_coordinates()[this_trp_CD2][0]
            coords_NE1 = self.receptor.information.get_coordinates()[this_trp_NE1][0]
            coords_CE2 = self.receptor.information.get_coordinates()[this_trp_CE2][0]
            coords_CE3 = self.receptor.information.get_coordinates()[this_trp_CE3][0]
            coords_CZ2 = self.receptor.information.get_coordinates()[this_trp_CZ2][0]
            coords_CZ3 = self.receptor.information.get_coordinates()[this_trp_CZ3][0]
            coords_CH2 = self.receptor.information.get_coordinates()[this_trp_CH2][0]
            ring1_coords = [coords_CG, coords_CD1, coords_CD2, coords_CE2, coords_NE1]
            ring2_coords = [coords_CD2, coords_CE2, coords_CE3, coords_CZ2, coords_CZ3, coords_CH2]
            #this_trp_ring_atoms = self.receptor.selections.select_atoms({"resseq":trp_res_num})
            #print 'trp_res_num', trp_res_num, 'this_trp_res_atoms', this_trp_CG, this_trp_CD1, this_trp_CD2,\
            #                                                        this_trp_CE1, this_trp_CE2, this_trp_CZ
            #for coords in [ring1_coords, ring2_coords]:
            if self.receptor.geometry.is_planar(ring1_coords[0], ring1_coords[1], ring1_coords[3], ring1_coords[4]):
                #print "Planar!"
                center_x = sum([float(i[0])/len(ring1_coords) for i in ring1_coords])
                center_y = sum([float(i[1])/len(ring1_coords) for i in ring1_coords])
                center_z = sum([float(i[2])/len(ring1_coords) for i in ring1_coords])
                center = [center_x, center_y, center_z]
                center_pt = point(center)

                norm_vector = numpy.cross(ring1_coords[0]-center, ring1_coords[2]-center)
                norm_vector_pt = point(norm_vector)
                radius = 3 # for now
                #print this_trp_CG, ring1_coords[0], ring1_coords[1], norm_vector
                self.features['ar_rings'][this_trp_CG[0]] = {'coordinates':center_pt,
                                         'norm_vector':norm_vector_pt, "radius":radius}
                    #print self.receptor.information.get_atom_information()[receptor_atom_index]
            if self.receptor.geometry.is_planar(ring2_coords[0], ring2_coords[1], ring2_coords[3], ring2_coords[4]):
                #print "Planar!"
                center_x = sum([float(i[0])/len(ring2_coords) for i in ring2_coords])
                center_y = sum([float(i[1])/len(ring2_coords) for i in ring2_coords])
                center_z = sum([float(i[2])/len(ring2_coords) for i in ring2_coords])
                center = [center_x, center_y, center_z]
                center_pt = point(center)
                norm_vector = numpy.cross(ring2_coords[0]-center, ring2_coords[2]-center)
                norm_vector_pt = point(norm_vector)
                radius = 3 # for now
                self.features['ar_rings'][this_trp_CD2[0]] = {'coordinates':center_pt,
                                         'norm_vector':norm_vector_pt, "radius":radius}
                    #print self.receptor.information.get_atom_information()[receptor_atom_index]

    def characterize_occupancies(self):
        #Numbers taken from POVME2
        #From wikipedia, which is from "Bondi's compilation" (1964)
        radii = {'H':1.2, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'P':1.8, 'S':1.8}
        #Same as above, but setting H to 0.0
        #radii = {'H':0.0, 'C':1.7, 'N':1.55, 'O':1.52, 'F':1.47, 'P':1.8, 'S':1.8}
        for this_atom_info in self.receptor.get_atom_information():
            #print this_atom_info
            this_serial = this_atom_info['serial']
            this_element = this_atom_info['element_stripped']
            receptor_atom_index = self.receptor.selections.select_atoms({'serial':this_serial})
            coordinates_pt = point(self.receptor.get_coordinates()[receptor_atom_index][0])
            if this_element in list(radii.keys()):
                self.features['occupancies'][this_serial] = {'coordinates': coordinates_pt, 'radius':radii[this_element]}
            else:
                print("WARNING: Atom type %s is not C, H, O, N P, S, or F. Assuming atomic radius of 2 angstroms" %(this_element))
                self.features['occupancies'][this_atom_info['serial']] = {'coordinates':coordinates_pt, 'radius':2.0}

    def create_atom_dictionary(self):
        #atomDictionary = {}
        
        #List neighboring atoms (pymolecule version)
        for main_atom_index in self.receptor.selections.select_all():
            
            main_type = self.receptor.information.get_atom_information()[main_atom_index][13]
            type_list = main_type
            
            #Get element type of neighboring atoms
            for atom in self.receptor.selections.select_all_atoms_bound_to_selection([main_atom_index]):
                atom_type = self.receptor.information.get_atom_information()[atom][13]
                type_list = type_list + atom_type

            key = type_list
            key = key[0] + ''.join(sorted(key[1:]))
            
            #Check if key already exists
            if key in self.atom_dictionary:
                self.atom_dictionary[key].append(main_atom_index)
            else:
                self.atom_dictionary[key] = [main_atom_index] 
            
        #return atomDictionary




    def draw_hydrophilics(self, parameters):
        self.vmdScript += ['mol new']
        self.vmdScript += ['mol rename top hydrophilic']
        self.vmdScript += ['draw material Glass1']
        self.vmdScript += ['draw color blue']
        for hydrophilic in self.features['hydrophilics']:
            coords = self.features['hydrophilics'][hydrophilic]['coordinates']
            self.vmdScript += ['draw sphere {%f %f %f} radius %f ' % (coords.x, coords.y, coords.z, parameters['hydroph_gaussian_cutoff'])]


    def draw_hydrophobics(self, parameters):
        self.vmdScript += ['mol new']
        self.vmdScript += ['mol rename top hydrophobic']
        self.vmdScript += ['draw material Glass1']
        self.vmdScript += ['draw color orange']
        for hydrophobic in self.features['hydrophobics']:
            coords = self.features['hydrophobics'][hydrophobic]['coordinates']
            self.vmdScript += ['draw sphere {%f %f %f} radius %f ' % (coords.x, coords.y, coords.z, parameters['hydroph_gaussian_cutoff'])]



    def draw_hb_acceptors(self):
        self.vmdScript += ['mol new']
        self.vmdScript += ['mol rename top hb_acceptors']
        self.vmdScript += ['draw material Glass2']
        self.vmdScript += ['draw color purple']
        for hb_acceptor in self.features['hb_acceptors']:
            this_acceptor = self.features['hb_acceptors'][hb_acceptor]['coordinates']
            self.vmdScript += ['draw sphere {%f %f %f} radius %f' % (this_acceptor.x, this_acceptor.y, this_acceptor.z, 1)]


    def draw_hb_donors(self, parameters):
        self.vmdScript += ['mol new']
        self.vmdScript += ['mol rename top hb_donors']
        self.vmdScript += ['draw material Glass1']
        self.vmdScript += ['draw color cyan']
        for hb_donor in self.features['hb_donors']:
            this_donor = self.features['hb_donors'][hb_donor]
            coords = this_donor['coordinates']
            unit_vector = self.functions.return_normalized_vector(this_donor['hbd_vector'])
            vector = self.functions.vector_scalar_multiply(unit_vector, parameters['hydrogen_bond_dist_cutoff'])
            base_radius = math.fabs(math.tan(math.radians(parameters['hydrogen_bond_angle_cutoff'])) * parameters['hydrogen_bond_dist_cutoff'])
            self.vmdScript += ['draw cone {%f %f %f} {%f %f %f} radius %f' % (coords.x + vector.x, coords.y + vector.y,
                    coords.z + vector.z,
                    coords.x, coords.y, coords.z,
                    base_radius)]


    def draw_aromatic_rings(self):
        self.vmdScript += ['mol new']
        self.vmdScript += ['mol rename top ar_rings']
        self.vmdScript += ['draw color tan']
        use_triangles = False
        if use_triangles:
            toroid_triangles = toroidPoints.triangles
            for ar_ring in self.features['ar_rings']:
                center = self.features['ar_rings'][ar_ring]['coordinates']
                for triangle in toroid_triangles:
                    self.vmdScript += ['draw triangle {%f %f %f} {%f %f %f} {%f %f %f}' % (triangle[0][0] + center[0], triangle[0][1] + center[1], triangle[0][2] + center[2], triangle[1][0] + center[0], triangle[1][1] + center[1], triangle[1][2] + center[2], triangle[2][0] + center[0], triangle[2][1] + center[1], triangle[2][2] + center[2])]

        use_spheres = True
        if use_spheres:
            for ar_ring in self.features['ar_rings']:
                center = self.features['ar_rings'][ar_ring]['coordinates']
                norm_vector = self.features['ar_rings'][ar_ring]['norm_vector'].coords() / numpy.linalg.norm(self.features['ar_rings'][ar_ring]['norm_vector'].coords())
                #print 'norm_vector:', norm_vector
                #The ring is initially oriented with the Z axis running through the opening
                rot_axis = numpy.cross([0, 0, 1], norm_vector)
                #print 'rot_axis', rot_axis
                angle = numpy.arccos(numpy.dot(norm_vector, [0, 0, 1]) / numpy.linalg.norm(norm_vector)) # Using formula cos(theta) = a * b / (|a| |b|)
                r = self.functions.axisangle_to_q(rot_axis, angle)
                n = 20
                #print angle, r
                for i in numpy.arange(0, 2 * numpy.pi, numpy.pi / n):
                    sphere_x = numpy.sin(i)
                    sphere_y = numpy.cos(i)
                    sphere_z = 0
                    sphere_vector = sphere_x, sphere_y, sphere_z
                    sphere_vector = self.functions.qv_mult(r, sphere_vector) #print r, sphere_vector
                    self.vmdScript += ['draw sphere {%f %f %f} radius 0.5' % (center[0] + sphere_vector[0], center[1] + sphere_vector[1], center[2] + sphere_vector[2])]

    def load_pymolecule(self, filename, bonds_by_distance=True, serial_reindex=True, resseq_reindex=False):
        #Not really sure how to implement this.
        loadedAtom = pymolecule.Molecule()
        loadedAtom.fileio.load_pdb_into(filename, bonds_by_distance, serial_reindex, resseq_reindex)
        
        return loadedAtom

    def write_vmd_script_file(self, vmdScriptName):
        text = '\n'.join(self.vmdScript)
        fo = open(vmdScriptName, 'w')
        fo.write(text)
        fo.close()

    def write_vmd_script(self, vmdScriptName, parameters):
        self.vmdScript = []
        #Make hydrophobic representation
        self.draw_hydrophobics(parameters)

        #Make hydrophilic representation
        self.draw_hydrophilics(parameters)

        #hb acceptors
        self.draw_hb_acceptors()

        #hb donors
        self.draw_hb_donors(parameters)

        #aromantic rings
        self.draw_aromatic_rings()

        self.write_vmd_script_file(vmdScriptName)

    def rotate(self, quaternion, center='average'):
        if center == 'average':
            coords = []
            for feature in list(self.features.keys()):
                    for key in list(self.features[feature].keys()):
                        coords.append(self.features[feature][key]['coordinates'].coords())
            coords = numpy.array(coords)
            center = numpy.array([numpy.mean(coords[:,0]), numpy.mean(coords[:,1]), numpy.mean(coords[:,2])])

        #Move whole thing to be vcentered around 0,0,0
        self.translate(-1.0*center)
        for feature in list(self.features.keys()):
                for index in list(self.features[feature].keys()):
                    oldCoordsList = self.features[feature][index]['coordinates'].coords()
                    oldMagnitude = self.features[feature][index]['coordinates'].magnitude()
                    #rotation = self.functions.axisangle_to_q(quaternion[1:], numpy.arcsin(quaternion[0])*2)
                    rotation = self.functions.axisangle_to_q(quaternion[1:], numpy.arccos(quaternion[0])*2)
                    newCoordsList = numpy.array(self.functions.qv_mult(rotation, oldCoordsList)) * oldMagnitude
                    #coordsList = [coords.x, coords.y, coords.z]
                    self.features[feature][index]['coordinates'] = point(newCoordsList)
                    for key in self.features[feature][index]:
                        if 'vector' in key:
                            coordsList = self.features[feature][index][key].coords()
                            #coordsList = [coords.x, coords.y, coords.z]
                            self.features[feature][index][key] = point(self.functions.qv_mult(rotation, coordsList))
        #Move whole shebang back to where it was originally
        self.translate(center)


    def translate(self, vector):
        for feature in list(self.features.keys()):
            for key in list(self.features[feature].keys()):
                oldCoords = self.features[feature][key]['coordinates'].coords()
                newCoords = oldCoords + numpy.array(vector)
                self.features[feature][key]['coordinates'] = point(newCoords)

    def suggest_borders(self, padding = 6.0):
        '''Finds the maximum extent of all the features in x, y, and z, and then adds a padding constant onto it.'''
        coords = []
        borders = numpy.array([99999.,-99999.,99999.,-99999.,99999.,-99999.])

        for feature in self.features:
            for index in self.features[feature]:

                coordsList = self.features[feature][index]['coordinates'].coords()
                if coordsList[0] < borders[0]:
                    borders[0] = coordsList[0]
                if coordsList[0] > borders[1]:
                    borders[1] = coordsList[0]
                if coordsList[1] < borders[2]:
                    borders[2] = coordsList[1]
                if coordsList[1] > borders[3]:
                    borders[3] = coordsList[1]
                if coordsList[2] < borders[4]:
                    borders[4] = coordsList[2]
                if coordsList[2] > borders[5]:
                    borders[5] = coordsList[2]

        return borders + numpy.array([-padding, padding, -padding, padding, -padding, padding])
    
    def suggest_borders_atoms(self, padding = 6.0):
        '''Finds the maximum extent of all the features in x, y, and z, and then adds a padding constant onto it.'''
        coords = []
        borders = numpy.array([99999.,-99999.,99999.,-99999.,99999.,-99999.])

        for key in self.atom_dictionary:
            for index in self.atom_dictionary[key]:
                coordsList = self.receptor.information.get_coordinates()[index]
                coordsList = point(coordsList)
                if coordsList[0] < borders[0]:
                    borders[0] = coordsList[0]
                if coordsList[0] > borders[1]:
                    borders[1] = coordsList[0]
                if coordsList[1] < borders[2]:
                    borders[2] = coordsList[1]
                if coordsList[1] > borders[3]:
                    borders[3] = coordsList[1]
                if coordsList[2] < borders[4]:
                    borders[4] = coordsList[2]
                if coordsList[2] > borders[5]:
                    borders[5] = coordsList[2]

        return borders + numpy.array([-padding, padding, -padding, padding, -padding, padding])

    def __init__(self, receptor, parameters, isLigand=False):

        self.atomDictionaryBinned = False
        self.hbondAcceptorsBinned = False
        self.hbondDonorsBinned = False
        self.aromaticsBinned = False
        self.hydrophobicsBinned = False
        self.hydrophilicsBinned = False
        self.occupanciesBinned = False
        self.receptor=receptor
        #self.receptor = pymolecule.Molecule()

        #self.receptor.fileio.load_pdb_into(receptor_pdbqt_filename, serial_reindex=False, resseq_reindex=False)
        #self.receptor.fileio.load_pdb_into(receptor_pdbqt_filename, serial_reindex=True, resseq_reindex=False)

        #receptor.assign_secondary_structure()
        self.parameters = parameters
        self.functions = MathFunctions()

        # Make dictionaries of hbonds, hydrophobic, pi, and charge regions
        self.features = {}
        #HB_DONOR FORMAT: location, phi, psi
        self.features['hb_donors'] = {}
        #HB ACCEPTOR FORMAT: location
        self.features['hb_acceptors'] = {}
        #HYDROPHOBIC FORMAT: location
        self.features['hydrophobics'] = {}
        #HYDROPHOBIC FORMAT: location
        self.features['hydrophilics'] = {}
        #PI FORMAT: location (ring center), radius, normal
        self.features['ar_rings'] = {}
        #CHARGED FORMAT: location, positive
        self.features['chargeds'] = {}
        #OCCUPANCY FORMAT: location, radius
        self.features['occupancies'] = {}
        if isLigand:
            self.characterize_ligand_features()
        else:
            self.characterize_features()
        
        self.atom_dictionary = {}
        self.create_atom_dictionary()
        
        print("Peeled!")


class command_line_parameters:

    params = {}

    def is_num(self, num):
        try:
            t = float(num)
            return t
        except ValueError:
            return num

    def __init__(self, parameters):

        # first, set defaults
        self.params['close_contacts_dist1_cutoff'] = 2.5
        self.params['close_contacts_dist2_cutoff'] = 4.0
        self.params['electrostatic_dist_cutoff'] = 4.0
        self.params['active_site_flexibility_dist_cutoff'] = 4.0
        self.params['hydrophobic_dist_cutoff'] = 4.0
        self.params['hydrogen_bond_dist_cutoff'] = 4.0
        self.params['hydrogen_bond_angle_cutoff'] = 40.0
        self.params['pi_padding_dist'] = 0.75
        self.params['pi_pi_interacting_dist_cutoff'] = 7.5
        self.params['pi_stacking_angle_tolerance'] = 30.0
        self.params['T_stacking_angle_tolerance'] = 30.0
        self.params['T_stacking_closest_dist_cutoff'] = 5.0
        self.params['cation_pi_dist_cutoff'] = 6.0
        self.params['salt_bridge_dist_cutoff'] = 5.5
        self.params['receptor'] = ''
        #self.params['ligand'] = ''
        self.params['output_dir'] = ''
        self.params['output_file'] = ''

        # now get user inputed values

        for index in range(len(parameters)):
            item = parameters[index]
            if item[:1] == '-': # so it's a parameter key value
                key = item.replace('-','')
                value = self.is_num(parameters[index+1])
                if key in list(self.params.keys()):
                    self.params[key] = value
                    parameters[index] = ""
                    parameters[index + 1] = ""

        # make a list of all the command-line parameters not used
        self.error = ""
        for index in range(1,len(parameters)):
            item = parameters[index]
            if item != "": self.error = self.error + item + " "

        # Make sure the output directory, if specified, ends in a /
        if self.params['output_dir'] != "":
            if self.params['output_dir'][-1:] != os.sep:
                self.params['output_dir'] = self.params['output_dir'] + os.sep

        # If an output directory is specified but a log file isn't, set a default logfile
        if self.params['output_dir'] != "" and self.params['output_file'] == '': self.params['output_file'] = self.params['output_dir'] + 'output.pdb'

    def okay_to_proceed(self):
        # at the very least, you need the ligand and the receptor
        if self.params['receptor'] != '':
            return True
        else:
            return False

def intro():
    version = "1.1.0"
    lines = []
    lines.append("")
    lines.append("BINANA " + version)
    lines.append("============")
    lines.append("   BINANA is released under the GNU General Public License (see http://www.gnu.org/licenses/gpl.html). If you have any questions, comments, or suggestions, please don't hesitate to contact me, Jacob Durrant, at jdurrant [at] ucsd [dot] edu. If you use BINANA in your work, please cite [REFERENCE HERE].")
    lines.append("")
    lines.append("Introduction")
    lines.append("============")
    lines.append("   BINANA (BINding ANAlyzer) is a python-implemented algorithm for analyzing ligand binding. The program identifies key binding characteristics like hydrogen bonds, salt bridges, and pi interactions. As input, BINANA accepts receptor and ligand files in the PDBQT format. PDBQT files can be generated from the more common PDB file format using the free converter provided with AutoDockTools, available at http://mgltools.scripps.edu/downloads")
    lines.append("   As output, BINANA describes ligand binding. Here's a simple example of how to run the program:")
    lines.append("")
    lines.append("python binana.py -receptor /path/to/receptor.pdbqt -ligand /path/to/ligand.pdbqt")
    lines.append("")
    lines.append("   To create a single PDB file showing the different binding characteristics with those characteristics described in the PDB header:")
    lines.append("")
    lines.append("python binana.py -receptor /path/to/receptor.pdbqt -ligand /path/to/ligand.pdbqt -output_file /path/to/output.pdb")
    lines.append("")
    lines.append("   Note that in the above example, errors and warnings are not written to the output file. To save these to a file, try:")
    lines.append("")
    lines.append("python binana.py -receptor /path/to/receptor.pdbqt -ligand /path/to/ligand.pdbqt -output_file /path/to/output.pdb > errors.txt")
    lines.append("")
    lines.append("   You can also send the program output to a directory, which will be created if it does not already exist. If a directory is specified, the program automatically separates the output PDB file into separate files for each interaction analyzed, and a description of the interactions is written to a file called 'log.txt'. Additionally, a VMD state file is created so the results can be easily visualized in VMD, a free program available for download at http://www.ks.uiuc.edu/Development/Download/download.cgi?PackageName=VMD Again, to save warnings and errors, append something like \"> errors.txt\" to the end of your command:")
    lines.append("")
    lines.append("python binana.py -receptor /path/to/receptor.pdbqt -ligand /path/to/ligand.pdbqt -output_dir /path/to/output/directory/ > errors.txt")
    lines.append("")
    lines.append("   Though we recommend using program defaults, the following command-line tags can also be specified: -close_contacts_dist1_cutoff -close_contacts_dist2_cutoff -electrostatic_dist_cutoff -active_site_flexibility_dist_cutoff -hydrophobic_dist_cutoff -hydrogen_bond_dist_cutoff -hydrogen_bond_angle_cutoff -pi_padding_dist -pi_pi_interacting_dist_cutoff -pi_stacking_angle_tolerance -T_stacking_angle_tolerance -T_stacking_closest_dist_cutoff -cation_pi_dist_cutoff -salt_bridge_dist_cutoff")
    lines.append("   For example, if you want to tell BINANA to detect only hydrogen bonds where the donor and acceptor are less than 3.0 angstroms distant, run:")
    lines.append("")
    lines.append("python binana.py -receptor /path/to/receptor.pdbqt -ligand /path/to/ligand.pdbqt -hydrogen_bond_dist_cutoff 3.0")
    lines.append("")
    lines.append("   What follows is a detailed description of the BINANA algorithm and a further explaination of the optional parameters described above. Parameter names are enclosed in braces.")
    lines.append("")
    lines.append("Close Contacts")
    lines.append("==============")
    lines.append("   BINANA begins by identifying all ligand and protein atoms that come within {close_contacts_dist1_cutoff} angstroms of each other. These close-contact atoms are then characterized according to their respective AutoDock atom types, without regard for the receptor or ligand. The number of each pair of close-contact atoms of given AutoDock atom types is then tallied. For example, the program counts the number of times a hydrogen-bond accepting oxygen atom (atom type OA), either on the ligand or the receptor, comes within {close_contacts_dist1_cutoff} angstroms of a polar hydrogen atom (atom type HD) on the corresponding binding partner, be it the receptor or the ligand. A similar list of atom-type pairs is tallied for all ligand and receptor atoms that come within {close_contacts_dist2_cutoff} angstroms of each other, where {close_contacts_dist2_cutoff} > {close_contacts_dist1_cutoff}.")
    lines.append("")
    lines.append("Electrostatic Interactions")
    lines.append("==========================")
    lines.append("   For each atom-type pair of atoms that come within {electrostatic_dist_cutoff} angstroms of each other, as described above, a summed electrostatic energy is calculated using the Gasteiger partial charges assigned by AutoDockTools.")
    lines.append("")
    lines.append("Binding-Pocket Flexibility")
    lines.append("==========================")
    lines.append("   BINANA also provides useful information about the flexibility of a binding pocket. Each receptor atom that comes with {active_site_flexibility_dist_cutoff} angstroms of any ligand atom is characterized according to whether or not it belongs to a protein side chain or backbone. Additionally, the secondary structure of the corresponding protein residue of each atom, be it alpha helix, beta sheet, or other, is also determined. Thus, there are six possible characterizations for each atom: alpha-sidechain, alpha-backbone, beta-sidechain, beta-backbone, other-sidechain, other-backbone. The number of close-contact receptor atoms falling into each of these six categories is tallied as a metric of binding-site flexibility.")
    lines.append("   All protein atoms with the atom names \"CA,\" \"C,\" \"O,\" or \"N\" are assumed to belong to the backbone. All other receptor atoms are assigned side-chain status. Determining the secondary structure of the corresponding residue of each close-contact receptor atom is more difficult. First, preliminary secondary-structure assignments are made based on the phi and psi angles of each residue. If phi in (-145, -35) and psi in (-70, 50), the residue is assumed to be in the alpha-helix conformation. If phi in [-180, -40) and psi in (90,180], or phi in [-180,-70) and psi in [-180, -165], the residue is assumed to be in the beta-sheet conformation. Otherwise, the secondary structure of the residue is labeled \"other.\"")
    lines.append("   Inspection of actual alpha-helix structures revealed that the alpha carbon of an alpha-helix residue i is generally within 6.0 angstroms of the alpha carbon of an alpha-helix residue three residues away (i + 3 or i - 3). Any residue that has been preliminarily labeled \"alpha helix\" that fails to meet this criteria is instead labeled \"other.\" Additionally, the residues of any alpha helix comprised of fewer than four consecutive residues are also labeled \"other,\" as these tended belong to be small loops rather than genuine helices.")
    lines.append("   True beta strands hydrogen bond with neighboring beta strands to form beta sheets. Inspection of actual beta strands revealed that the Calpha of a beta-sheet residue, i, is typically within 6.0 angstroms of the Calpha of another beta-sheet residue, usually on a different strand, when the residues [i - 2, i + 2] are excluded. Any residue labeled \"beta sheet\" that does not meet this criteria is labeled \"other\" instead. Additionally, the residues of beta strands that are less than three residues long are likewise labeled \"other,\" as these residues typically belong to loops rather than legitimate strands.")
    lines.append("")
    lines.append("Hydrophobic Contacts")
    lines.append("====================")
    lines.append("   To identify hydrophobic contacts, BINANA simply tallies the number of times a ligand carbon atom comes within {hydrophobic_dist_cutoff} angstroms of a receptor carbon atom. These hydrophobic contacts are categorized according to the flexibility of the receptor carbon atom. There are six possible classifications: alpha-sidechain, alpha-backbone, beta-sidechain, beta-backbone, other-sidechain, other-backbone. The total number of hydrophobic contacts is simply the sum of these six counts.")
    lines.append("")
    lines.append("Hydrogen Bonds")
    lines.append("==============")
    lines.append("   BINANA allows hydroxyl and amine groups to act as hydrogen-bond donors. Oxygen and nitrogen atoms can act as hydrogen-bond acceptors. Fairly liberal cutoffs are implemented in order to accommodate low-resolution crystal structures. A hydrogen bond is identified if the hydrogen-bond donor comes within {hydrogen_bond_dist_cutoff} angstroms of the hydrogen-bond acceptor, and the angle formed between the donor, the hydrogen atom, and the acceptor is no greater than {hydrogen_bond_angle_cutoff} degrees. BINANA tallies the number of hydrogen bonds according to the secondary structure of the receptor atom, the side-chain/backbone status of the receptor atom, and the location (ligand or receptor) of the hydrogen bond donor. Thus there are twelve possible categorizations: alpha-sidechain-ligand, alpha-backbone-ligand, beta-sidechain-ligand, beta-backbone-ligand, other-sidechain-ligand, other-backbone-ligand, alpha-sidechain-receptor, alpha-backbone-receptor, beta-sidechain-receptor, beta-backbone-receptor, other-sidechain-receptor, other-backbone-receptor.")
    lines.append("")
    lines.append("Salt Bridges")
    lines.append("============")
    lines.append("   BINANA also seeks to identify possible salt bridges binding the ligand to the receptor. First, charged functional groups are identified and labeled with a representative point to denote their location. For non-protein residues, BINANA searches for common functional groups or atoms that are known to be charged. Atoms containing the following names are assumed to be metal cations: MG, MN, RH, ZN, FE, BI, AS, AG. The identifying coordinate is centered on the metal cation itself. Sp3-hybridized amines (which could pick up a hydrogen atom) and quaternary ammonium cations are also assumed to be charged; the representative coordinate is centered on the nitrogen atom. Imidamides where both of the constituent amines are primary, as in the guanidino group, are also fairly common charged groups. The representative coordinate is placed between the two constituent nitrogen atoms. ")
    lines.append("   Carboxylate groups are likewise charged; the identifying coordinate is placed between the two oxygen atoms. Any group containing a phosphorus atom bound to two oxygen atoms that are themselves bound to no other heavy atoms (i.e., a phosphate group) is also likely charged; the representative coordinate is centered on the phosphorus atom. Similarly, any group containing a sulfur atom bound to three oxygen atoms that are themselves bound to no other heavy atoms (i.e., a sulfonate group) is also likely charged; the representative coordinate is centered on the sulfur atom. Note that while BINANA is thorough in its attempt to identify charged functional groups on non-protein residues, it is not exhaustive. For example, one could imagine a protonated amine in an aromatic ring that, though charged, would not be identified as a charged group.")
    lines.append("   Identifying the charged functional groups of protein residues is much simpler. Functional groups are identified based on standardized protein atom names. Lysine residues have an amine; the representative coordinate is centered on the nitrogen atom. Arginine has a guanidino group; the coordinate is centered between the two terminal nitrogen atoms. Histadine is always considered charged, as it could pick up a hydrogen atom. The representative charge is placed between the two ring nitrogen atoms. Finally, glutamate and aspartate contain charged carboxylate groups; the representative coordinate is placed between the two oxygen atoms.")
    lines.append("   Having identified the location of all charged groups, BINANA is ready to predict potential salt bridges. First, the algorithm identifies all representative charge coordinates within {salt_bridge_dist_cutoff} angstroms of each other. Next, it verifies that the two identified coordinates correspond to charges that are opposite. If so, a salt bridge is detected. These salt bridges are characterized and tallied by the secondary structure of the associated protein residue: alpha helix, beta sheet, or other.")
    lines.append("")
    lines.append("pi Interactions")
    lines.append("===============")
    lines.append("   A number of interactions are known to involve pi systems. In order to detect the aromatic rings of non-protein residues, a recursive subroutine identifies all five or six member rings, aromatic or not. The dihedral angles between adjacent ring atoms, and between adjacent ring atoms and the first atom of ring substituents, are checked to ensure that none deviate from planarity by more than 15 degrees. Planarity establishes aromaticity. For protein residues, aromatic rings are identified using standardized protein-atom names. Phenylalanine, tyrosine, and histidine all have aromatic rings. Tryptophan is assigned two aromatic rings. ")
    lines.append("   Once an aromatic ring is identified, it must be fully characterized. First, a plane is defined that passes through three ring atoms, preferably the first, third, and fifth atoms. The center of the ring is calculated by averaging the coordinates of all ring atoms, and the radius is given to be the maximum distance between the center point and any of those atoms. From this information, a ring disk can be defined that is centered on the ring center point, oriented along the ring plane, and has a radius equal to that of the ring plus a small buffer ({pi_padding_dist} angstroms).")
    lines.append("   Having identified and characterized all aromatic rings, the algorithm next attempts to identify pi-pi stacking interactions. First, every aromatic ring of the ligand is compared to every aromatic ring of the receptor. If the centers of two rings are within {pi_pi_interacting_dist_cutoff} angstroms of each other, the angle between the two vectors normal to the planes of each ring is calculated. If this angle suggests that the two planes are within {pi_stacking_angle_tolerance} degrees of being parallel, then each of the ring atoms is projected onto the plane of the opposite ring by identifying the nearest point on that plane. If any of these projected points fall within the ring disk of the opposite ring, the two aromatic rings are said to participate in a pi-pi stacking interaction. We note that it is not sufficient to simply project the ring center point onto the plane of the opposite ring because pi-pi stacking interactions are often off center.")
    lines.append("   To detect T-stacking or edge-face interactions, every aromatic ring of the ligand is again compared to every aromatic ring of the receptor. If the centers of two rings are within {pi_pi_interacting_dist_cutoff} angstroms of each other, the angle between the two vectors normal to the planes of each ring is again calculated. If this angle shows that the two planes are within {T_stacking_angle_tolerance} degrees of being perpendicular, a second distance check is performed to verify that the two rings come within {T_stacking_closest_dist_cutoff} angstroms at their nearest point. If so, each of the coordinates of the ring center points is projected onto the plane of the opposite ring by identifying the nearest point on the respective plane. If either of the projected center points falls within the ring disk of the opposite ring, the two aromatic rings are said to participate in a T-stacking interaction.")
    lines.append("   Finally, BINANA also detects cation-pi interactions between the ligand and the receptor. Each of the representative coordinates identifying charged functional groups is compared to each of the center points of the identified aromatic rings. If the distance between any of these pairs is less than {cation_pi_dist_cutoff} angstroms, the coordinate identifying the charged functional group is projected onto the plane of the aromatic ring by identifying the nearest point on that plane. If the projected coordinate falls within the ring disk of the aromatic ring, a cation-pi interaction is identified. ")
    lines.append("   pi-pi stacking, T-stacking, and cation-pi interactions are tallied according to the secondary structure of the receptor residue containing the associated aromatic ring or charged functional group: alpha helix, beta sheet, or other.")
    lines.append("")
    lines.append("Ligand Atom Types and Rotatable Bonds")
    lines.append("=====================================")
    lines.append("   BINANA also tallies the number of atoms of each AutoDock type present in the ligand as well as the number of ligand rotatable bonds identified by the AutoDockTools scripts used to generate the input PDBQT files.")
    lines.append("")
    lines.append("                            [- END INTRO -]")
    lines.append("")
    wrapped = []
    for line in lines:
        if line == "":
            wrapped.append("")
        elif "python binana.py" in line:
            wrapped.append(line)
        else:
            wrapped.extend(textwrap.wrap(line, 80))


    print("              ..I?+?.                                       ")
    print("             ..$7$7.                                        ")
    print("             ..$$Z,                                         ")
    print("            ..I$$O?.                                        ")
    print("            :??7$$$.                                        ")
    print("          .=????7III                                        ")
    print("        ..I?????I???..                                      ")
    print("      ...???????I????.                                      ")
    print("      ..?+??????I?????                                      ")
    print("      .??I??????II????+.                                    ")
    print("      ~?II??????II+???=+..                                  ")
    print("     .??I+??????II????+=+.                                  ")
    print("     .+?I,???????I?????+~:,                                 ")
    print("     .+I?~???????+I????++~~:~..                             ")
    print("     .+??==???????=+?+++++:~~:~~:.....                      ")
    print("     .++?+.???????+~+?+++++,:~:~~~,:,:::,........ .......   ")
    print("     .?=++=~++??+++:=??++++++~::,:~~::::~~:~=~~~~~~~~~?+=.. ")
    print("      :+====?+++++++.=??++++++????=~~=~~::,:,,,,:::~==++=   ")
    print("      .+==~=+++++++++,=I?+++++++..:?77I????????????I7?...   ")
    print("      .:====~??+++++?++:I?+++++++~. ..................      ")
    print("      ..++=~=:+?++++++++??I?++++++??..                      ")
    print("      ..?~==~=~.??++++++++?7??++++++++,.                    ")
    print("        ~=====~,.=?++++++++??7?+++++++++~                   ")
    print("        .==~~:==...=?+?+++++???7??+++++++?.                 ")
    print("        .?===~=:~   ..??++?+++???7?I?++?++?,.               ")
    print("        .?=~=~~:=..     .7???++?+??7$7I?++III,.             ")
    print("        .==~~~:~:..        ..,777I???7$ODDMDNM,             ")
    print("        .,~~~~~=::.               .~$II7ZOONMN.             ")
    print("        ..:=~~~=:?.               .......,=I7I.             ")
    print("        ..=+=~==::.                      . ....             ")
    print("        ..?+===~~=.                                         ")
    print("          ~+==~~=~.                                         ")
    print("          .+?==~=:.                                         ")
    print("          ..I7II++.                                         ")
    print("            .=I$$$.                                         ")
    print("               ....                                         ")

    for i in wrapped:
        print(i)

if __name__=='__main__':
    intro()

    cmd_params = command_line_parameters(sys.argv[:])

    if cmd_params.okay_to_proceed() == False:
        print("Error: You need to specify the ligand and receptor PDBQT files to analyze using\nthe -receptor and -ligand tags from the command line.\n")
        sys.exit(0)

    if cmd_params.error != "":
        print("Warning: The following command-line parameters were not recognized:")
        print("   " + cmd_params.error + "\n")

    #lig = cmd_params.params['ligand']
    rec = cmd_params.params['receptor']

    #ligand = PDB()
    receptor = PDB()

    #ligand.LoadPDB(lig)
    receptor.LoadPDB(rec)

    parameters = cmd_params.params
    my_peel = peel(rec, parameters)
    #Prepare visualization
    #self.receptor.SavePDB('receptor.pdb') #(original version)
    my_peel.receptor.fileio.save_pdb('receptor.pdb')

    my_peel.write_vmd_script('visualize.vmd', parameters)
    print("Done!")

