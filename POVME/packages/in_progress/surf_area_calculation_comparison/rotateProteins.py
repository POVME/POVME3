import glob
import MDAnalysis as MDA
import numpy

#targets = glob.glob('../CSAR_FULL_RELEASE_29NOVEMBER2012/*/SETUP_DOCKING_FILES/PROTEIN_ALONE/*pdb')
targets = glob.glob('../CSAR_FULL_RELEASE_29NOVEMBER2012/*/SETUP_DOCKING_FILES/COMPLEX/*pdb')



for target in targets:
    print(target)
    for rotation in range(5):
        u = MDA.Universe(target)
        atom1 = u.atoms[numpy.random.randint(0,len(u.atoms))]
        atom2 = u.atoms[numpy.random.randint(0,len(u.atoms))]
        angle = numpy.random.randint(0,360)
        u.atoms.rotateby(angle, (atom1, atom2))

        rotFilename = target.replace('.pdb','_rot%i.pdb' %(rotation+1))
        writer = MDA.Writer(rotFilename)
        writer.write(u)
