import MDAnalysis as MDA
import glob
import numpy
import os

povmeTemplate = '''GridSpacing            1.0
PointsInclusionSphere  !X! !Y! !Z! !RAD!
SavePoints            true
PDBFileName                 !TARGET!
DistanceCutoff              1.09
ConvexHullExclusion         each
ContiguousPocketSeedSphere  !X! !Y! !Z! 4.0
ContiguousPointsCriteria    3
CalculateSurfaceArea      true
NumProcessors               1
UseDiskNotMemory            false
OutputFilenamePrefix          ./!ID!/!ID!_
SaveIndividualPocketVolumes   true  
SavePocketVolumesTrajectory   false
SavePocketVolumesNumpy        false
OutputEqualNumPointsPerFrame  true
SaveTabbedVolumeFile          true
SaveVolumetricDensityDX       false
SaveVolumetricDensityNpy      false
SaveColoredMap                true
CompressOutput                false
'''

targets = glob.glob('../../../../CSAR_FULL_RELEASE_29NOVEMBER2012/*/SETUP_DOCKING_FILES/PROTEIN_ALONE/*pdb')

bSiteRes = {'CDK2':[10,81], 'CHK1':[15,80], 'ERK2':[25,101],
            'LPXC':[197,212], 'UROKINASE':[91,195]}

povmeShellScript = ''
mdPocketShellScript = ''
dPocketShellScript = ''
for target in targets:
    complexTarget = target.replace('PROTEIN_ALONE','COMPLEX').replace('_pro','')
    protName = target.split('/')[-4]
    pdbName = target.split('/')[-1]
    protID = pdbName.replace('.pdb','')
    if protName in bSiteRes:
        print(target, 'IN LIST UNDER TITLE', protName)
        u = MDA.Universe(target)
        protein = u.selectAtoms('not resname LIG')
        bSiteResAtoms = u.selectAtoms(' or '.join(['(resnum %i and name CA)'%(i) for i in bSiteRes[protName]]))
        coords = bSiteResAtoms.coordinates()
        bSiteCenter = numpy.mean(coords, axis=0)
        
        povmeScript = povmeTemplate.replace(
              '!X!',str(bSiteCenter[0])).replace(
              '!Y!',str(bSiteCenter[1])).replace(
              '!Z!',str(bSiteCenter[2])).replace(
              '!TARGET!','../'+target).replace(
              '!RAD!','10').replace(
              '!ID!',protID)
        povmeScriptName ='%s_povmeScript.ini' %(protID)
        with open('povme/%s' %(povmeScriptName),'w') as fo:
            fo.write(povmeScript)
        povmeShellScript += '../../../../arun python ../../../POVME2.py %s\n' %(povmeScriptName)

        
        os.system('cp %s mdpocket/' %(target))
        os.system('cp %s dpocket/' %(complexTarget))
        with open('mdpocket/just_%s_pdbList' %(protID),'w') as fo:
            fo.write(pdbName)
        
        mdPocketShellScript += 'mkdir %s \nmdpocket -L just_%s_pdbList -f ../povme/%s/%s_point_field.pdb -o %s/mdpout\n' %(protID, protID, protID, protID, protID)
        complexPdbName = pdbName.replace('_pro','')
        ligandName = [line for line in open(complexTarget).readlines() if line[:4]=='ATOM'][-1].split()[3]
        with open('dpocket/just_%s_pdbList' %(protID),'w') as fo:
            fo.write('%s %s' %(complexPdbName, ligandName))
        dPocketShellScript +=  'mkdir %s \ndpocket -v 20000 -f just_%s_pdbList -o %s/dpout\n' %(protID, protID, protID)
        
        

with open('povme/runAll.sh','w') as fo:
    fo.write(povmeShellScript)
with open('mdpocket/runAll.sh','w') as fo:
    fo.write(mdPocketShellScript)
with open('dpocket/runAll.sh','w') as fo:
    fo.write(dPocketShellScript)
        
            

    
