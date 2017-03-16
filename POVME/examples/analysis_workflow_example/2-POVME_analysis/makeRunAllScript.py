
import glob
systemNames = [i.split('/')[-1].replace('_aligned.pdb','') for i in glob.glob('../1-trajectories/*_aligned.pdb')]


povmeInputTemplate = '''
PDBFileName	../1-trajectories/!!!SYSTEM NAME!!!_aligned.pdb
GridSpacing	1.0
#InclusionSphere	25.400 37.100 34.200 12.0
InclusionSphere		38.2 -47.5 63.7 15
#SeedSphere	25.400 37.100 34.200 4.0 
#SeedBox	27.472 19.570 40.038 6.0 6.0 6.0
#SeedBox	13.487 23.673 36.818 6.0 6.0 6.0
ContiguousPointsCriteria 3
DistanceCutoff	1.09
ConvexHullExclusion	none
OutputFilenamePrefix	./!!!SYSTEM NAME!!!/!!!SYSTEM NAME!!!_
CompressOutput	false
NumProcessors	5


'''

for systemName in systemNames:
    scriptText = povmeInputTemplate.replace('!!!SYSTEM NAME!!!',systemName)
    with open('povme_%s_input.ini' %(systemName), 'wb') as of:
        of.write(scriptText)
        
with open('runAll.sh', 'wb') as of:
    for systemName in systemNames:
        newLine = 'POVME3.py povme_%s_input.ini  \n' %(systemName)
        of.write(newLine)
