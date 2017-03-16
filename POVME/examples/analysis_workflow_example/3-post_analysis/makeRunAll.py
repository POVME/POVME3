import os
import shutil
import glob


systemNames = [i.split('/')[-1].replace('_aligned.pdb','') for i in glob.glob('../1-trajectories/*_aligned.pdb')]
dashTList = ' '.join(['-t !!!SYSTEM NAME!!!_:../../1-trajectories/!!!SYSTEM NAME!!!_aligned.pdb'.replace('!!!SYSTEM NAME!!!',systemName) for systemName in systemNames])

runAllText = []


systemName = 'ALL'
if os.path.exists(systemName):
    shutil.rmtree(systemName)
os.mkdir(systemName)
ALL_frameList = ['../../2-POVME_analysis/%s/%s_frameInfo/%s_frame_*[0-9].npy' %(i,i,i) for i in systemNames]
runAllScript = '''

binding_site_overlap.py -f !!ALL FRAME LIST!!

cluster.py -n 5 -m tanimoto_matrix.npy -i indexMapToFrames.csv !!DASH T LIST!! 
'''
runAllScript = runAllScript.replace('!!ALL FRAME LIST!!',' '.join(ALL_frameList))
runAllScript = runAllScript.replace('!!DASH T LIST!!', dashTList)
with open('ALL/runClustering_ALL.sh','wb') as of:
    of.write(runAllScript)

with open('ALL/bsubRunScript.sh','wb') as of:
    of.write('''#ln -s ../../../../packages/clustering/binding_site_overlap.py 
#ln -s ../../../../packages/clustering/cluster.py
. ./runClustering_ALL.sh 
''')
    
    
