
cd 2-POVME_analysis
python makeRunAllScript.py
. runAll.sh
#wait
cd ../

cd 3-post_analysis
python makeRunAll.py
cd ALL
. bsubRunScript.sh
python ../analyzeClusterMembershipInGroups.py
pocketPointsPca.py
cd ../
cd ../

