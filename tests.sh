cd POVME/tests
#cd consistent_convex_hull_first
#../../arun python ../../POVME3.py sample_input.ini
#cd ../
python ./runRegTests.py
# Capture return code
RET=$?
cd ../../
# Pass return code of runRegTests.py to be the exit code for this script
$(exit $RET)