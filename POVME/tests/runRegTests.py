import glob
import contextlib
import os
#import filecmp
import difflib
import re
import subprocess
import time
import argparse
import sys


@contextlib.contextmanager
def Chdir(directory):
    cwd = os.getcwd()
    os.chdir(directory)
    yield
    os.chdir(cwd)

rawRegTests = open('regTests').readlines()
#Ignore commented out lines
rawRegTests = [i for i in rawRegTests if not i.strip()[0]=='#']


regTests = []
results = {}
arunPath = os.getcwd()+'/../arun'


regexFileName = "file_comparison_ignore_regex"
regexes_to_ignore = [i.strip() for i in open(regexFileName).readlines()]

def remove_regex_lines(list_of_lines):
    lines_to_remove = set([])
    for i, line in enumerate(list_of_lines):
        for regex in regexes_to_ignore:
            if re.search(regex,line) != None:
                lines_to_remove.add(i)
    for i in sorted(list(lines_to_remove), reverse=True):
        list_of_lines.pop(i)
    return list_of_lines


def compareFile(origFile, args):
    origFileData = open(origFile).readlines()
    origFileData = remove_regex_lines(origFileData)
    newFile = origFile.replace('.orig','')
    if not os.path.exists(newFile):
        print('File %s does not exist - Unable to perform file comparison' %(newFile))
        passed = False
        return passed
    newFileData = open(newFile).readlines()
    newFileData = remove_regex_lines(newFileData)
    #diff = difflib.compare(newFileData, origFileData)
    #if difflib.ndiff(newFileData, origFileData)
    #if filecmp.cmp(newFile, origFile) == True:
    if newFileData == origFileData:
        #passedFiles.append(origFile)
        passed = True
        print("Files %s and %s match!" %(newFile, origFile))
    else:
        #failedFiles.append(origFile)
        passed = False
        print('File %s DOES NOT MATCH' %(origFile))
        if args.compare == True:
            validChoice = False
            while validChoice == False:
                choice = input('Files %s and %s differ. View differences (y,n,v)? ' %(newFile, origFile))
                if choice == 'y':
                    ignoreStr = ' '.join(['-I %s' %(i) for i in regexes_to_ignore])
                    os.system('diff %s %s %s | less' %(ignoreStr, newFile, origFile))
                    validChoice = True
                elif choice == 'v':
                    os.system('tkdiff %s %s' %(newFile, origFile))
                    validChoice = True
                elif choice == 'n':
                    validChoice = True
    return passed

def runTests(args):
    for line in rawRegTests:
        linesp = line.split()
        title = linesp[0]
        results[title]={}
        directory = os.getcwd()+'/'+linesp[1]+'/'
        script = ' '.join(linesp[2:])
        regTests.append([title, directory, script])
        passedFiles = []
        failedFiles = []

        with Chdir(directory):
            print() 
            origFiles = glob.glob('*orig')
            if args.remove_old == True:
                files_to_remove = [i.replace('.orig','') for i in origFiles]
                for file_to_remove in files_to_remove:
                    print('Removing', file_to_remove)
                    os.system('rm %s' %(file_to_remove))
            print("RUNNING TEST %s" %(linesp))
            runCommand = '%s python %s > output' %(arunPath, script)
            print("Run command: %s" %(runCommand))
            #os.system('%s python %s > output' %(arunPath,script))
            #print '%s python %s/%s' %(arunPath, directory, script)
            start = time.time()
            p=subprocess.Popen(runCommand,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               cwd = directory, shell=True)

            out, error = p.communicate()
            print('Test ran in %s seconds' %(time.time()-start))
            #print out, error

            results[title]['testPassed'] = True
            if p.returncode == 0:
                print("Exit status: Completed")
                results[title]['exit_status'] = True
            else:
                print("Exit status: Failed")
                print('error %s' %(error))
                results[title]['exit_status'] = False
                results[title]['testPassed'] = False

            for origFile in origFiles:
                passed = compareFile(origFile, args)
                if passed == True:
                    passedFiles.append(origFile)
                else:
                    failedFiles.append(origFile)
                    results[title]['testPassed'] = False
        results[title]['passedFiles'] = passedFiles
        results[title]['failedFiles'] = failedFiles



    nPassed = 0
    nFailed = 0
    print()
    print('===============')
    print("RESULTS SUMMARY")
    print('===============')
    for test in regTests:
        print("----Test %s   Exit Status = %s ----" %(test[0], results[test[0]]['exit_status']))
        print("%i file comparisons succeeded: %s" %(len(results[test[0]]['passedFiles']),
                                                    results[test[0]]['passedFiles']))
        print("%i file comparisons failed: %s" %(len(results[test[0]]['failedFiles']), results[test[0]]['failedFiles']))
        if results[test[0]]['testPassed'] == True:
            nPassed += 1
        else:
            nFailed += 1
    print()
    print('%i tests passed, %i tests failed' %(nPassed, nFailed))
    if nFailed > 0:
        print("TEST SUITE FAILED")
        sys.exit(1)


        
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run regression tests for this module")
    parser.add_argument('-c', '--compare', 
                        help='Offer comparison options if file fails diff',
                        action='store_true')
    parser.add_argument('-r', '--remove_old', 
                        help='Remove existing files with .orig complements before running',
                        action='store_true')

    args = parser.parse_args()
    
    runTests(args)
