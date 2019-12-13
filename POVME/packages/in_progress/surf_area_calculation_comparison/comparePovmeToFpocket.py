import glob
import numpy
import pylab
import random
import copy
targets = glob.glob('povme/*')
#print targets

normalize = True
allAsSubplots = True
comparePovmeTo = 'dpocket'
#comparePovmeTo = 'mdpocket'

if comparePovmeTo == 'mdpocket':
    targets = [i.split('/')[1] for i in targets if (not ('.' in i) and not('rot' in i))]
if comparePovmeTo == 'dpocket':
    targets = [i.split('/')[1] for i in targets if (not ('.' in i) and not('rot' in i) and (not 'APO' in i))]
#targets = [i.split('/')[1] for i in targets if (not ('.' in i))]


#comparePovmeTo = 'mdpocket'
#plotDifference = True

povmeVolumes = dict((target,[]) for target in targets)
povmeSAs = dict((target,[]) for target in targets)
povmePols = dict((target,[]) for target in targets)
mdpocketVolumes = dict((target,[]) for target in targets)
mdpocketSAs = dict((target,[]) for target in targets)
mdpocketPols = dict((target,[]) for target in targets)
#toPlot = []
povmeVolAvList = []
povmeSAAvList = []
povmePolAvList = []
mdpocketVolAvList = []
mdpocketSAAvList = []
mdpocketPolAvList = []


povmeVolStdList = []
povmeSAStdList = []
povmePolStdList = []
mdpocketVolStdList = []
mdpocketSAStdList = []
mdpocketPolStdList = []
toPlotColors = []
toPlotLabels = []


toPlotVolStds = []
toPlotSAStds = []
toPlotPolStds = []



for target in targets:
    for rotation in ['','_rot1','_rot2','_rot3','_rot4','_rot5']:
        povmeVolData = open('povme/%s%s/%s%s_frame_1.pdb' %(target, rotation, target, rotation)).readlines()
        povmeVolume = float(povmeVolData[1].split()[3])
        povmeVolumes[target].append(povmeVolume)
        #povmeVolList.append(povmeVolume)
        povmeSAData = open('povme/%s%s/%s%s_frame_1_surface.pdb' %(target, rotation, target, rotation)).readlines()
        povmeSA = float(povmeSAData[1].split()[5])
        povmeSAs[target].append(povmeSA)
        #povmeSAList.append(povmeSA)

        povmePol = float(povmeSAData[2].split()[4])
        povmePols[target].append(povmePol)
        #povmePolList.append(povmePol)

        if comparePovmeTo == 'mdpocket':
            mdpocketData = open('mdpocket/%s%s/mdpout_descriptors.txt' %(target, rotation)).readlines()
            volumeIndex = mdpocketData[0].split().index('pock_volume')
            mdpocketVolume = float(mdpocketData[1].split()[volumeIndex])
            SAIndex = mdpocketData[0].split().index('pock_asa')
            mdpocketSA = float(mdpocketData[1].split()[SAIndex])
            hydIndex = mdpocketData[0].split().index('hydrophobicity_score') # NOT USED
            polIndex = mdpocketData[0].split().index('polarity_score')
            mdpocketPol = float(mdpocketData[1].split()[polIndex])
            
            mdpocketVolumes[target].append(mdpocketVolume)
            mdpocketSAs[target].append(mdpocketSA)
            mdpocketPols[target].append(mdpocketPol)
        if comparePovmeTo == 'dpocket':
            ## Uncomment this to use the unbiased fpocket output
            mdpocketData = open('dpocket/%s%s/dpout_fp.txt' %(target, rotation)).readlines()
            ## Uncomment this to use the dpocket ligand-based pocket output

            #mdpocketData = open('dpocket/%s%s/dpout_exp.txt' %(target, rotation)).readlines()
            if len(mdpocketData) == 1: # If there are no pockets detected
                print("No data for %s %s. Skipping!" %(target, rotation))
                continue
            else:
                volumeIndex = mdpocketData[0].split().index('pock_vol')
                mdpocketVolume = float(mdpocketData[1].split()[volumeIndex])
                #SAIndex = mdpocketData[0].split().index('pock_asa')
                #mdpocketSA = float(mdpocketData[1].split()[SAIndex])
                mdpocketSA = 0.
                hydIndex = mdpocketData[0].split().index('hydrophobicity_score') # NOT USED
                polIndex = mdpocketData[0].split().index('polarity_score_norm')
                mdpocketPol = float(mdpocketData[1].split()[hydIndex])
            
                mdpocketVolumes[target].append(mdpocketVolume)
                mdpocketSAs[target].append(mdpocketSA)
                mdpocketPols[target].append(mdpocketPol)
    povmeVolAvList.append(numpy.mean(povmeVolumes[target]))
    povmeSAAvList.append(numpy.mean(povmeSAs[target]))
    povmePolAvList.append(numpy.mean(povmePols[target]))
    
    mdpocketVolAvList.append(numpy.mean(mdpocketVolumes[target]))
    mdpocketSAAvList.append(numpy.mean(mdpocketSAs[target]))
    mdpocketPolAvList.append(numpy.mean(mdpocketPols[target]))
    
    
    povmeVolStdList.append(numpy.std(povmeVolumes[target]))
    povmeSAStdList.append(numpy.std(povmeSAs[target]))
    povmePolStdList.append(numpy.std(povmePols[target]))
    
    mdpocketVolStdList.append(numpy.std(mdpocketVolumes[target]))
    mdpocketSAStdList.append(numpy.std(mdpocketSAs[target]))
    mdpocketPolStdList.append(numpy.std(mdpocketPols[target]))
    
    #toPlotVolStds = toPlotVolStds + [numpy.std(povmeVolumes[target]), numpy.std(mdpocketVolumes[target]),0]
    #toPlotSAStds = toPlotSAStds + [numpy.std(povmeSAs[target]), numpy.std(mdpocketSAs[target]),0]
    #toPlotPolStds = toPlotPolStds + [numpy.std(povmePols[target]), numpy.std(mdpocketPols[target]),0]
    
    #toPlot = toPlot + [povmeVolume, mdpocketVolume,0.0]
    toPlotColors = toPlotColors + [[0,0,1],[0,1,0],[1,0,0]]
    #toPlotLabels = toPlotLabels + ['povme_'+target, 'mdpocket_'+target,'']
    toPlotLabels = toPlotLabels + ['',target.replace('_pro',''),'']
    #print toPlot

def analyze_difference(list1, list2, stdList1 = None, stdList2 = None, normalize=True, difference=True):
    povmeVolumeArray = numpy.array(list1)
    mdpocketVolumeArray = numpy.array(list2)
    if stdList1 == None:
        povmeStdArray = numpy.zeros((len(list1)))
    else:
        povmeStdArray = numpy.array(stdList1)
    if stdList2 == None:
        mdpocketStdArray = numpy.zeros((len(list2)))
    else:
        mdpocketStdArray = numpy.array(stdList2)
        
    if normalize:
        povmeNormConst = numpy.mean(povmeVolumeArray)
        povmeVolumeArray /= povmeNormConst
        mdpocketNormConst = numpy.mean(mdpocketVolumeArray)
        mdpocketVolumeArray /= mdpocketNormConst
        povmeStdArray /= povmeNormConst
        mdpocketStdArray /= mdpocketNormConst

    toPlot = list(zip(povmeVolumeArray,mdpocketVolumeArray, numpy.zeros((len(targets)))))
    toPlot = [i for t in toPlot for i in t]
    if stdList1 != None and stdList2 != None:
        toPlotStds = list(zip(povmeStdArray,mdpocketStdArray, numpy.zeros((len(targets))) ))
        toPlotStds = [i for t in toPlotStds for i in t]
    
    differences = []
    #if difference:
    for i in range(len(toPlot)/3):
        thisDifference = abs(toPlot[i*3] - toPlot[1+(i*3)])
        # Or, if we want to dot-product the two vectors
        #thisDifference = abs(toPlot[i*3] * toPlot[1+(i*3)])
        toPlot[2+(i*3)] = thisDifference
        differences.append(thisDifference)
    if stdList1 != None and stdList2 != None:
        return toPlot, toPlotStds, differences
    else:
        return toPlot, differences
    #else:
    #    if stdList1 != None and stdList2 != None:
    #        return toPlot, toPlotStds
    #    else:
    #        return toPlot


def bootstrap_analysis_sum_difference(list1, list2):
    nBootstraps = 30000

    fpvl = copy.deepcopy(list1)
    pvl = numpy.array(copy.deepcopy(list2))
    sumDifferences = []
    for i in range(nBootstraps):
        newIndices = numpy.random.randint(0,len(pvl),(len(pvl)))
        newPvl = pvl[newIndices]
        #random.shuffle(pvl)
        toPlot, difference = analyze_difference(newPvl,fpvl, normalize=normalize)
        sumDifferences.append(sum(difference))
    hist = pylab.hist(sumDifferences, bins=50)
    realDifference = sum(realDifferences)

    betterPop = 0
    for pop,border in zip(hist[0], hist[1][1:]):
        if border > realDifference:
            betterPop += pop
    
    pylab.axvline(realDifference)
    betterThanPercent = 100*float(betterPop)/nBootstraps
    pylab.text(realDifference, 0.5*max(hist[0]),'%.1g' %(betterThanPercent))
    print('Real difference is better than %r percent of bootstrap permutations' %(betterThanPercent))

def linear_regression(xData, yData, xErr = None, yErr = None, fit=True, plot=True):
    
    
    pylab.scatter(xData, yData)
    pylab.errorbar(xData, yData, xerr=xErr, yerr=yErr, linestyle='None')
    pylab.xlabel('POVME')
    pylab.ylabel('mdpocket')       
    if fit:
        coefs = numpy.lib.polyfit(xData, yData, deg=1, full=True)
        #print self.coefs
        fit_y = numpy.lib.polyval(coefs[0], xData)
              
        #calculate R^2 (from http://stackoverflow.com/questions/893657/
        #how-do-i-calculate-r-squared-using-python-and-numpy)
        #p=numpy.poly1d(self.coefs[0])
        #yhat = p(xData)
        ybar = numpy.sum(yData) / len(yData)
        ssreg = numpy.sum((fit_y-ybar)**2)
        sstot = numpy.sum((yData - ybar)**2)
        Rsq = ssreg/sstot
        if plot:
            pylab.plot(xData, fit_y)      
            pylab.text(xData[-1], fit_y[-1],
                       'slope = %.3f\nRsq = %.3f' %(coefs[0][0], Rsq))
            
    return Rsq
        
def subplotOrShow(subplotInd):#, title='', xLabel='', yLabel=''):
    #pylab.title = title
    #pylab.xaxis = xLabel
    #pylab.yaxis = yLabel
    if allAsSubplots:
        pylab.subplot(subplotInd)
    else:
        pylab.show()

print('ANALYZING POCKET VOLUMES')

if allAsSubplots:
    pylab.subplot(331)
toPlotVols, bar_toPlotVolStds, realDifferences = analyze_difference(povmeVolAvList, mdpocketVolAvList, stdList1=povmeVolStdList, stdList2=mdpocketVolStdList, normalize=normalize)
#print len(toPlotVolStds)
pylab.bar(list(range(len(toPlotVols))), toPlotVols, color=toPlotColors,yerr=bar_toPlotVolStds)
pylab.xticks(list(range(len(toPlotVols))),toPlotLabels,rotation = 90)
pylab.title('Volume Comparison')
pylab.ylabel('Volume'+' (normalized)'*normalize)

subplotOrShow(332)
linear_regression(povmeVolAvList, mdpocketVolAvList, xErr=povmeVolStdList, yErr = mdpocketVolStdList)
pylab.title('Volume Comparison')

subplotOrShow(333)
bootstrap_analysis_sum_difference(povmeVolAvList, mdpocketVolAvList)
pylab.title('Bootstrap resampling summed error')
pylab.xlabel('summed error')
pylab.ylabel('frequency')


if sum(mdpocketSAAvList) != 0:
    print('ANALYZING SURFACE AREAS')
    subplotOrShow(334)
    toPlotSAs, bar_toPlotSAStds, realDifferences = analyze_difference(povmeSAAvList, mdpocketSAAvList,stdList1=povmeSAStdList, stdList2=mdpocketSAStdList, normalize=normalize)

    pylab.bar(list(range(len(toPlotSAs))), toPlotSAs, color=toPlotColors, yerr=bar_toPlotSAStds)
    pylab.xticks(list(range(len(toPlotSAs))),toPlotLabels,rotation = 90)
    pylab.title('Surface Area Comparison')
    pylab.ylabel('Surf. Area'+' (normalized)'*normalize)

    subplotOrShow(335)
    linear_regression(povmeSAAvList, mdpocketSAAvList, xErr=povmeSAStdList, yErr = mdpocketSAStdList)
    pylab.title('Surf. Area Comparison')

    subplotOrShow(336)
    bootstrap_analysis_sum_difference(povmeSAAvList, mdpocketSAAvList)
    pylab.title('Bootstrap resampling summed error')
    pylab.xlabel('summed error')
    pylab.ylabel('frequency')


print('ANALYZING POLARITY SCORES')
subplotOrShow(337)
toPlotPols, bar_toPlotPolStds, realDifferences = analyze_difference(povmePolAvList, mdpocketPolAvList,stdList1=povmePolStdList, stdList2=mdpocketPolStdList, normalize=normalize)
pylab.bar(list(range(len(toPlotPols))), toPlotPols, color=toPlotColors, yerr=bar_toPlotPolStds)
pylab.xticks(list(range(len(toPlotPols))),toPlotLabels,rotation = 90)
pylab.title('Polarity Comparison')
pylab.ylabel('Polarity'+' (normalized)'*normalize)

subplotOrShow(338)
linear_regression(povmePolAvList, mdpocketPolAvList)
pylab.title('Polarity Comparison')

subplotOrShow(339)
bootstrap_analysis_sum_difference(povmePolAvList, mdpocketPolAvList)
pylab.title('Bootstrap resampling summed error')
pylab.xlabel('summed error')
pylab.ylabel('frequency')

if allAsSubplots:
    pylab.show()
#pylab.hist(differences)
#pylab.show()bootstrap = True
