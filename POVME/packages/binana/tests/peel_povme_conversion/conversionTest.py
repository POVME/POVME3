import numpy
import POVME.packages.binana.peel as peel
import POVME3

### Loads a numpy map from POVME output, converts it to a featureMap,
### then converts it back and sees if the data is the same

# Load the data
rawData = numpy.load('POVME_frame0.npy')
# Convert it to all_pts format
rawDict = {}
for point in rawData:
    rawDict[tuple(point)] = 1
minX = numpy.min(rawData[:,0])
maxX = numpy.max(rawData[:,0])
minY = numpy.min(rawData[:,1])
maxY = numpy.max(rawData[:,1])
minZ = numpy.min(rawData[:,2])
maxZ = numpy.max(rawData[:,2])
resolution = 1.0
xDim = int(numpy.round(maxX-minX))
yDim = int(numpy.round(maxY-minY))
zDim = int(numpy.round(maxZ-minZ))
data = numpy.zeros((xDim * yDim * zDim, 4))
c = 0
for x in numpy.arange(minX, maxX):
    for y in numpy.arange(minY, maxY):
        for z in numpy.arange(minZ, maxZ):
            data[c,0] = x
            data[c,1] = y
            data[c,2] = z
            data[c,3] = rawDict.get((x,y,z), 0.0)
            c += 1

#print data

# Have POVME's output function write it
povmeOutputName = 'POVME_generated.dx'
POVME3.dx_freq(data, {'OutputFilenamePrefix':'test_',
                      'SaveVolumetricDensityDX':True,
                      'CompressOutput':False})

# Convert it to a featuremap and have that write it
featureMap = peel.featureMap.fromPovmeList(data, 1.0)
featureMap.write_dx_file('peel_generated.dx')

# Convert it back to a POVME list and write it
data_2 = featureMap.toPovmeList()
POVME3.dx_freq(data, {'OutputFilenamePrefix':'test2_',
                      'SaveVolumetricDensityDX':True,
                      'CompressOutput':False})


# Compare outputs
dxFile_povme1 = open('test_volumetric_density.dx').read()
dxFile_peel = open('peel_generated.dx').read()
dxFile_povme2 = open('test2_volumetric_density.dx').read()

if dxFile_povme1 != dxFile_peel:
    print("dxFile_povme1 differs from dxFile_peel")
if dxFile_povme2 != dxFile_peel:
    print("dxFile_povme2 differs from dxFile_peel")
if dxFile_povme1 != dxFile_povme2:
    print("dxFile_povme1 differs from dxFile_povme2")
