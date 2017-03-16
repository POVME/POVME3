import POVME.packages.binana.peel as peel
## Tests to make sure the file conversion functionality works for reading/writing dx files
my_fm = peel.featureMap.fromDxFile('exampleInput.dx')
my_fm.write_dx_file('exampleOutput.dx')
