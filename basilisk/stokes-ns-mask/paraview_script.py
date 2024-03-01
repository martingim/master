# trace generated using paraview version 5.12.0-RC3
#import paraview
#paraview.compatibility.major = 5
#paraview.compatibility.minor = 12

# to automatically get the filenames
import glob
from os import path
#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()


Filenames = [path.abspath(name) for name in glob.glob('vtu/TIME*.vtu')]
Filenames.sort()

# create a new 'XML Unstructured Grid Reader'
tIME00 = XMLUnstructuredGridReader(registrationName='TIME-00*', FileName=Filenames)
tIME00.CellArrayStatus = ['f', 'p', 'u.x']

# Properties modified on tIME00
tIME00.TimeArray = 'None'
vTU1 = CreateExtractor('VTU', tIME00, registrationName='VTU1')
# trace defaults for the extractor.
vTU1.Trigger = 'TimeStep'

# Properties modified on vTU1.Writer
vTU1.Writer.FileName = '{timestep:06d}.pvtu'
vTU1.Writer.DataMode = 'Ascii'
vTU1.Writer.UseSubdirectory = 0

# save extracts
SaveExtracts(ExtractsOutputDirectory='/home/martin/Documents/master/basilisk/stokes-ns-mask/vtu/ascii',
    GenerateCinemaSpecification=0)
