# trace generated using paraview version 5.9.0

#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# create a new 'OpenFOAMReader'
c0p04_Nb2_Lambda5_Re1e5foam = OpenFOAMReader(registrationName='c0p04_Nb2_Lambda5_Re1e5.foam', FileName='/media/clemenco1q/My Passport/19TRIBINE/simulation/performed_cases/bi_axiale_turbine/reynolds/c0p04_Nb2_Lambda5_Re1e5/c0p04_Nb2_Lambda5_Re1e5.foam')

# get active view
renderView1 = GetActiveViewOrCreate('RenderView')

# Properties modified on renderView1
renderView1.CameraParallelProjection = 1

# Properties modified on c0p04_Nb2_Lambda5_Re1e5foam
c0p04_Nb2_Lambda5_Re1e5foam.CaseType = 'Decomposed Case'

# show data in view
c0p04_Nb2_Lambda5_Re1e5foamDisplay = Show(c0p04_Nb2_Lambda5_Re1e5foam, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
c0p04_Nb2_Lambda5_Re1e5foamDisplay.Representation = 'Surface'

# reset view to fit data
renderView1.ResetCamera()

# show color bar/color legend
c0p04_Nb2_Lambda5_Re1e5foamDisplay.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'p'
pLUT = GetColorTransferFunction('p')

# get opacity transfer function/opacity map for 'p'
pPWF = GetOpacityTransferFunction('p')

# get animation scene
animationScene1 = GetAnimationScene()

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

# create a new 'Threshold'
threshold1 = Threshold(registrationName='Threshold1', Input=c0p04_Nb2_Lambda5_Re1e5foam)

# Properties modified on c0p04_Nb2_Lambda5_Re1e5foam
#c0p04_Nb2_Lambda5_Re1e5foam.CellArrays = ['H', 'HbyA', 'U', 'U_0', 'V0', 'cellTypes', 'k', 'k_0', 'nut', 'omega', 'omega_0', 'p', 'rAU', 'yPlus', 'yPlusExtremum', 'zoneID']

# Properties modified on threshold1
threshold1.Scalars = ['CELLS', 'cellTypes']
threshold1.ThresholdRange = [0.0, 1.0]

# show data in view
threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold1Display.Representation = 'Surface'

# hide data in view
Hide(c0p04_Nb2_Lambda5_Re1e5foam, renderView1)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# create a new 'Threshold'
threshold2 = Threshold(registrationName='Threshold2', Input=threshold1)

# Properties modified on threshold2
threshold2.Scalars = ['CELLS', 'zoneID']
threshold2.ThresholdRange = [1.0, 1000.0]

# show data in view
threshold2Display = Show(threshold2, renderView1, 'UnstructuredGridRepresentation')

# trace defaults for the display properties.
threshold2Display.Representation = 'Surface'

# hide data in view
Hide(threshold1, renderView1)

# show color bar/color legend
threshold2Display.SetScalarBarVisibility(renderView1, True)

# update the view to ensure updated data information
renderView1.Update()

# set active source
SetActiveSource(threshold1)

# show data in view
threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# Properties modified on animationScene1
animationScene1.AnimationTime = 28.54875756600115

# get the time-keeper
timeKeeper1 = GetTimeKeeper()

# set active source
SetActiveSource(threshold2)

# set active source
SetActiveSource(threshold1)

# set active source
SetActiveSource(c0p04_Nb2_Lambda5_Re1e5foam)

# set active source
SetActiveSource(threshold1)

# rescale color and/or opacity maps used to exactly fit the current data range
threshold1Display.RescaleTransferFunctionToDataRange(False, True)

# set scalar coloring
ColorBy(threshold1Display, ('POINTS', 'vorticity', 'Magnitude'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
threshold1Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# get color transfer function/color map for 'vorticity'
vorticityLUT = GetColorTransferFunction('vorticity')

# get opacity transfer function/opacity map for 'vorticity'
vorticityPWF = GetOpacityTransferFunction('vorticity')

# set scalar coloring
ColorBy(threshold1Display, ('POINTS', 'vorticity', 'Z'))

# rescale color and/or opacity maps used to exactly fit the current data range
threshold1Display.RescaleTransferFunctionToDataRange(False, False)

# Update a scalar bar component title.
UpdateScalarBarsComponentTitle(vorticityLUT, threshold1Display)

# Rescale transfer function
vorticityLUT.RescaleTransferFunction(-80.0, 80.0)

# Rescale transfer function
vorticityPWF.RescaleTransferFunction(-80.0, 80.0)

# set active source
SetActiveSource(threshold2)

# set scalar coloring
ColorBy(threshold2Display, ('POINTS', 'vorticity', 'Z'))

# Hide the scalar bar for this color map if no visible data is colored by it.
HideScalarBarIfNotNeeded(pLUT, renderView1)

# rescale color and/or opacity maps used to include current data range
threshold2Display.RescaleTransferFunctionToDataRange(True, False)

# show color bar/color legend
threshold2Display.SetScalarBarVisibility(renderView1, True)

# Rescale transfer function
vorticityLUT.RescaleTransferFunction(-80.0, 80.0)

# Rescale transfer function
vorticityPWF.RescaleTransferFunction(-80.0, 80.0)

# Properties modified on threshold2Display
threshold2Display.Ambient = 0.8

# Properties modified on threshold2Display
threshold2Display.Diffuse = 0.5

# set active source
SetActiveSource(threshold1)

# Properties modified on threshold1Display
threshold1Display.Ambient = 0.8

# Properties modified on threshold1Display
threshold1Display.Diffuse = 0.5

# hide color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, False)

# get layout
layout1 = GetLayout()

# layout/tab size in pixels
layout1.SetSize(842, 811)

# current camera placement for renderView1
renderView1.CameraPosition = [0.03235402075457012, -0.19480470887005935, 94.34480138967297]
renderView1.CameraFocalPoint = [0.03235402075457012, -0.19480470887005935, 0.0]
renderView1.CameraParallelScale = 0.44588352770268275
renderView1.CameraParallelProjection = 1

# save screenshot
SaveScreenshot('/media/clemenco1q/My Passport/19TRIBINE/simulation/performed_cases/bi_axiale_turbine/reynolds/c0p04_Nb2_Lambda5_Re1e5/voritcity2.png', renderView1, ImageResolution=[842, 811])

# layout/tab size in pixels
layout1.SetSize(842, 811)

# current camera placement for renderView1
renderView1.CameraPosition = [-0.044749293306441756, 0.0896561719980696, 94.34480138967297]
renderView1.CameraFocalPoint = [-0.044749293306441756, 0.0896561719980696, 0.0]
renderView1.CameraParallelScale = 0.08019607173305583
renderView1.CameraParallelProjection = 1

# save screenshot
SaveScreenshot('/media/clemenco1q/My Passport/19TRIBINE/simulation/performed_cases/bi_axiale_turbine/reynolds/c0p04_Nb2_Lambda5_Re1e5/voritcity2_zoom.png', renderView1, ImageResolution=[842, 811])

# set active source
SetActiveSource(threshold2)

# set active source
SetActiveSource(threshold1)

# hide data in view
Hide(threshold2, renderView1)

# show data in view
threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
threshold1Display.SetScalarBarVisibility(renderView1, True)

# destroy threshold2
Delete(threshold2)
del threshold2

# set active source
SetActiveSource(c0p04_Nb2_Lambda5_Re1e5foam)

# hide data in view
Hide(threshold1, renderView1)

# show data in view
c0p04_Nb2_Lambda5_Re1e5foamDisplay = Show(c0p04_Nb2_Lambda5_Re1e5foam, renderView1, 'UnstructuredGridRepresentation')

# show color bar/color legend
c0p04_Nb2_Lambda5_Re1e5foamDisplay.SetScalarBarVisibility(renderView1, True)

# destroy threshold1
Delete(threshold1)
del threshold1

# destroy c0p04_Nb2_Lambda5_Re1e5foam
Delete(c0p04_Nb2_Lambda5_Re1e5foam)
del c0p04_Nb2_Lambda5_Re1e5foam

# update animation scene based on data timesteps
animationScene1.UpdateAnimationUsingDataTimeSteps()

#================================================================
# addendum: following script captures some of the application
# state to faithfully reproduce the visualization during playback
#================================================================

#--------------------------------
# saving layout sizes for layouts

# layout/tab size in pixels
layout1.SetSize(842, 811)

#-----------------------------------
# saving camera placements for views

# current camera placement for renderView1
renderView1.CameraPosition = [-0.044749293306441756, 0.0896561719980696, 94.34480138967297]
renderView1.CameraFocalPoint = [-0.044749293306441756, 0.0896561719980696, 0.0]
renderView1.CameraParallelScale = 0.08019607173305579
renderView1.CameraParallelProjection = 1

#--------------------------------------------
# uncomment the following to render all views
# RenderAllViews()
# alternatively, if you want to write images, you can use SaveScreenshot(...).