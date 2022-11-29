#!/usr/bin/env python



#### import the simple module from the paraview
from paraview.simple import *

import sys
import os
# By default, we try to read reference_cases on the cluster:
if os.path.isdir("/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axial_turbine/py/func"):
    sys.path.append("/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axial_turbine/py/func") 
    print('connected to cluster')
    isOnCluster = True
# if we are only in local, we read reference_cases locally:
else:
    sys.path.append("/home/users/clemenco1q/Documents/19TRIBINE_local_confinement2/simulation/19tribine_ref_cases/bi_axial_turbine/py/func")
    print('Lecture en local')
    isOnCluster = False

from postprocess_function_def import read_time_dir, \
                                     x_coord, \
                                     y_coord

from preprocess_function_def import y_coord, x_coord



from general_function_def import read_simulation_properties





def save_screenShots(foam_file_name, path_wd, path_to_foam, Time, Field, Scale, isVectorField, Screens, Ambient = 0.8, Diffuse = 0.5):
    """save a list of screenshot defined in Screens.
    
    inputs:
    - foam_file_name: str, '.foam' file name.
    - path_wd: str, path to OF case
    - path_to_foam: str, path to '.foam' file
    - Time: float, simulation time to plot
    - Field: str, name of the field
    - Scale: [min, max], values of the 
    - isVectorField: bool, True is the field is a VectorField
    - Screens: list of dict
    - Ambient: float, setting of lighting
    - Diffuse: float, setting of lighting"""
    #### disable automatic camera reset on 'Show'
    paraview.simple._DisableFirstRenderCameraReset()

    # create a new 'OpenFOAMReader'
    my_OF_reader = OpenFOAMReader(registrationName=foam_file_name, FileName=path_to_foam)

    # get active view
    renderView1 = GetActiveViewOrCreate('RenderView')

    # Properties modified on renderView1
    renderView1.CameraParallelProjection = 1

    # Properties modified on my_OF_reader
    my_OF_reader.CaseType = 'Decomposed Case'

    # show data in view
    my_OF_readerDisplay = Show(my_OF_reader, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    my_OF_readerDisplay.Representation = 'Surface'

    # reset view to fit data
    renderView1.ResetCamera()

    # show color bar/color legend
    my_OF_readerDisplay.SetScalarBarVisibility(renderView1, False)

    # get animation scene
    animationScene1 = GetAnimationScene()

    # update animation scene based on data timesteps
    animationScene1.UpdateAnimationUsingDataTimeSteps()

    # create a new 'Threshold'
    threshold1 = Threshold(registrationName='Threshold1', Input=my_OF_reader)

    # Properties modified on threshold1
    threshold1.Scalars = ['CELLS', 'cellTypes']
    threshold1.ThresholdRange = [0.0, 1.0]

    # show data in view
    threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

    # trace defaults for the display properties.
    threshold1Display.Representation = 'Surface'

    # hide data in view
    Hide(my_OF_reader, renderView1)

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

    # update the view to ensure updated data information
    renderView1.Update()

    # set active source
    SetActiveSource(threshold1)

    # show data in view
    threshold1Display = Show(threshold1, renderView1, 'UnstructuredGridRepresentation')

    # Properties modified on animationScene1
    animationScene1.AnimationTime = Time

    # get the time-keeper
    timeKeeper1 = GetTimeKeeper()

    # set active source
    SetActiveSource(threshold2)

    # set active source
    SetActiveSource(threshold1)

    # set active source
    SetActiveSource(my_OF_reader)

    # set active source
    SetActiveSource(threshold1)

    # set scalar coloring
    if isVectorField:
        ColorBy(threshold1Display, ('POINTS', Field[0],  Field[1]))
    else:
        ColorBy(threshold1Display, ('POINTS', Field[0]))

    # get color transfer function/color map for Field
    myfieldLUT = GetColorTransferFunction(Field[0])

    # get opacity transfer function/opacity map for Field
    myfieldPWF = GetOpacityTransferFunction(Field[0])

    # Rescale transfer function
    myfieldLUT.RescaleTransferFunction(Scale[0], Scale[1])

    # Rescale transfer function
    myfieldPWF.RescaleTransferFunction(Scale[0], Scale[1])

    # set active source
    SetActiveSource(threshold2)

    # set scalar coloring
    if isVectorField:
        ColorBy(threshold2Display, ('POINTS', Field[0],  Field[1]))
    else:
        ColorBy(threshold2Display, ('POINTS', Field[0]))
        

    # Rescale transfer function
    myfieldLUT.RescaleTransferFunction(Scale[0], Scale[1])

    # Rescale transfer function
    myfieldPWF.RescaleTransferFunction(Scale[0], Scale[1])

    # Properties modified on threshold2Display
    threshold2Display.Ambient = Ambient

    # Properties modified on threshold2Display
    threshold2Display.Diffuse = Diffuse

    # set active source
    SetActiveSource(threshold1)

    # Properties modified on threshold1Display
    threshold1Display.Ambient = Ambient

    # Properties modified on threshold1Display
    threshold1Display.Diffuse = Diffuse

    # hide color bar/color legend
    threshold1Display.SetScalarBarVisibility(renderView1, False)
    
    
    # hide color bar/color legend
    threshold2Display.SetScalarBarVisibility(renderView1, False)
    
    # get layout
    layout1 = GetLayout()



                    
                    
    for screen in Screens:
        # boudin box
        xmin, xmax = screen['bbox'][0][0], screen['bbox'][0][1]
        ymin, ymax = screen['bbox'][1][0], screen['bbox'][1][1]

        # layout/tab size in pixels
        image_height = ymax - ymin
        image_width = xmax - xmin
        pixel_width = screen['pixel_width']
        pixel_height = int(pixel_width * image_height / image_width)

        layout1.SetSize(pixel_width, pixel_height)

        # compute camera scale:
        cameraScale = 0.5 * max(image_width, image_height)

        # current camera placement for renderView1

        renderView1.CameraPosition = [0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 1]
        renderView1.CameraFocalPoint = [0.5 * (xmin + xmax), 0.5 * (ymin + ymax), 0]
        renderView1.CameraParallelScale = cameraScale
        renderView1.CameraParallelProjection = 1

        # save screenshot
        SaveScreenshot(screen['save_path'], renderView1, ImageResolution=[pixel_width, pixel_height])



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
    SetActiveSource(my_OF_reader)

    # hide data in view
    Hide(threshold1, renderView1)

    # show data in view
    my_OF_readerDisplay = Show(my_OF_reader, renderView1, 'UnstructuredGridRepresentation')

    # show color bar/color legend
    my_OF_readerDisplay.SetScalarBarVisibility(renderView1, True)

    # destroy threshold1
    Delete(threshold1)
    del threshold1

    # destroy my_OF_reader
    Delete(my_OF_reader)
    del my_OF_reader






        
