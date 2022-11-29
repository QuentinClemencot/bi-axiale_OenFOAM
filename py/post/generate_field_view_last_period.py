#!/usr/bin/env python

'''the purpose of this script is to save set of png representing field of last period
of revolution of a bi-axiale turbine.
inputs:
- doVorticity, doPressure, doUx: bools, put 'True' to save vorticity, pressure and Ux'''


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


from pv_function_def import save_screenShots


# read file call's arguments:

doVorticity = sys.argv[1]  # bool, True if you want vorticity field to be sreenshot
doPressure = sys.argv[2]    # bool, idem for pressure
doUx = sys.argv[3]          # bool, idem for velocity

# Read screen shot parameter:
# scene view settings:
Ambient = read_simulation_properties('Ambient', path = '.', file = 'fieldViewParameter.py')
Diffuse = read_simulation_properties('Diffuse', path = '.', file = 'fieldViewParameter.py')

# scale fo field:
Scale_vorticity = read_simulation_properties('Scale_vorticity', path = '.', file = 'fieldViewParameter.py')
Scale_p = read_simulation_properties('Scale_p', path = '.', file = 'fieldViewParameter.py')
Scale_Ux = read_simulation_properties('Scale_Ux', path = '.', file = 'fieldViewParameter.py')

# field to read:
isVectorField_vorticity = read_simulation_properties('Ambient', path = '.', file = 'fieldViewParameter.py')
isVectorField_p = read_simulation_properties('isVectorField_p', path = '.', file = 'fieldViewParameter.py')
isVectorField_Ux = read_simulation_properties('isVectorField_Ux', path = '.', file = 'fieldViewParameter.py')
Field_vorticity = read_simulation_properties('Field_vorticity', path = '.', file = 'fieldViewParameter.py')
Field_p = read_simulation_properties('Field_p', path = '.', file = 'fieldViewParameter.py')
Field_Ux = read_simulation_properties('Field_Ux', path = '.', file = 'fieldViewParameter.py')

# image definition
pixel_width = read_simulation_properties('pixel_width', path = '.', file = 'fieldViewParameter.py')

# image demention
rotor_height = read_simulation_properties('rotor_height', path = '.', file = 'fieldViewParameter.py')
rotor_width  = read_simulation_properties('rotor_width', path = '.', file = 'fieldViewParameter.py')
zoom_height  = read_simulation_properties('zoom_height', path = '.', file = 'fieldViewParameter.py')
zoom_width  = read_simulation_properties('zoom_width', path = '.', file = 'fieldViewParameter.py')

# read turbine parameters:
c = read_simulation_properties('c')
h = read_simulation_properties('h')
r = read_simulation_properties('r')
T = read_simulation_properties('T')
d = h + 2 * r


# define path
path_wd = '.'
foam_file_name = os.path.split(os.getcwd())[1] + '.foam' # current directory.foam file for pv opening
path_to_foam = os.path.join(path_wd, foam_file_name)

# read time dir in processor0
list_time = read_time_dir(os.path.join(path_wd, 'processor0'))

# keep only the times of the last period:
time_max = max([float(t) for t in list_time])
time_limit = time_max - T
list_time = [t for t in list_time if float(t) > time_limit]

# compute bounding box of rotor
xmin_rotor, xmax_rotor = -rotor_width * d / 2, rotor_width * d / 2
top_space = (rotor_height - 1) * d / 2 # space between rotor top and bbox top
										# same as bottom_space                    
ymin_rotor, ymax_rotor = -h -r - top_space, top_space + r


for time in list_time:
    # compute s parameter:
    timeAdim = float(time) / T
    s = timeAdim % 1
    # compute coordonate of A point:
    x, y = x_coord(s = s, h = h, r = r), y_coord(s = s, h = h, r = r)

    # compute bounding box of zoom on blade
    xmin_blade, xmax_blade = - zoom_width * c / 2 + x, zoom_width * c / 2 + x
    ymin_blade, ymax_blade = - zoom_height * c / 2 + y, zoom_height * c / 2 + y 

    print(f"Create field png for Time: {time}s")
    # time for png name
    time_for_saving = str(round(float(time), 8))
    time_for_saving = time_for_saving.replace(".", "p")
    
    # s for png name
    s_for_saving = str(round(float(s), 6))
    s_for_saving = s_for_saving.replace(".", "p")

    if doVorticity:
        # save vorticity images:
        save_screenShots(foam_file_name = foam_file_name,
                         path_wd = path_wd,
                         path_to_foam = path_to_foam,
                         Time = float(time), 
                         Field = Field_vorticity, 
                         Scale = Scale_vorticity, 
                         isVectorField = isVectorField_vorticity, 
                         #Screens = [{'save_path': os.path.join(path_to_save, 'vorticity_' + time_for_saving + '_rotor.png'),
                         Screens = [{'save_path': './images/vorticity_s'+s_for_saving+'_t'+time_for_saving+'_rotor.png',
                                     'pixel_width': pixel_width,
                                     'bbox': ((xmin_rotor, xmax_rotor ),(ymin_rotor, ymax_rotor))}, 
                                    
                                    {'save_path': './images/vorticity_s'+s_for_saving+'_t'+time_for_saving+'_blade.png',
                                     'pixel_width': pixel_width,
                                     'bbox': ((xmin_blade, xmax_blade ),(ymin_blade, ymax_blade))} ],
                          Ambient = Ambient,
                          Diffuse = Diffuse)
             
    if doPressure:
        # save pressure images:
        save_screenShots(foam_file_name = foam_file_name,
                         path_wd = path_wd,
                         path_to_foam = path_to_foam,
                         Time = float(time), 
                         Field = Field_p, 
                         Scale = Scale_p, 
                         isVectorField = isVectorField_p, 
                         Screens = [{'save_path': './images/p_s'+s_for_saving+'_t'+time_for_saving+'_rotor.png',
                                     'pixel_width': pixel_width,
                                     'bbox': ((xmin_rotor, xmax_rotor ),(ymin_rotor, ymax_rotor))}, 
                                    
                                    {'save_path': './images/p_s'+s_for_saving+'_t'+time_for_saving+'_blade.png',
                                     'pixel_width': pixel_width,
                                     'bbox': ((xmin_blade, xmax_blade ),(ymin_blade, ymax_blade))} ],
                          Ambient = Ambient,
                          Diffuse = Diffuse)

    if doUx:
        # save Ux images:
        save_screenShots(foam_file_name = foam_file_name,
                         path_wd = path_wd,
                         path_to_foam = path_to_foam,
                         Time = float(time), 
                         Field = Field_Ux, 
                         Scale = Scale_Ux, 
                         isVectorField = isVectorField_Ux, 
                         Screens = [{'save_path': './images/Ux_s'+s_for_saving+'_t'+time_for_saving+'_rotor.png',
                                     'pixel_width': pixel_width,
                                     'bbox': ((xmin_rotor, xmax_rotor ),(ymin_rotor, ymax_rotor))}, 
                                    
                                    {'save_path': './images/Ux_s'+s_for_saving+'_t'+time_for_saving+'_blade.png',
                                     'pixel_width': pixel_width,
                                     'bbox': ((xmin_blade, xmax_blade ),(ymin_blade, ymax_blade))} ],
                          Ambient = Ambient,
                          Diffuse = Diffuse)
