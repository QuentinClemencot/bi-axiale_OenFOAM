#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
sys.path.append("/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axe_multi_blade/py/func")
sys.path.append("/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axe_multi_blade/py/post")

from preprocess_function_def import (gammaAngle, \
                                     x_coord, \
                                     y_coord, \
                                     thetaLawType3)

from postprocess_function_def import (read_time_dir, reorder_points) 

from general_function_def import (closest_rank_element, \
                                  read_data, \
                                  read_simulation_properties, \
                                  read_OF_scalar_data, \
                                  read_OF_vector_data, \
                                  rot_array)
import os
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.animation as animation
from math import pi
import copy
'''
The purpose of this script is the visualisation of power coefficient of 
multi blade bi-axial turbine.

It use the folder : './postProcessing'


This scipt must be launch from an OpenFoam working directory.

It python version is 3.8.2



'''
##################### Read usefull parameters ##############
Nb = read_simulation_properties('Nb') #number of blade
T = read_simulation_properties('T')
r = read_simulation_properties('r')
h = read_simulation_properties('h')
c = read_simulation_properties('c')
xA = read_simulation_properties('xl')
xB = read_simulation_properties('xc')
###useless:

N = read_simulation_properties('N')

Nstep = read_simulation_properties('Nstep') #number of time step per revoluiton


b = read_simulation_properties('b') #blade span [m]
Uref = read_simulation_properties('Uref')
Lambda = read_simulation_properties('Lambda')
beta0 = read_simulation_properties('beta0')
lawType = read_simulation_properties('lawType')

rho = read_simulation_properties('rho')

pdyn = 0.5 * rho * (Uref**2) #dynamique pressure
d = 2*r + h #height of the area swept by the blades



##################### determine fields that will be print ##############

printCp = False
##################### Post-process Input directory ##############

path_surfaces = './postProcessing/surfaces'

if not os.path.exists(path_surfaces):
    raise ValueError('No file ' + path_surfaces)

p_file_name = 'scalarField/p'
wss_file_name = 'vectorField/wallShearStress'
faceCenters_file_name = 'faceCentres'


if printCp:
    path_powerCoeff = './py_postprocess_data' 
    if not os.path.exists(path_powerCoeff):
        raise ValueError('No file ' + path_powerCoeff)
    powerCoeff_file_name = 'CpAllTime.txt'
##################### Define post-process Output directory ##############

path_postprocess_output='./py_postprocess_data'
try: 
    os.makedirs(path_postprocess_output)
except OSError:
    if not os.path.isdir(path_postprocess_output):
        raise ValueError('No dir ' + path_postprocess_output)


#################### Create dictionnary containing data ######################

data = {f"blade{i}":{} for i in range(1, Nb+1)}
#data is a dictionnary containing sub-dict keys "blade1", "blade2", ect..  


time_list = read_time_dir(path_surfaces)
#time_list is an ascending order list of time directory (str type) stored in path_surfaces


for time in time_list:
    for blade in data:
        path = os.path.join(path_surfaces, time, blade)
        print(f"is going to read path {path}")
        if not os.path.exists(os.path.join(path, p_file_name)):
            raise ValueError('No file ' + os.path.join(path, p_file_name) )
        (lenght_p, p) = read_OF_scalar_data(path, p_file_name)
        
        if not os.path.exists(os.path.join(path, wss_file_name)):
            raise ValueError('No file ' + os.path.join(path, wss_file_name) )
        (lenght_wss, wss_xo, wss_yo, useless) = read_OF_vector_data(path, wss_file_name)

        if not os.path.exists(os.path.join(path, faceCenters_file_name)):
            raise ValueError('No file ' + os.path.join(path, faceCenters_file_name) )  
        (lenght_faceCenter, x, y, useless) = read_OF_vector_data(path, faceCenters_file_name)
		

        data[blade][time] = {}
        data[blade][time]['p'] = p
        data[blade][time]['wss'] = {'x0': wss_xo, 'y0': wss_yo}
        data[blade][time]['faceCenter'] = {'x0': np.array(x), 'y0': np.array(y)}

print("data containing in ./postProcessing/surfaces has been read.")


##################### add information about blade position ##############
for time in time_list:
    for blade in data:
        # First: determine position parameter of blades for each sample: 
        rank = int(blade.replace('blade','')) #rank == 1 for blade1, ect..
        s0 = (rank - 1) / Nb                  #position parameter at t=0
        s = ((float(time) + s0 * T) % T) / T #lets s3(t) the posiion of the rank 3 blade.
						                     #s3(t=0)=s0=s1(t=s0*T)
        data[blade][time]['s'] = s

        # Second: determine position of A point:
        data[blade][time]['pointA'] = {'x0': x_coord(s, h, r), 'y0': y_coord(s, h, r)}
        
        # Third: determine theta angle between xo direction and the chord of the blade:
        # theta is given in radian. 
        tmp_values= thetaLawType3(c, \
                                  [data[blade][time]['pointA']['x0']], \
                                  [data[blade][time]['pointA']['y0']], \
                                  [gammaAngle(s, h, r)], \
                                  h, \
                                  r, \
                                  [s], \
                                  1, \
                                  xA, \
                                  xB)

        data[blade][time]['gamma'] = tmp_values[0][0]
        data[blade][time]['pointB'] = {'x0': tmp_values[1][0], 'y0': tmp_values[2][0]}            
##################### reorder blade profil point in a correct order ##############
for time in time_list:
    for blade in data:
        data[blade][time] = reorder_points(data[blade][time])
            

#data['blade1'][time_list[0]] = reorder_points(data['blade1'][time_list[0]])


##################### read power coefficient #################
if printCp:
    data_temporal = read_data(path_powerCoeff, powerCoeff_file_name)

    #extract last period
    data_temporal['s'] = [t % 1 for t in data_temporal['TimeAdim']]

    tAdimMax = data_temporal['TimeAdim'][-1]
    extractMax = int(tAdimMax)
    extractMin = extractMax - 1

    rankMax = closest_rank_element(data_temporal['TimeAdim'], extractMax)
    rankMin = closest_rank_element(data_temporal['TimeAdim'], extractMin)

    data_temporal_lastPeriod = copy.deepcopy(data_temporal)

    for key in data_temporal_lastPeriod:
        data_temporal_lastPeriod[key] = data_temporal_lastPeriod[key][rankMin:rankMax]


##################### Point A position during a revolution ##############
trajectory_s = np.linspace(0, 1, 1000)
trajectory_x = np.array([x_coord(s, h, r) for s in trajectory_s])
trajectory_y = np.array([y_coord(s, h, r) for s in trajectory_s])



######################Create animation ########################

fig = plt.figure(figsize=(10, 10))
ax1 = fig.add_subplot(2, 4, 1) #up left
ax1bis = fig.add_subplot(2, 4, 5) #down left
ax2 = fig.add_subplot(2, 4, 2) #up left
ax3 = fig.add_subplot(2, 4, 6) #down left
ax4 = fig.add_subplot(2, 4, 3) #down left
ax5 = fig.add_subplot(2, 4, 7) #down left
ax6 = fig.add_subplot(2, 4, 4) #down left
ax7 = fig.add_subplot(2, 4, 8) #down left

ax1.plot(trajectory_x, trajectory_y, color='k')
ax1.set_xlabel("x [m]")
ax1.set_ylabel("y [m]")
ax1.set_aspect('equal', 'box')
blade1_Ax_t0 = data['blade1'][time_list[0]]['pointA']['x0']
blade1_Ay_t0 = data['blade1'][time_list[0]]['pointA']['y0']

if Nb >= 2:
    blade2_Ax_t0 = data['blade2'][time_list[0]]['pointA']['x0']
    blade2_Ay_t0 = data['blade2'][time_list[0]]['pointA']['y0']
    scat2 = ax1.scatter(blade2_Ax_t0, blade2_Ay_t0, s=75., color='r')
scat1 = ax1.scatter(blade1_Ax_t0, blade1_Ay_t0, s=75., color='b')


#power coefficient :
if printCp:
    ax1bis.plot(data_temporal_lastPeriod['s'], data_temporal_lastPeriod['CpTotal'], color='k')
    ax1bis.plot(data_temporal_lastPeriod['s'], data_temporal_lastPeriod['blade1'], color='b')
    if Nb >= 2: ax1bis.plot(data_temporal_lastPeriod['s'], data_temporal_lastPeriod['blade2'], color='r')
    ax1bis.set_xlabel("t/T")
    ax1bis.set_ylabel("power coefficient")
    ax1bis.grid()

    i0 = closest_rank_element(data_temporal_lastPeriod['s'], (float(time_list[0]) % T) / T)
    scat_Cp_blade1 = ax1bis.scatter(data_temporal_lastPeriod['s'][i0], data_temporal_lastPeriod['blade1'][i0], s=75., color='b')
    if Nb >= 2: scat_Cp_blade2 = ax1bis.scatter(data_temporal_lastPeriod['s'][i0], data_temporal_lastPeriod['blade2'][i0], s=75., color='r')
    scat_Cp_total = ax1bis.scatter(data_temporal_lastPeriod['s'][i0], data_temporal_lastPeriod['CpTotal'][i0], s=75., color='k')

#blades orientations :
init_blade1_x = data['blade1'][time_list[0]]['faceCenter']['x0']
init_blade1_y = data['blade1'][time_list[0]]['faceCenter']['y0']
blade1_Ax2_t0 = data['blade1'][time_list[0]]['pointA']['x0']
blade1_Ay2_t0 = data['blade1'][time_list[0]]['pointA']['y0']
blade1_Bx2_t0 = data['blade1'][time_list[0]]['pointB']['x0']
blade1_By2_t0 = data['blade1'][time_list[0]]['pointB']['y0']

ax2.set_aspect('equal')
ax2.set_xlim(-c, c)
ax2.set_ylim(-c, c)
lineBlade1Traj, = ax2.plot(trajectory_x, trajectory_y, color='k')
lineBlade1, = ax2.plot(init_blade1_x,init_blade1_y)
scat1_ax2 = ax2.scatter(blade1_Ax2_t0, blade1_Ay2_t0, s=75., color='b')
scat2_ax2 = ax2.scatter(blade1_Bx2_t0, blade1_By2_t0, s=75., color='g')
if Nb >= 2:
    init_blade2_x = data['blade2'][time_list[0]]['faceCenter']['x0']
    init_blade2_y = data['blade2'][time_list[0]]['faceCenter']['y0']
    ax3.set_aspect('equal')
    ax3.set_xlim(-c, c)
    ax3.set_ylim(-c, c)
    lineBlade2, = ax3.plot(init_blade2_x,init_blade2_y)

#pressure coefficient : 
init_blade1_xadim = data['blade1'][time_list[0]]['x_adim']
init_blade1_Cp = data['blade1'][time_list[0]]['p']
ax4.set_xlabel("relative position on chord")
ax4.set_ylabel("Pressure Coeff")
lineBlade4, = ax4.plot(init_blade1_xadim,init_blade1_Cp)
ax4.grid()

if Nb >= 2:
    init_blade2_xadim = data['blade2'][time_list[0]]['x_adim']
    init_blade2_Cp = data['blade2'][time_list[0]]['p']
    ax5.set_xlabel("relative position on chord")
    ax5.set_ylabel("Pressure Coeff")
    lineBlade5, = ax5.plot(init_blade2_xadim,init_blade2_Cp)
    ax5.grid()

for ax, color in zip([ax2, ax3, ax4, ax5], ['blue', 'red', 'blue', 'red']):
    plt.setp(ax.spines.values(), color=color)
    plt.setp([ax.get_xticklines(), ax.get_yticklines()], color=color)

for ax in (ax4, ax5):
    ax.set_ylim(-30, 25)

def animate(i, data, data_temporal_lastPeriod = None):
    #ax1 update:
    tmp_scat1 = [data['blade1'][time_list[i]]['pointA']['x0'], \
                 data['blade1'][time_list[i]]['pointA']['y0']]
    scat1.set_offsets(tmp_scat1)
    if Nb >= 2:
        tmp_scat2 = [data['blade2'][time_list[i]]['pointA']['x0'], \
                     data['blade2'][time_list[i]]['pointA']['y0']]
        scat2.set_offsets(tmp_scat2)

    
   
    #ax1bis update:
    if printCp:
        i0 = closest_rank_element(data_temporal_lastPeriod['s'], (float(time_list[i]) % T) / T)
        tmp_scat_Cp_blade1 = [data_temporal_lastPeriod['s'][i0], \
                              data_temporal_lastPeriod['blade1'][i0]]

        tmp_scat_CpTotal = [data_temporal_lastPeriod['s'][i0], \
                          data_temporal_lastPeriod['CpTotal'][i0]]
        scat_Cp_blade1.set_offsets(tmp_scat_Cp_blade1)

        scat_Cp_total.set_offsets(tmp_scat_CpTotal)   
        if Nb >= 2:
            tmp_scat_Cp_blade2 = [data_temporal_lastPeriod['s'][i0], \
                                 data_temporal_lastPeriod['blade2'][i0]]
            scat_Cp_blade2.set_offsets(tmp_scat_Cp_blade2)

 
    #ax2 update:    
    blade1_x = data['blade1'][time_list[i]]['faceCenter']['x0'] - data['blade1'][time_list[i]]['pointA']['x0']
    blade1_y = data['blade1'][time_list[i]]['faceCenter']['y0'] - data['blade1'][time_list[i]]['pointA']['y0']
    blade1B_x = data['blade1'][time_list[i]]['pointB']['x0'] - data['blade1'][time_list[i]]['pointA']['x0']
    blade1B_y = data['blade1'][time_list[i]]['pointB']['y0'] - data['blade1'][time_list[i]]['pointA']['y0']
    print(f"time = {time_list[i]}")
    print(f"T = {T}")
    print(f"s = {data['blade1'][time_list[i]]['s']}")
    print(f"max yo={max(data['blade1'][time_list[i]]['faceCenter']['y0'])}")
    print(f"min yo={min(data['blade1'][time_list[i]]['faceCenter']['y0'])}")
    print(f"yA = {data['blade1'][time_list[i]]['pointA']['y0']}")
    lineBlade1Traj.set_data(trajectory_x - data['blade1'][time_list[i]]['pointA']['x0'], trajectory_y - data['blade1'][time_list[i]]['pointA']['y0'])
    lineBlade1.set_data(blade1_x, blade1_y)
    scat1_ax2.set_offsets([0, 0])
    scat2_ax2.set_offsets([blade1B_x, blade1B_y])
    
    #ax3 update:  
    if Nb >= 2:     
        blade2_x = data['blade2'][time_list[i]]['faceCenter']['x0'] - data['blade2'][time_list[i]]['pointA']['x0']
        blade2_y = data['blade2'][time_list[i]]['faceCenter']['y0'] - data['blade2'][time_list[i]]['pointA']['y0']
        lineBlade2.set_data(blade2_x, blade2_y)

    #ax4 update:    
    init_blade1_xadim = data['blade1'][time_list[i]]['x_adim']
    init_blade1_Cp = data['blade1'][time_list[i]]['p']
    lineBlade4.set_data(init_blade1_xadim, init_blade1_Cp)

    #ax5 update:
    if Nb >= 2:  
        init_blade2_xadim = data['blade2'][time_list[i]]['x_adim']
        init_blade2_Cp = data['blade2'][time_list[i]]['p']
        lineBlade5.set_data(init_blade2_xadim, init_blade2_Cp)
    
    if Nb >= 2:  return  scat1, scat2, lineBlade1, lineBlade2, lineBlade4, lineBlade5
    else : return  scat1,  lineBlade1, lineBlade4


if printCp : 
    ani = animation.FuncAnimation(
        fig, animate, frames=range(0, len(time_list)), interval=400, save_count=50, fargs = (data, data_temporal_lastPeriod))
else:
    ani = animation.FuncAnimation(
        fig, animate, frames=range(0, len(time_list)), interval=400, save_count=50, fargs = (data, None))
# To save the animation, use e.g.
#
# ani.save("movie.mp4")
#
# or
#
# from matplotlib.animation import FFMpegWriter
# writer = FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# ani.save("movie.mp4", writer=writer)

plt.show()



