#!/usr/bin/env python
# -*- coding: utf-8 -*-

# load home made function for case initialization:
import sys
sys.path.append("/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/\
bi_axial_turbine/py/func")

from preprocess_function_def import init_multiBlades_data, \
                                    add_CofR, \
                                    add_gamma, \
                                    add_theta, \
                                    add_settingAngle, \
                                    add_Time, \
                                    write_motionData, \
                                    write_dynamicMeshDict, \
                                    write_init_pos
                                    
from general_function_def import read_simulation_properties, \
                                 make_it_continuous
# load usefull libs:                                      
import os
import numpy as np
from math import pi

'''
The purpose of this script is the creation of the dynamicMeshDict and
blade tabulated position files

This scipt must be launch from an OpenFoam working directory.

It python version is 3.7.1



'''
# Read usefull parameters ##############

T = read_simulation_properties('T')
N = read_simulation_properties('nRevRegular') + \
    read_simulation_properties('nRevProbes')
Nstep = max(read_simulation_properties('nStep_per_revRegular'), \
            read_simulation_properties('nStep_per_revProbes'))
Nb = read_simulation_properties('Nb')
r = read_simulation_properties('r')
h = read_simulation_properties('h')
Uref = read_simulation_properties('Uref')
Lambda = read_simulation_properties('Lambda')
beta0 = read_simulation_properties('beta0')
lawType = read_simulation_properties('lawType')
xb = read_simulation_properties('xb')
xa = read_simulation_properties('xa')
c = read_simulation_properties('c')


# Define Output directory ##############
path_wd = os.getcwd()  # path working directory
path_output = './constant'

try:
    os.makedirs(path_output)
except OSError:
    if not os.path.isdir(path_output):
        raise ValueError('No dir ' + path_output)


# Define position parameter for every blade ##############

data = {}


# s is a dimensionless position parameter
# between two successive time step, s increases by ds
(s, ds) = np.linspace(0, 1, num=Nstep + 1, retstep=True)


data = init_multiBlades_data(s, Nb)


# add postion of the center of rotation #####

for key in data:
    data[key] = add_CofR(data[key], h, r)


# add angle between x direction and CofR trajectory's tangent ########
'''
-pi<gamma<=pi
'''

for key in data:
    data[key] = add_gamma(data[key], h, r)
    data[key]['gammaTot'] = make_it_continuous(data[key]['gamma'], -pi, pi)


# add angle between x direction and chord ########
'''
-pi<theta<=pi
'''
for key in data:
    data[key] = add_theta(beta0, c, xb, xa, data[key], h, lawType, r)
    data[key]['thetaTot'] = make_it_continuous(data[key]['theta'], -pi, pi)


# add setting angle beta=theta-gamma ###########
for key in data:
    data[key] = add_settingAngle(data[key])


# concatenate data N + 1 times and add 'Time' key in a new dict ###########
'''
tabulated data will contain position for N+1 revolution (N revolution will be
simulated). Therefore, you are sure that even for the last time step simulated,
 position are store in tabulated data
'''

for key in data:
    data[key] = add_Time(data[key], T, N + 1, Nstep)


# write motion data #######################

write_motionData(data, path_output)

# write dynamicMeshDict #######################

write_dynamicMeshDict(data, path_output, path_wd)

# write initial position of blade #######################

write_init_pos(data, './constant', 'blade_init_pos.dat')


