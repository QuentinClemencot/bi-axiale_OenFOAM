#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
sys.path.append("/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axe_multi_blade/py/func")

from postprocess_function_def import *
from general_function_def import *
import os
import matplotlib.pyplot as plt
import numpy as np


'''
The purpose of this script is the visualisation of data contain in 
postProcessing/solverInfo of an oscilating airfoil simulation.

It generate a figure containing the evolution of residual through time

This scipt must be launch from an OpenFoam working directory.

It python version is 3.7.1



'''



##################### Define post-process Output directory ##############

path_postprocess_output='./py_postprocess_data'
try: 
 os.makedirs(path_postprocess_output)
except OSError:
 if not os.path.isdir(path_postprocess_output):
 	    raise ValueError('No dir ' + path_postprocess_output)

##################### Define Post-process Input directory ##############

path_postprocess_input='./postProcessing/solverInfo'

if not os.path.exists(path_postprocess_input):
    raise ValueError('No file ' + path_postprocess_input)

input_data_file_name='solverInfo.dat'


##################### List  time subdirectories in  Post-process Input directory ##############


time_list=read_time_dir(path_postprocess_input)


##################### Read  all time files and concatenate data ##############

data = read_OF_data(path_postprocess_input,time_list,OF_version='1912plus')




#######################  Figure Generation   ###################
fig, ax = plt.subplots(figsize=(10,10))


ax.plot(data['Time'],data['Ux_initial'],'b',label='Ux_initial')
ax.plot(data['Time'],data['Ux_final'],'b-',label='Ux_final')
ax.plot(data['Time'],data['Uy_initial'],'g',label='Uy_initial')
ax.plot(data['Time'],data['Uy_final'],'g-',label='Uy_final')
ax.plot(data['Time'],data['p_initial'],'r',label='p_initial')
ax.plot(data['Time'],data['p_final'],'r-',label='p_final')
ax.plot(data['Time'],data['k_initial'],'c',label='k_initial')
ax.plot(data['Time'],data['k_final'],'c-',label='k_final')
ax.plot(data['Time'],data['omega_initial'],'y',label='omega_initial')
ax.plot(data['Time'],data['omega_final'],'y-',label='omega_final')
ax.set(xlabel='t (s)', ylabel='residuals')
ax.grid()
ax.legend(loc='best')


plt.savefig(os.path.join(path_postprocess_output,'residuals.pdf'))

plt.yscale("log")
plt.show()
plt.close(fig)


























