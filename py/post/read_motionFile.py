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
The purpose of this script is the visualisation of constant/motion_blade*.dat files.

It use the files : 'constant/motion_blade*.dat' and './simulationProperties.py'


This scipt must be launch from an OpenFoam working directory.

It python version is 3.7.1



'''
##################### Read usefull parameters ##############
r=read_simulation_properties('r')
h=read_simulation_properties('h')
Nb=read_simulation_properties('Nb')
T=read_simulation_properties('T')


##################### Define folder containing motion files ##############

path_motion_folder='./constant/'

#'./constant/' contains files named motion_blade*.dat. 
#example : if Nb=3, there is motion_blade1.dat, motion_blade2.dat, motion_blade3.dat.

if not os.path.exists(path_motion_folder):
    raise ValueError('No file ' + path_motion_folder)


##################### Define post-process Output directory ##############

path_postprocess_output='./py_postprocess_data'
try: 
 os.makedirs(path_postprocess_output)
except OSError:
 if not os.path.isdir(path_postprocess_output):
 	    raise ValueError('No dir ' + path_postprocess_output)



#################### Create dictionnary containing motion data ######################
data={}

for i in range(Nb):

	key='blade'+str(i+1)
	data[key]={}

	#read motion file :
	data[key] = read_motion(path_motion_folder,'motion_blade'+str(i+1)+'.dat',OF_version='1912plus')
	
	#create TimeAdim=t/T key
	data[key]['TimeAdim']=data[key]['Time']/T

	#create s key
	size=len(data[key]['TimeAdim'])
	data[key]['s']=data[key]['TimeAdim']+np.ones(size)*i/Nb



##################### figure generation ##############

for key in data :
	#trajectory readen from motion_blade*.dat:
	fig, axs = plt.subplots(1,4,figsize=(15,4))
	fig.suptitle('Trajectory of the center of rotation of '+key+' for h='+str(h)+', r='+str(r))
	plt.subplots_adjust(wspace = 0.3)

	axs[0].plot(data[key]['x-x0'],data[key]['y-y0'],'.')
	axs[0].set(xlabel='x - x(t=0)', ylabel='y-y(t=0)')
	axs[0].grid()

	axs[1].plot(data[key]['s'],data[key]['theta-theta0'],'.')
	axs[1].set(xlabel='t/T', ylabel='theta-theta(t=0)')
	axs[1].grid()


	axs[2].plot(data[key]['s'],data[key]['y-y0'],'.')
	axs[2].set(xlabel='t/T', ylabel='y - y(t=0)')
	axs[2].grid()


	axs[3].plot(data[key]['s'],data[key]['x-x0'],'.')
	axs[3].set(xlabel='t/T', ylabel='x - x(t=0)')
	axs[3].grid()

	plt.show()
	plt.close(fig)


	#trajectory readen from motion.dat
	fig, axs = plt.subplots(3,3,figsize=(15,15))
	fig.suptitle('Position, velocity and acceleration of the CofR and theta angle of '+key+' for h='+str(h)+', r='+str(r))
	plt.subplots_adjust(wspace = 0.3)
	plt.subplots_adjust(hspace = 0.3)
	'''
	axs[0,0].plot(data[key]['x-x0'],data[key]['y-y0'])
	axs[0,0].set(xlabel='x - x(t=0)', ylabel='y-y(t=0)')
	axs[0,0].set_xlim([-1,1])
	axs[0,0].set_ylim([-1,1])
	axs[0,0].grid()
	'''
	axs[0,0].plot(data[key]['s'],data[key]['theta-theta0'])
	axs[0,0].set(xlabel='t/T', ylabel='theta-theta(t=0)')
	axs[0,0].grid()


	axs[0,1].plot(data[key]['s'],data[key]['y-y0'])
	axs[0,1].set(xlabel='t/T', ylabel='y - y(t=0)')
	axs[0,1].grid()


	axs[0,2].plot(data[key]['s'],data[key]['x-x0'])
	axs[0,2].set(xlabel='Time', ylabel='x - x(t=0)')
	axs[0,2].grid()

	axs[1,0].plot(data[key]['s'],data[key]['dtheta'])
	axs[1,0].set(xlabel='Time', ylabel='dtheta/dt')
	dMax=max(data[key]['dtheta'])
	dMin=min(data[key]['dtheta'])
	if dMax-dMin<0.001:
		axs[1,0].set_ylim([dMin-1,dMax+1])
	axs[1,0].grid()


	axs[1,1].plot(data[key]['s'],data[key]['dy'])
	axs[1,1].set(xlabel='Time', ylabel='dy/dt')
	axs[1,1].grid()


	axs[1,2].plot(data[key]['s'],data[key]['dx'])
	axs[1,2].set(xlabel='t/T', ylabel='dx/dt')
	axs[1,2].grid()

	axs[2,0].plot(data[key]['s'],data[key]['d2theta'])
	axs[2,0].set(xlabel='t/T', ylabel='d²theta/dt²')
	dMax=max(data[key]['d2theta'])
	dMin=min(data[key]['d2theta'])
	if dMax-dMin<0.001:
		axs[2,0].set_ylim([dMin-1,dMax+1])
	axs[2,0].grid()


	axs[2,1].plot(data[key]['s'],data[key]['d2y'])
	axs[2,1].set(xlabel='t/T', ylabel='d²y/dt²')
	axs[2,1].grid()


	axs[2,2].plot(data[key]['s'],data[key]['d2x'])
	axs[2,2].set(xlabel='t/T', ylabel='d²x/dt²')
	axs[2,2].grid()

	plt.show()
	plt.close(fig)







