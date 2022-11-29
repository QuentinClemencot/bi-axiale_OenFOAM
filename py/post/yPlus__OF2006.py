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
The purpose of this script is the visualisation of power coefficient of 
multi blade bi-axial turbine.

It use the folder : './postProcessing'


This scipt must be launch from an OpenFoam working directory.

It python version is 3.7.1



'''
##################### Read usefull parameters ##############

T=read_simulation_properties('T')


##################### Post-process Input directory ##############


path_postprocess_input='./postProcessing/yPlusExtremum'
#print(path_postprocess_input)
if not os.path.exists(path_postprocess_input):
	raise ValueError('No file ' + path_postprocess_input)


##################### read yPlus min max mean information   ##############

#List  time subdirectories in  Post-process Input directory :
time_list=read_time_dir(path_postprocess_input)


datayPlus=read_OF_data(path_postprocess_input,time_list,OF_version='1912plus')

yPlusMax=max(datayPlus['max'])
yPlusMin=min(datayPlus['min'])
print('\nyPlus max='+str(yPlusMax)+'\tyPlusMin='+str(yPlusMin)+'\n')


datayPlus['TimeAdim']=datayPlus['Time']/T




##################### figure generation ##############


#plot yPlus=f(alpha)
fig, axs = plt.subplots(figsize=(10,5))
fig.suptitle('yPlus for blade 1')

if (len(datayPlus['TimeAdim'])>10):
	axs.plot(datayPlus['TimeAdim'][10:],datayPlus['min'][10:],'.',color='green',label='min')
	axs.plot(datayPlus['TimeAdim'][10:],datayPlus['max'][10:],'.',color='red',label='max')
	axs.plot(datayPlus['TimeAdim'][10:],datayPlus['average'][10:],'.',color='blue',label='average')
	axs.set(xlabel='t/T', ylabel='y+')
	axs.grid()
	axs.legend(loc='best')

else :
	axs.plot(datayPlus['TimeAdim'],datayPlus['min'],'.',color='green',label='min')
	axs.plot(datayPlus['TimeAdim'],datayPlus['max'],'.',color='red',label='max')
	axs.plot(datayPlus['TimeAdim'],datayPlus['average'],'.',color='blue',label='average')
	axs.set(xlabel='t/T', ylabel='y+')
	axs.grid()
	axs.legend(loc='best')



plt.show()
plt.close(fig)



