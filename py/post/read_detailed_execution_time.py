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
The purpose of this script is the visualisation of the repartition of execution 
time for simulation case solve with overPimpleDyMTimeRecordFoam.

overPimpleDyMTimeRecordFoam is an home made solver, source code : 
/fsnet/project/nrj/2019/19TRIBINE/simulation/OpenFOAM/clemenco1q-v2006/
application/overPimpleDyMTimeRecordFoam

It is a copy of overPimpleDyMFoam OF v2006 with extra execution time infomation.
It permits to known mesh update and pimple loop execution time.



This scipt must be launch from an OpenFoam working directory.

It python version is 3.7.1



'''

##################### Usefull library ##############

import numpy as np 



##################### Post-process Input directory ##############

wd='./'
log_file_name='log'


##################### Read usefull parameters ##############
isSimProperties=True

path_in_simuProperties = os.path.join(wd, 'simulationProperties.py')
if not os.path.exists(path_in_simuProperties):
	isSimProperties=False

if isSimProperties: T=read_simulation_properties('T')


#################### Create dictionnary containing data ######################

data = read_execution_time(wd,log_file_name,OF_version='2006plus')
	
sizeTime=len(data['Time'])
sizeMeshUp=len(data['MeshUp'])
sizePimpleLoop=len(data['PimpleLoop'])

meanMeshup=np.mean(data['MeshUp'])
meanPimpleLoop=np.mean(data['PimpleLoop'])

if isSimProperties: 
	data['TimeAdim']=data['Time']/T
	sizeTimeAdim=len(data['TimeAdim'])

print('Results of reading the log file : \n')
print('number of time step readen : '+str(sizeTime))

if isSimProperties: 
	print('from t='+str(round(data['Time'][0],5))+'s (t/T='+str(round(data['TimeAdim'][0],5))+')'+
		' to t='+str(round(data['Time'][-1],5))+'s (t/T='+str(round(data['TimeAdim'][-1],5))+')')
else : 
	print('from t='+str(round(data['Time'][0],5))+'s '+
		' to t='+str(round(data['Time'][-1],5))+'s')

print('number of mesh update time readen : '+str(sizeMeshUp))
print('number of pimple loop time readen : '+str(sizePimpleLoop))

print('Average mesh update time : '+str(round(meanMeshup,3)))
print('Average pimple loop time : '+str(round(meanPimpleLoop,3)))

#######################  Figure Generation   ###################
if isSimProperties:
	if sizeTime==sizeMeshUp+1 and sizeMeshUp==sizePimpleLoop : 
		data['TimeAdim']=np.delete(data['TimeAdim'], 0)
		sizeTimeAdim=sizeTime-1

	if sizeTime==sizePimpleLoop+1 and sizeMeshUp==sizeTime : 
		data['TimeAdim']=np.delete(data['TimeAdim'], 0)
		sizeTimeAdim=sizeTime-1
		data['MeshUp']=np.delete(data['MeshUp'], 0)
		sizeMeshUp=sizeMeshUp-1

	if sizeTimeAdim==sizeMeshUp and sizeMeshUp==sizePimpleLoop :
		#force coeff=f(t) 
		fig, axs = plt.subplots(figsize=(10,5))

		axs.plot(data['TimeAdim'],data['MeshUp'],label='Mesh update')
		axs.plot(data['TimeAdim'],data['PimpleLoop'],label='Pimple loop')
		axs.set(xlabel='t/T', ylabel='execution time (ms)')
		axs.grid()
		axs.legend(loc='best')
		plt.show()
		plt.close(fig)

else :
	if sizeTime==sizeMeshUp+1 and sizeMeshUp==sizePimpleLoop : 
		data['Time']=np.delete(data['Time'], 0)
		sizeTime=sizeTime-1

	if sizeTime==sizeMeshUp and sizeMeshUp==sizePimpleLoop :
		#force coeff=f(t) 
		fig, axs = plt.subplots(figsize=(10,5))

		axs.plot(data['Time'],data['MeshUp'],label='Mesh update')
		axs.plot(data['Time'],data['PimpleLoop'],label='Pimple loop')
		axs.set(xlabel='t simulated (s)', ylabel='execution time (ms)')
		axs.grid()
		axs.legend(loc='best')
		plt.show()
		plt.close(fig)






