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
N=read_simulation_properties('N')
Nb=read_simulation_properties('Nb')
Nstep=read_simulation_properties('Nstep')
r=read_simulation_properties('r')
h=read_simulation_properties('h')
c=read_simulation_properties('c')
Uref=read_simulation_properties('Uref')
Lambda=read_simulation_properties('Lambda')
beta0=read_simulation_properties('beta0')
lawType=read_simulation_properties('lawType')
xl=read_simulation_properties('xl')
xc=read_simulation_properties('xc')

##################### Post-process Input directory ##############

path_postProcessing='./postProcessing/'

#'./postProcessing/' contains folders named forceCoeffsBlade*. 
#example : if Nb=3, there is forceCoeffsBlade1, forceCoeffsBlade2, forceCoeffsBlade3.

if not os.path.exists(path_postProcessing):
    raise ValueError('No file ' + path_postProcessing)

input_data_file_name='coefficient.dat'


##################### Define post-process Output directory ##############

path_postprocess_output='./py_postprocess_data'
try: 
 os.makedirs(path_postprocess_output)
except OSError:
 if not os.path.isdir(path_postprocess_output):
 	    raise ValueError('No dir ' + path_postprocess_output)


#################### Create dictionnary containing data ######################
data={}

for i in range(Nb):
	blade='blade'+str(i+1)
	data[blade]={}


	path_postprocess_input=path_postProcessing+'forceCoeffsBlade'+str(i+1)
	print(path_postprocess_input)
	if not os.path.exists(path_postprocess_input):
		raise ValueError('No file ' + path_postprocess_input)

	#List  time subdirectories in  Post-process Input directory :
	time_list=read_time_dir(path_postprocess_input)

	#Read  all time files and concatenate data :
	data[blade] = read_OF_data(path_postprocess_input,time_list,OF_version='1912plus')

	#force coeff readen in postprocessing do not represent anymore lift and drag forces
	#because the CofR is moving throught time. they are rename as : 
	if 'Cd' in data[blade]  : data[blade] ['Cx'] = data[blade] .pop('Cd')
	else :
		raise ValueError('No key Cd in data file ' + path_postprocess_input)

	if 'Cl' in data[blade] : data[blade] ['Cy'] = data[blade] .pop('Cl')
	else:
		raise ValueError('No key Cl in data file ' + path_postprocess_input)

	'''
	CmPitch give the pitching moment calculated at reference frame origine point O(0,0)
	It value change if we calculate it at an other point (following the "BABAR" rule).

	So we define a key 'CmO' which is a subDict containing 3 keys :
	 'x':float=0,'y':float=0, 'value':float

	We will define later the key 'CmA' : the pitching moment coeff at the blade rotation point A,
	deduct from CmO with the BABAR rules
	'''
	data[blade]['CmO']={}
	if 'CmPitch' in data[blade] : data[blade]['CmO']['value'] = data[blade].pop('CmPitch')
	size=len(data[blade]['CmO']['value'])
	print('size='+str(size))
	data[blade]['CmO']['x']=np.zeros(size)
	data[blade]['CmO']['y']=np.zeros(size)

	#create TimeAdim=t/T key
	data[blade]['TimeAdim']=data[blade]['Time']/T

	print(str(len(data[blade]['Time']))+' time steps have been read from t='+str(data[blade]['Time'][0])+
	's to t='+str(data[blade]['Time'][-1]))





##################### add to data the s position parameter ##############
for blade in data :
	rank=blade.replace('blade','')
	rank=int(rank)
	print('read information of '+blade+' , the blade at position '+str(rank))
	data[blade] = add_posParam(data[blade],h,r,T,Nb,rank)



print('max s='+str(max(data[blade]['s']))+', min s='+str(min(data[blade]['s'])))


################## add postion of the center of rotation #####
for blade in data :
	data[blade]=add_CofR(data[blade],h,r)
	print('max x='+str(max(data[blade]['CofR']['x']))+', min x='+str(min(data[blade]['CofR']['x'])))
	print('max y='+str(max(data[blade]['CofR']['y']))+', min y='+str(min(data[blade]['CofR']['y'])))







################## add relative velocity of the fluid with respect to the blade ########

#data=add_Urel(data,h,r,Uref,Lambda)



################## add angle between x direction and CofR trajectory's tangent ########
'''
-pi<gamma<=pi
'''

for blade in data : 
	data[blade] = add_gamma(data[blade],h,r)

	data[blade]['gammaTot']=make_it_continuous(data[blade]['gamma'],-pi,pi)


################## add angle between x direction and chord ########
'''
-pi<theta<=pi
'''
for blade in data : 
	data[blade] = add_theta(beta0,c,xl,xc,data[blade],h,lawType,r)

	data[blade]['thetaTot']=make_it_continuous(data[blade]['theta'],-pi,pi)

################## add setting angle beta=theta-gamma ###########

for blade in data : 
	data[blade]=add_settingAngle(data[blade])


################## add rotation speed ###########

for blade in data : 
	data[blade]=add_omega(data[blade]) #omega=angular velovity=d/dt(totRot) [rad/s]



######### add CmA=pitching moment coefficient calculated at blade rotation point ########
for blade in data : 
	data[blade]=add_CmA(data[blade],c)


################## add tangential and normal force coefficient ########

'''
let xo be the (1,0,0) vector in the fixed referentiel 
	yo be the (0,1,0)
let t(s) be  vector tangential to the trajectory at position s
let n(s) be  vector normal to the trajectory at position s, with :

t(s=0)=-yo
n(s=0)=xo
t(s=0.5)=yo
n(s=0.5)=-xo

Cx is the force coefficient in the x direction. 
Cy is the force coefficient in the y direction.

Ct is the force coefficient in the tangential to the trajectory direction.
Cn is the force coefficient in the normal to the trajectory direction.

'''
for blade in data : 
	data[blade]=add_Ct_Cn(data[blade])


################### add power coefficient ############################

for blade in data :
	print('len(data[blade][CmA][value])=')
	print(len(data[blade]['CmA']['value']))

	print('len(data[Ct])=')
	print(len(data[blade]['Ct']))

	print('len(data[omega])=')
	print(len(data[blade]['omega']))

	data[blade]=add_Cp_Cpt_Cpr(data[blade],c,h,r,Uref,Lambda)

	'''other way to calculate power coeff : CmO*omega
	pertinent only for Darrius (with fixed setting angle?)
	'''
	data[blade]=add_CpO(data[blade],c,h,r,Uref,Lambda)




################## return Cp every nT time #########################
'''
extract Cp mean on every cycle for every blade to see convergence cycle to cycle
'''

Cp_evol={}
if data['blade1']['Time'][-1]>2*T:
	for blade in data :
		Cp_evol[blade]={}
		Cp_evol[blade]['value']=extract_Cp_evolution(data[blade],T)
		#['value'] key contains the average Cp of each cycle
		#['variation'] key contains pourcentage of change of Cp average from one
		# cycle to another
		size=len(Cp_evol['blade1']['value'])
		Cp_evol[blade]['variation']=np.zeros(size)
		for i in range(size-1):
			var=Cp_evol[blade]['value'][i+1]-Cp_evol[blade]['value'][i]
			var=100*var/Cp_evol[blade]['value'][i]
			Cp_evol[blade]['variation'][i+1]=var

	size=len(Cp_evol['blade1']['value'])
	Cp_evol_tot={}
	Cp_evol_tot['value']=np.zeros(size)
	Cp_evol_tot['variation']=np.zeros(size)
	for blade in data :
		Cp_evol_tot['value']=Cp_evol_tot['value']+Cp_evol[blade]['value']
	for i in range(size-1):
		var=Cp_evol_tot['value'][i+1]-Cp_evol_tot['value'][i]
		var=100*var/Cp_evol_tot['value'][i]
		Cp_evol_tot['variation'][i+1]=var
	
	Cp_evol['total']={}
	Cp_evol['total']=copy.deepcopy(Cp_evol_tot)




##################### extract only last period from  data ##############
'''
we check if the simulation has runned for more than one period :
'''
data_lastPer={}

if data['blade1']['Time'][-1]>T:
	meanGlobalCp=0
	for blade in data :
		data_lastPer[blade]=extract_period(data[blade],T,position='last')

		#we reorder data_lastPer by s ascending order:
		#data_lastPer[blade]=buble_sort_dict(data_lastPer[blade],'s')

		#time step being constant, the average Cp is simply the mean of Cp(t) 
		#value during one periode
		data_lastPer[blade]['meanCP']=np.mean(data_lastPer[blade]['Cp'])
		meanGlobalCp+=data_lastPer[blade]['meanCP']

	
	#Here we create data_lastPer['global']['Cp'] which is the sum of Cp of every
	#blade with repect to time

	size=len(data_lastPer['blade1']['s'])
	Cp_tot=np.zeros(size)
	for blade in data :
		Cp_tot=Cp_tot+data_lastPer[blade]['Cp']

	data_lastPer['global']={}
	data_lastPer['global']['TimeAdim']=copy.deepcopy(data_lastPer['blade1']['TimeAdim'])
	data_lastPer['global']['s']=copy.deepcopy(data_lastPer['blade1']['s'])
	data_lastPer['global']['Cp']=copy.deepcopy(Cp_tot)	




#####################  extract only penultimate period from  data ##############

'''
we check if the simulation has runned for more than one period :
'''
data_penultimPer={}

if data['blade1']['Time'][-1]>2*T:
	meanPenultimCp=0
	for blade in data :
		data_penultimPer[blade]=extract_period(data[blade],T,position='penultimate')

		#we reorder data_lastPer by s ascending order:
		#data_lastPer[blade]=buble_sort_dict(data_lastPer[blade],'s')

		#time step being constant, the average Cp is simply the mean of Cp(t) 
		#value during one periode
		data_penultimPer[blade]['meanCP']=np.mean(data_penultimPer[blade]['Cp'])
		meanPenultimCp+=data_lastPer[blade]['meanCP']

	meanPenultimCp=meanPenultimCp
	
	#Here we create data_penultimPer['global']['Cp'] which is the sum of Cp of 
	#every blade with repect to time

	size=len(data_penultimPer['blade1']['s'])
	Cp_tot=np.zeros(size)
	for i in range(size) :
		for blade in data :
			Cp_tot[i]+=data_penultimPer[blade]['Cp'][i]

	data_penultimPer['global']={}
	data_penultimPer['global']['s']=copy.deepcopy(data_penultimPer['blade1']['s'])
	data_penultimPer['global']['TimeAdim']=copy.deepcopy(data_penultimPer['blade1']['TimeAdim'])
	data_penultimPer['global']['Cp']=copy.deepcopy(Cp_tot)	








##################### figure generation ##############




#plot Cp
for blade in data :
	#force coeff=f(t) 
	fig, axs = plt.subplots(figsize=(10,5))

	axs.plot(data[blade]['TimeAdim'],data[blade]['Cp'],'.',label='Cp')
	axs.plot(data[blade]['TimeAdim'],data[blade]['Cpt'],label='Cpt')
	axs.plot(data[blade]['TimeAdim'],data[blade]['Cpr'],label='Cpr')
	if h==0:
		axs.plot(data[blade]['TimeAdim'],data[blade]['CpO'],label='CpO')
	#axs.plot(data[blade]['TimeAdim'],data[blade]['omega'],label='omega')
	#axs.plot(data[blade]['TimeAdim'],data[blade]['CmA']['value'],label='CmA')
	#axs.plot(data[blade]['TimeAdim'],data[blade]['CmO']['value'],label='CmO')
	axs.set(xlabel='t/T', ylabel='Cp')
	#axs.set_ylim(-200,200)

	axs.grid()
	axs.legend(loc='best')
	plt.show()
	plt.close(fig)


#plot Cp variation cycle to cycle
if data['blade1']['Time'][-1]>2*T:
	fig, axs = plt.subplots(1,2,figsize=(10,5))
	for blade in Cp_evol :
		#force coeff=f(t) 
		axs[0].plot(Cp_evol[blade]['value'],label=blade)
		axs[1].plot(Cp_evol[blade]['variation'],label=blade)

	axs[0].set(xlabel='cycle number', ylabel='cycle average Cp')
	axs[0].grid()
	axs[0].legend(loc='best')

	axs[1].set(xlabel='cycle number', ylabel='variation of Cp cycle-to-cycle')
	axs[1].grid()
	axs[1].legend(loc='best')

	plt.savefig(os.path.join(path_postprocess_output,'Cp_variation.pdf'))
	plt.show()
	plt.close(fig)


#Comparison of two last period

if 'global' in data_penultimPer :

	#force coeff=f(t) 
	fig, axs = plt.subplots(figsize=(10,5))

	for key in data_penultimPer:
		if not key in data_lastPer :
			raise ValueError('key ' +key+' present in meanPenultimCp'+
                  ' but not in data_lastPer')

		axs.plot(data_lastPer[key]['s'],data_lastPer[key]['Cp'],label=key+' last')
		axs.plot(data_penultimPer[key]['s'],data_penultimPer[key]['Cp'],label=key+' penultimate')


	axs.set(xlabel='t/T', ylabel='Cp')
	axs.set_ylim(-50,50)

	axs.grid()
	axs.legend(loc='best')
	plt.show()
	plt.close(fig)



#plot Cp_tot on last period

if 'global' in data_lastPer :
	#force coeff=f(t) 
	fig, axs = plt.subplots(figsize=(10,5))

	for key in data_lastPer:
		if key=='global':
			axs.plot(data_lastPer[key]['TimeAdim'],data_lastPer[key]['Cp'],
				label='global (Cp='+str(round(meanGlobalCp,7))+')')
		else :
			axs.plot(data_lastPer[key]['TimeAdim'],data_lastPer[key]['Cp'],label=key)

	axs.set(xlabel='t/T', ylabel='Cp')
	axs.set_ylim(-50,50)

	axs.grid()
	axs.legend(loc='best')
	plt.savefig(os.path.join(path_postprocess_output,'Cp_lastPeriod.pdf'))
	plt.show()

	plt.close(fig)



advanceComparison=True

if  advanceComparison:
	translationOnePer=[] #1D array will contain starting and ending adim 
			# time of pure translation (ie A and B point and translation part)
			#the first element is starting time of translation, the second is the ending time
			#the third is the next starting time, ect...
			#ie even indice are starting time, odd indice are ending time.
	rotationOnePer=[] #same but for pure translation

	l = 2*pi*r + 2*h
	cbis = c *(xl - xc)
	arc=2*asin(0.5*cbis/r)*r #arc lenght between A and B when both point are 
								#on turn part

	trans1_start=cbis/l
	trans1_end=h/l

	rot1_start=(h+arc)/l
	rot1_end=0.5

	trans2_start=0.5+cbis/l
	trans2_end=0.5+h/l

	rot2_start=0.5+(h+arc)/l
	rot2_end=1

	translationOnePer.append(trans1_start)
	translationOnePer.append(trans1_end)
	translationOnePer.append(trans2_start)
	translationOnePer.append(trans2_end)

	rotationOnePer.append(rot1_start)
	rotationOnePer.append(rot1_end)
	rotationOnePer.append(rot2_start)
	rotationOnePer.append(rot2_end)


#plot CmA on last period
for blade in data_lastPer :
	#force coeff=f(t) 
	if blade != 'global':
		fig, axs = plt.subplots(figsize=(10,5))

		axs.plot(data_lastPer[blade]['s'][1:-1],data_lastPer[blade]['CmA']['value'][1:-1])
		axs.set(xlabel='t/T', ylabel='CmA')
		#axs.set_ylim(-200,200)

		#Add band of color for pure translation and pure rotation part
		axs.axvspan(rotationOnePer[0], rotationOnePer[1], alpha=0.15, color='blue')
		axs.axvspan(rotationOnePer[2], rotationOnePer[3], alpha=0.15, color='blue')
		axs.axvspan(translationOnePer[0], translationOnePer[1], alpha=0.15, color='red')
		axs.axvspan(translationOnePer[2], translationOnePer[3], alpha=0.15, color='red')

		axs.grid()
		axs.legend(loc='best')
		plt.savefig(os.path.join(path_postprocess_output,f"CmA_c{str(c).replace('.','p')}_Nb{Nb}_lambda{Lambda}.pdf"))
		plt.show()
		plt.close(fig)

#plot force coeff
for blade in data_lastPer :
	#force coeff=f(t) 
	if blade != 'global':
		fig, axs = plt.subplots(figsize=(10,5))

		axs.plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['Cx'],label='Cx')
		axs.plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['Cy'],label='Cy')
		axs.plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['Cn'],label='Cn')
		axs.plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['Ct'],label='Ct')
		axs.plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['Cpt'],label='Cpt')
		axs.plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['Cp'],label='Cp')
		axs.set(xlabel='t/T', ylabel='force coeff')
		#axs.set_ylim(-200,200)



		axs.grid()
		axs.legend(loc='best')
		plt.show()
		plt.close(fig)

##################### write data ##############
wdata={}

wdata['Time']=copy.deepcopy(data['blade1']['Time'])
wdata['TimeAdim']=copy.deepcopy(data['blade1']['TimeAdim'])

size=len(wdata['Time'])
wdata['CpTotal']=np.zeros(size)

for blade in data:
	if len(data[blade]['Cp'])!=size:
		raise ValueError('len(Cp) from key ' +blade+'='+
			str(len(data[blade]['Cp']))+' not egal to len(Time)='+str(size))
	wdata[blade]=copy.deepcopy(data[blade]['Cp'])
	wdata['CpTotal']=wdata['CpTotal']+data[blade]['Cp']


write_data(wdata,path_postprocess_output,'CpAllTime.txt',8)


wdata={}

for blade in Cp_evol:
	wdata[blade+'_value']=Cp_evol[blade]['value']
	wdata[blade+'_variation']=Cp_evol[blade]['variation']

write_data(wdata,path_postprocess_output,'CpEvol_CycleToCycle.txt',8)














