#!/usr/bin/env python
# -*- coding: utf-8 -*-


import sys
sys.path.append("/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axe_multi_blade/py/func")
sys.path.append("/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axe_multi_blade/py/post")

from postprocess_function_def import *
from general_function_def import *
from blade_dynamique_carac import xG_position_OO12, area_moment_inertia_OO12, section_OO12
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

T = read_simulation_properties('T')
N = read_simulation_properties('N')
Nb = read_simulation_properties('Nb') #number of blade
Nstep = read_simulation_properties('Nstep') #number of time step per revoluiton
r = read_simulation_properties('r')
h = read_simulation_properties('h')
c = read_simulation_properties('c')
b = read_simulation_properties('b') #blade span [m]
Uref = read_simulation_properties('Uref')
Lambda = read_simulation_properties('Lambda')
beta0 = read_simulation_properties('beta0')
lawType = read_simulation_properties('lawType')
xA = read_simulation_properties('xl')
xB = read_simulation_properties('xc')
rho = read_simulation_properties('rho')

pdyn = 0.5 * rho * (Uref**2) #dynamique pressure
d = 2*r + h #height of the area swept by the blades

all_figure = False #if you want all figure to be print

##################### Calculate dynamic parameters ##############

xG_dim = xG_position_OO12(c) #[m] position of the center of gravity,
								#xG_dim=O leading edge position
								#xG_dim=c trailing edge position
xG = xG_dim/c #[]relative  position of the center of gravity,
								#xG_dim=O leading edge position
								#xG_dim=1 trailing edge position

JG = area_moment_inertia_OO12(xG_dim, c)[0] #[m^4] area moment of inertia along axe (G ez)

rhoBlade = 10000 #[kg/m³] density of blade material

IG = JG * rhoBlade * b #[kg.m²] moment of inertia along axe (G ez)

section = section_OO12(c) #[m²] blade profil section area
m = section * b * rhoBlade #[kg] mass of one blade

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
data = {}

for i in range(Nb):
	blade = 'blade' + str(i+1)
	data[blade] = {}

	path_postprocess_input = path_postProcessing + 'forceCoeffsBlade' + str(i+1)
	print(f"\nIs going to read {path_postprocess_input} containing {blade} data.\n")
	if not os.path.exists(path_postprocess_input):
		raise ValueError(f"No file {path_postprocess_input}")

	#List  time subdirectories in  Post-process Input directory :
	time_list = read_time_dir(path_postprocess_input)

	#Read  all time files and concatenate data :
	data[blade] = read_OF_data(path_postprocess_input, time_list, OF_version = '1912plus')

	#force coeff readen in postprocessing do not represent anymore lift and drag forces
	#because the CofR is moving throught time. they are rename as : 
	if 'Cd' in data[blade]: 
		data[blade]['Cx'] = data[blade].pop('Cd')
		data[blade]['fx'] = pdyn * b * c * data[blade]['Cx']
	else :
		raise ValueError(f"No key Cd in data file  {path_postprocess_input}")

	if 'Cl' in data[blade]:
		data[blade]['Cy'] = data[blade].pop('Cl')
		data[blade]['fy'] = pdyn * b * c * data[blade]['Cy']
	else:
		raise ValueError(f"No key Cl in data file  {path_postprocess_input}")

	'''
	CmPitch give the pitching moment calculated at reference frame origine point O(0,0)
	It value change if we calculate it at an other point (following the "BABAR" rule).

	So we define a key 'CmO' which is a subDict containing 3 keys :
	 'x':float=0,'y':float=0, 'value':float

	We will define later the key 'CmA' : the pitching moment coeff at the blade rotation point A,
	deduct from CmO with the BABAR rules
	'''

	if 'CmPitch' in data[blade] : 
		data[blade]['CmO'] = {}
		data[blade]['CmO']['value'] = data[blade].pop('CmPitch')
		size = len(data[blade]['CmO']['value'])
	#print('size='+str(size))
		data[blade]['CmO']['x'] = np.zeros(size)
		data[blade]['CmO']['y'] = np.zeros(size)
	else:
		raise ValueError(f"No key Cl in data file  {path_postprocess_input}")
	#create TimeAdim=t/T key
	data[blade]['TimeAdim'] = data[blade]['Time']/T

	print(f"{len(data[blade]['Time'])} time steps have been read from t="+
          f"{data[blade]['Time'][0]:4f}s to t={data[blade]['Time'][-1]:4f}" )





##################### add to data the s position parameter ##############
for blade in data :
	rank = int(blade.replace('blade','')) #rank == 1 for blade1, ect..
	#print('read information of '+blade+' , the blade at position '+str(rank))
	data[blade] = add_posParam(data[blade],h,r,T,Nb,rank)



#print('max s='+str(max(data[blade]['s']))+', min s='+str(min(data[blade]['s'])))


################## add postion of A point #####
for blade in data :
	data[blade] = add_CofR(data[blade], h, r) # the A point 

    #rename 'CofR' as 'A'
	#if 'CofR' in data[blade]: data[blade]['A'] = data[blade].pop('CofR')
	#else :
	#	raise ValueError(f"No key CofR in data file  {path_postprocess_input}")


	#print(f"max x={max(data[blade]['A']['x'])}, min x={min(data[blade]['A']['x'])}")
	#print(f"max y={max(data[blade]['A']['y'])}, min y={min(data[blade]['A']['y'])}")




################## add relative velocity of the fluid with respect to the blade ########

#data=add_Urel(data,h,r,Uref,Lambda)



################## add angle between x direction and A point trajectory's tangent ########
'''
-pi<gamma<=pi
'''

for blade in data : 
	data[blade] = add_gamma(data[blade],h,r)

	data[blade]['gammaTot'] = make_it_continuous(data[blade]['gamma'],-pi,pi)


################## add angle between x direction and chord ########
'''
-pi<theta<=pi
if lawType == 3 also add ['B']['x'] and ['B']['y'] to dict 
'''
for blade in data : 
	data[blade] = add_theta(beta0,c,xA,xB,data[blade],h,lawType,r)
	data[blade]['thetaTot']=make_it_continuous(data[blade]['theta'],-pi,pi)

	#print(data[blade]['B']['x'])
	#print(data[blade]['B']['y'])
	if not 'B' in data[blade] :
		raise ValueError(f"No key B in dict data")

################ add setting angle beta=theta-gamma ###########

for blade in data : 
	data[blade]=add_settingAngle(data[blade])


################## add rotation speed ###########

for blade in data : 
	data[blade]=add_omega(data[blade]) #omega=angular velovity=d/dt(totRot) [rad/s]

################## add sB ################
for blade in data : 
	data[blade] = add_sB(data[blade],h,r)
	#print(f"max sB={max(data[blade]['sB'])}, min s={min(data[blade]['sB'])}")

if all_figure:
	for blade in data :
		fig, axs = plt.subplots(figsize=(10,5))

		axs.plot(data[blade]['TimeAdim'],data[blade]['s'],label='s')
		axs.plot(data[blade]['TimeAdim'],data[blade]['sB'],label='sB')
		axs.set(xlabel='t/T', ylabel='position parameter')
		#axs.set_ylim(-200,200)

		axs.grid()
		axs.legend(loc='best')
		plt.show()
		plt.close(fig)

################## add G ################
'''
G is the center of gravity of the blade

'''
for blade in data : 
	data[blade] = add_G(data[blade], xA, xB, xG)

if all_figure:
	for blade in data :
		fig, axs = plt.subplots(2,2,figsize=(10,10))

		axs[0,0].plot(data[blade]['s'],data[blade]['CofR']['x'],label='xA')
		axs[0,0].plot(data[blade]['s'],data[blade]['B']['x'],label='xB')
		axs[0,0].plot(data[blade]['s'],data[blade]['G']['x'],label='xG')
		axs[0,0].set(xlabel='s', ylabel='x')
		axs[0,0].grid()
		axs[0,0].legend(loc='best')

		axs[1,0].plot(data[blade]['s'],data[blade]['CofR']['y'],label='yA')
		axs[1,0].plot(data[blade]['s'],data[blade]['B']['y'],label='yB')
		axs[1,0].plot(data[blade]['s'],data[blade]['G']['y'],label='yG')
		axs[1,0].set(xlabel='s', ylabel='y')
		axs[1,0].grid()
		axs[1,0].legend(loc='best')

		axs[0,1].plot(data[blade]['CofR']['x'], data[blade]['CofR']['y'], label='yA')
		axs[0,1].plot(data[blade]['B']['x'], data[blade]['B']['y'], label='yB')
		axs[0,1].plot(data[blade]['G']['x'], data[blade]['G']['y'], label='yG')
		axs[0,1].set(xlabel='x', ylabel='y')
		axs[0,1].grid()
		axs[0,1].legend(loc='best')

		plt.show()
		plt.close(fig)

################## add G accelaration ################
'''
G is the center of gravity of the blade

'''
for blade in data : 
	data[blade] = add_d2G_dt2(data[blade])
	#print(f"max d2G_dt2(x)={max(data[blade]['d2G_dt2']['x'])}, \
    #        min d2G_dt2(x)={min(data[blade]['d2G_dt2']['x'])}")
	#print(f"max d2G_dt2(y)={max(data[blade]['d2G_dt2']['y'])}, \
    #        min d2G_dt2(y)={min(data[blade]['d2G_dt2']['y'])}")

################## add dynamique moment in G ################
'''
G is the center of gravity of the blade
dynMomentG = d/dt(Ig.Omega).ez
'''
for blade in data : 
	data[blade] = add_dynMomentG(data[blade], IG)

################## add dynamique moment in A ################
'''
A = CofR
dynMomentA = dynMomentG + (AG ^ m.d²/dt²OG).ez
'''
for blade in data : 
	data[blade] = add_dynMomentA(data[blade], m)

if all_figure:
	for blade in data :
		fig, axs = plt.subplots(figsize=(10,5))

		axs.plot(data[blade]['TimeAdim'],data[blade]['dynMomentA'],label='dynamic moment in A')
		axs.plot(data[blade]['TimeAdim'],data[blade]['dynMomentG'],label='dynamic moment in G')
		axs.set(xlabel='t/T', ylabel='dynamic moment')
		#axs.set_ylim(-200,200)

		axs.grid()
		axs.legend(loc='best')
		plt.show()
		plt.close(fig)

################## add gammaB angle ################
for blade in data : 
	data[blade] = add_gammaB(data[blade],h,r)

	data[blade]['gammaBTot'] = make_it_continuous(data[blade]['gammaB'],-pi,pi)

if all_figure:
	for blade in data :
		fig, axs = plt.subplots(2,figsize=(10,10))

		axs[0].plot(data[blade]['TimeAdim'],data[blade]['gamma'],label='gamma A')
		axs[0].plot(data[blade]['TimeAdim'],data[blade]['gammaB'],label='gamma B')
		axs[0].set(xlabel='t/T', ylabel='gamma ')
		axs[0].grid()
		axs[0].legend(loc='best')

		axs[1].plot(data[blade]['TimeAdim'],data[blade]['gammaTot'],label='gamma continue A')
		axs[1].plot(data[blade]['TimeAdim'],data[blade]['gammaBTot'],label='gamma continue B')
		axs[1].set(xlabel='t/T', ylabel='gamma')
		axs[1].grid()
		axs[1].legend(loc='best')

		plt.show()
		plt.close(fig)


######### add MO=pitching moment  exerced by fluid on blade calculated at O point ########
for blade in data :
	data[blade]['MO'] = copy.deepcopy(data[blade]['CmO']) 
	data[blade]['MO']['value'] = - pdyn * b * c**2 * data[blade]['CmO']['value']

'''
!!!!!!!!! Here the "-" would mean that CmO is the coeff moment of blade on the fluid ...
'''

######### add MA=pitching moment  exerced by fluid on blade calculated at A point ########

for blade in data : 
	data[blade]=add_MA(data[blade])



######### add MA_old "old fashion way" ########

for blade in data : 
	data[blade] = add_CmA(data[blade], c)
	data[blade]['CmA_old'] = copy.deepcopy(data[blade]['CmA'])
	data[blade]['MA_old'] = copy.deepcopy(data[blade]['CmA_old'])
	data[blade]['MA_old']['value'] = pdyn * b * (c**2) * data[blade]['MA_old']['value']

######### add CmA=pitching moment coefficient exerced by fluid on blade calculated at A point ########
for blade in data : 
	data[blade]['CmA']=copy.deepcopy(data[blade]['MA'])
	data[blade]['CmA']['value'] = (1/(pdyn * b * (c**2))) * data[blade]['CmA']['value']

######### add fnB = normal force exerced by the chain on the blade at B liaison ########

for blade in data: 
	data[blade]=add_fnB(data[blade])


######### add ftA = tangential force exerced by the chain on the blade at A liaison ########
######### add fnA = normal force exerced by the chain on the blade at A liaison ########
'''
obtain in writing theorem of dynamic resultant on blade projected on direction tangential to 
A point trajectory
'''
for blade in data : 
	data[blade]=add_ftA_fnA(data[blade], m)

######### add fxA = force exerced by the chain on the blade at A liaison projected in x direction
######### add fyA = force exerced by the chain on the blade at A liaison projected in y direction
######### add fxB = force exerced by the chain on the blade at B liaison projected in x direction
######### add fyB = force exerced by the chain on the blade at B liaison projected in y direction
for blade in data : 
	data[blade]=add_fxA_fyA_fxB_fyB(data[blade])

######### add Effectif power coefficient ########
'''
CpEffectif = power coefficient transmit to the blade = FtA * V / (Uref * q * S)
'''
for blade in data : 
	data[blade]=add_CpEffectif(data[blade], Uref, Lambda, pdyn, b, d)


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
	#print('len(data[blade][CmA][value])=')
	#print(len(data[blade]['CmA']['value']))

	#print('len(data[Ct])=')
	#print(len(data[blade]['Ct']))

	#print('len(data[omega])=')
	#print(len(data[blade]['omega']))

	data[blade]=add_Cp_Cpt_Cpr(data[blade],c,h,r,Uref,Lambda)

	'''other way to calculate power coeff : CmO*omega
	pertinent only for Darrius (with fixed setting angle?)
	'''
	data[blade]=add_CpO(data[blade],c,h,r,Uref,Lambda)




################## return Cp every nT time #########################
'''
extract Cp mean on every cycle for every blade to see convergence cycle to cycle
'''

Cp_evol = {}
if data['blade1']['Time'][-1] > 2*T:
	for blade in data :
		Cp_evol[blade] = {}
		Cp_evol[blade]['value'] = extract_Cp_evolution(data[blade], T)
		#['value'] key contains the average Cp of each cycle
		#['variation'] key contains pourcentage of change of Cp average from one
		# cycle to another
		size = len(Cp_evol['blade1']['value'])
		Cp_evol[blade]['variation'] = np.zeros(size)
		for i in range(size-1):
			var = Cp_evol[blade]['value'][i+1] - Cp_evol[blade]['value'][i]
			var = 100 * var / Cp_evol[blade]['value'][i]
			Cp_evol[blade]['variation'][i+1] = var

	size = len(Cp_evol['blade1']['value'])
	Cp_evol_tot = {}
	Cp_evol_tot['value'] = np.zeros(size)
	Cp_evol_tot['variation'] = np.zeros(size)
	for blade in data :
		Cp_evol_tot['value'] = Cp_evol_tot['value'] + Cp_evol[blade]['value']

	for i in range(size-1):
		var = Cp_evol_tot['value'][i+1] - Cp_evol_tot['value'][i]
		var = 100 * var / Cp_evol_tot['value'][i]
		Cp_evol_tot['variation'][i+1] = var
	
	Cp_evol['total'] = {}
	Cp_evol['total'] = copy.deepcopy(Cp_evol_tot)




##################### extract only last period from  data ##############
'''
we check if the simulation has runned for more than one period :
'''
data_lastPer={}

if data['blade1']['Time'][-1]>T:
	meanGlobalCp = 0
	meanGlobalCpEffectif = 0
	for blade in data :
		data_lastPer[blade] = extract_period(data[blade], T, position='last')

		#we reorder data_lastPer by s ascending order:
		#data_lastPer[blade]=buble_sort_dict(data_lastPer[blade],'s')

		#time step being constant, the average Cp is simply the mean of Cp(t) 
		#value during one periode
		data_lastPer[blade]['meanCP'] = np.mean(data_lastPer[blade]['Cp'])
		data_lastPer[blade]['meanCPEffectif'] = np.mean(data_lastPer[blade]['CpEffectif'])
		meanGlobalCp += data_lastPer[blade]['meanCP']
		meanGlobalCpEffectif += data_lastPer[blade]['meanCPEffectif']

	
	#Here we create data_lastPer['global']['Cp'] which is the sum of Cp of every
	#blade with repect to time

	size = len(data_lastPer['blade1']['s'])
	Cp_tot = np.zeros(size)
	CpEffectif_tot = np.zeros(size)
	for blade in data :
		Cp_tot = Cp_tot + data_lastPer[blade]['Cp']
		CpEffectif_tot = CpEffectif_tot + data_lastPer[blade]['CpEffectif']

	data_lastPer['global']={}
	data_lastPer['global']['TimeAdim']=copy.deepcopy(data_lastPer['blade1']['TimeAdim'])
	data_lastPer['global']['s']=copy.deepcopy(data_lastPer['blade1']['s'])
	data_lastPer['global']['Cp']=copy.deepcopy(Cp_tot)	
	data_lastPer['global']['CpEffectif']=copy.deepcopy(CpEffectif_tot)	



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

for blade in data :
	fig, axs = plt.subplots(3, 4, figsize=(15,15))

	MA_adim = (1/max(data_lastPer[blade]['MA']['value'])) * data_lastPer[blade]['MA']['value']
	MA_old_adim = (1/max(data_lastPer[blade]['MA_old']['value'])) * data_lastPer[blade]['MA_old']['value']
	MO_adim = (1/max(data_lastPer[blade]['MO']['value'])) * data_lastPer[blade]['MO']['value']
	fnB_adim = (1/max(data_lastPer[blade]['fnB'])) * data_lastPer[blade]['fnB']

	som_fx = data_lastPer[blade]['fxA'] + data_lastPer[blade]['fxB'] + data_lastPer[blade]['fx'] 
	som_fy = data_lastPer[blade]['fyA'] + data_lastPer[blade]['fyB'] + data_lastPer[blade]['fy'] 



	axs[0,0].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['MA']['value'], label = 'MA')
	axs[0,0].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['MA_old']['value'], label = 'MA_old')
	axs[0,0].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['MO']['value'], label = 'MO')
	axs[0,0].set(xlabel='t/T', ylabel = 'moment ')
	axs[0,0].grid()
	axs[0,0].legend(loc='best')



	axs[1,0].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['CpEffectif'],
				label=f'global Effectif (Cp={round(meanGlobalCpEffectif, 7)})')
	axs[1,0].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['Cp'],
				label=f'global old (Cp={round(meanGlobalCpEffectif, 7)})')
	axs[1,0].set(xlabel='t/T', ylabel = 'Cp')
	axs[1,0].grid()
	axs[1,0].legend(loc='best')





	axs[0,1].plot(data_lastPer[blade]['TimeAdim'], data_lastPer[blade]['fxA'], label = 'fxA')
	axs[0,1].plot(data_lastPer[blade]['TimeAdim'], data_lastPer[blade]['fxB'], label = 'fxB')
	axs[0,1].plot(data_lastPer[blade]['TimeAdim'], data_lastPer[blade]['fx'], label = 'fxFluid')
	axs[0,1].plot(data_lastPer[blade]['TimeAdim'], som_fx, label = 'somme')
	axs[0,1].plot(data_lastPer[blade]['TimeAdim'], m * data_lastPer[blade]['d2G_dt2']['x'], label = 'm*acc(G).x')
	axs[0,1].set(xlabel='t/T', ylabel='f [N]')
	axs[0,1].grid()
	axs[0,1].legend(loc='best')

	axs[1,1].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['fyA'], label = 'fyA')
	axs[1,1].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['fyB'], label = 'fyB')
	axs[1,1].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['fy'], label = 'fyFluid')
	axs[1,1].plot(data_lastPer[blade]['TimeAdim'],som_fy, label = 'somme')
	axs[1,1].plot(data_lastPer[blade]['TimeAdim'], m * data_lastPer[blade]['d2G_dt2']['y'], label = 'm*acc(G).y')
	axs[1,1].set(xlabel='t/T', ylabel='f [N]')
	axs[1,1].grid()
	axs[1,1].legend(loc='best')

	axs[2,1].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['ftA'], label = 'ftA')
	axs[2,1].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['fnA'], label = 'fnA')
	axs[2,1].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['fnB'], label = 'fnB')
	axs[2,1].set(xlabel='t/T', ylabel='f [N]')
	axs[2,1].grid()
	axs[2,1].legend(loc='best')

	axs[0,2].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['CofR']['x'], label = 'xA')
	axs[0,2].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['CofR']['y'], label = 'yA')
	axs[0,2].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['B']['x'], label = 'xB')
	axs[0,2].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['B']['y'], label = 'yB')
	axs[0,2].set(xlabel='t/T', ylabel='position')
	axs[0,2].grid()
	axs[0,2].legend(loc='best')

	axs[1,2].plot(data_lastPer[blade]['CofR']['x'],data_lastPer[blade]['CofR']['y'], label = 'A')
	axs[1,2].plot(data_lastPer[blade]['B']['x'],data_lastPer[blade]['B']['y'], label = 'B')
	axs[1,2].plot(data_lastPer[blade]['G']['x'],data_lastPer[blade]['G']['y'], label = 'G')
	axs[1,2].set(xlabel='x', ylabel='y')
	axs[1,2].grid()
	axs[1,2].legend(loc='best')

	axs[2,2].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['gamma'], label = 'gammaA')
	axs[2,2].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['gammaB'], label = 'gammaB')
	axs[2,2].set(xlabel='t/T', ylabel='gamma [rad]')
	axs[2,2].grid()
	axs[2,2].legend(loc='best')


	axs[0,3].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['Cx'], label = 'Cx')
	axs[0,3].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['Cy'], label = 'Cy')
	axs[0,3].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['CmO']['value'], label = 'CmO')
	axs[0,3].set(xlabel='t/T', ylabel='force and moment Coeff')
	axs[0,3].grid()
	axs[0,3].legend(loc='best')

	axs[1,3].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['G']['x'], label = 'x')
	axs[1,3].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['G']['y'], label = 'y')
	axs[1,3].set(xlabel='t/T', ylabel='position of G')
	axs[1,3].grid()
	axs[1,3].legend(loc='best')


	axs[2,3].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['d2G_dt2']['x'], label = 'x')
	axs[2,3].plot(data_lastPer[blade]['TimeAdim'],data_lastPer[blade]['d2G_dt2']['y'], label = 'y')
	axs[2,3].set(xlabel='t/T', ylabel='acceleration of G')
	axs[2,3].grid()
	axs[2,3].legend(loc='best')

	plt.show()
	plt.close(fig)



'''
#plot Cp
for blade in data :
	#force coeff=f(t) 
	fig, axs = plt.subplots(figsize=(10,5))

	axs.plot(data[blade]['TimeAdim'],data[blade]['Cp'],'.',label='Cp')
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

#plot Cp_tot and CpEffectif_tot on last period

if 'global' in data_lastPer :
	#force coeff=f(t) 
	fig, axs = plt.subplots(figsize=(10,5))

	axs.plot(data_lastPer['global']['TimeAdim'],data_lastPer['global']['Cp'],
				label=f'global (Cp={round(meanGlobalCp, 7)})')

	axs.plot(data_lastPer['global']['TimeAdim'],data_lastPer['global']['CpEffectif'],
				label=f'global (Cp={round(meanGlobalCpEffectif, 7)})')

	data_lastPer['global']['CpEffectif']

	axs.set(xlabel='t/T', ylabel='Cp')

	axs.grid()
	axs.legend(loc='best')
	plt.savefig(os.path.join(path_postprocess_output,'Cp_lastPeriod.pdf'))
	plt.show()

	plt.close(fig)
'''
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

#Cp last 
wdata={}

wdata['Time']=copy.deepcopy(data_lastPer['blade1']['Time'])
wdata['TimeAdim']=copy.deepcopy(data_lastPer['blade1']['TimeAdim'])

size=len(wdata['Time'])
wdata['CpTotal']=np.zeros(size)

for blade in data_lastPer:
	if len(data_lastPer[blade]['Cp'])!=size:
		raise ValueError('len(Cp) from key ' +blade+'='+
			str(len(data_lastPer[blade]['Cp']))+' not egal to len(Time)='+str(size))
	wdata[blade]=copy.deepcopy(data_lastPer[blade]['Cp'])
	wdata['CpTotal']=wdata['CpTotal']+data_lastPer[blade]['Cp']


write_data(wdata,path_postprocess_output,'Cp_lastPeriod.txt',8)

wdata={}

for blade in Cp_evol:
	wdata[blade+'_value']=Cp_evol[blade]['value']
	wdata[blade+'_variation']=Cp_evol[blade]['variation']

write_data(wdata,path_postprocess_output,'CpEvol_CycleToCycle.txt',8)














