#!/usr/bin/env python

from general_function_def import *
from preprocess_function_def import *
import os
import numpy as np
from math import *
import copy
import decimal



def read_time_dir(path):
	""" read the directory name contained in path.
	return an ascending order list of time directory (str type) stored in path.
	Args:
		path: str\n

	Returns:
		list_time_dir_str_sorted

	A way you might use me is:\n
		list_time = list_time('path_postProcessing_forceCoeffs')

	"""
	list_time_dir_str = [d for d in os.listdir(path)]
	list_time_dir_number = []

	for s in list_time_dir_str:
		if will_it_int(s): list_time_dir_number.append(int(s))
		elif will_it_float(s):
			list_time_dir_number.append(decimal.Decimal(s))
			# print("s=")
			# print(s)
			# print("float(s)=")
			# print(float(s))
		elif s == 'constant':
			pass
		else:
			raise ValueError('Directory name ' + s +
		                                 'can not be turn in a float or an integer')

	list_time_dir_number_sorted = sorted(list_time_dir_number)
	list_time_dir_str_sorted = [str(x) for x in list_time_dir_number_sorted]

	return (list_time_dir_str_sorted)



def read_OF_data(path, time_list, OF_version='1912plus'):
    """ reads and concatenates OF data (forceCoeff, residual...) contained in time directories
	return a dictionnary
	Args:
		path: str\n
        time_list: list[str]\n
        OF_version: str\n

	Returns:
		dictionnary of array, example :	\n
		{'T':np.array, 'Cl':np.array,'Cd':np.array}

	A way you might use me is:\n
		data = read_OF_data(path_postprocess_input,time_list,OF_version='1912plus')

    """
    data = {}

    if OF_version != '1912plus' and OF_version != '2006plus':
        raise ValueError('OF version ' + OF_version +
                         ' : read_OF_data do not deal with this version')
    if OF_version == '1912plus' or OF_version == '2006plus':
        for d in time_list:
            path_time_dir = os.path.join(path, d)

            # first, verify that there is only 1 file in d directory
            fd = os.listdir(path_time_dir)
            fd.sort() #sort in alphabet order
            if not len(fd) == 1:
                if len(fd) == 2:
                    if "force.dat" in fd[0] and "moment.dat" in fd[1]:
                    # if probes is "forces" instead of "forceCoeffs", 2 files are generated in time directory:
                              # "force.dat" and "moment.dat". We will read both files.
                        pass
                    elif '_' in fd[1] and '_' not in fd[0]:
                        fd[0] = fd[1]
                    elif '_' in fd[0] and '_' not in fd[1]:
                        pass

                    else:
                        print("oups!")
                        raise ValueError(
                            'Not excatly one file in' + path_time_dir + f': {fd} files detected')
                elif len(fd) == 4 and "force_" in fd[1] and "moment_" in fd[3]:
                    # we may be in a situration with fd[0] = 'forces.dat', fd[1] = 'forces_*.dat', 
                    # fd[2] = 'moment.dat', fd[3] = 'moment_*.dat'
                    # 'forces_*.dat' and 'moment_*.dat' files contains the complete set of data so
                    # we want them to be readen (so to be fd[0] and and  fd[1] )
                    fd[0], fd[1] = fd[1], fd[3]
                	
                else:
                    print("goes here bizarrement!")
                    raise ValueError('Not excatly one file in' +
                                     path_time_dir + f': {fd} files detected')
            '''
                if len(fd) == 2:
    	            if '_' in fd[1] and '_' not in fd[0]:
                        fd[0] = fd[1]
                    elif '_' in fd[0] and '_' not in fd[1]:
                        fd[0] = fd[0]
                    else:
                        raise ValueError(
                            'Not excatly one file in' + path_time_dir)

                else:
                    raise ValueError('Not excatly one file in' + path_time_dir)
             '''
            # second, read data in d directory
            # print(fd[0])
            # print(path_time_dir)
            
            if len(fd) in [2, 4]   and "force" in fd[0] and "moment" in fd[1]:
                d_data = read_data(path_time_dir, fd[0], delimiter=' ')
            	# we rename key in "force_key" because moment and force have the
            	# same column name
                tmp_key_to_change = []
                for key in d_data:
            	    if key != "Time":
                        tmp_key_to_change.append(key)
                for key in tmp_key_to_change:
                    d_data[f'force_{key}'] = d_data.pop(key)

                d_tmp = read_data(path_time_dir, fd[1], delimiter=' ')
                for key in d_tmp:
                    if key != "Time":
                        d_data[f'moment_{key}'] = copy.deepcopy(d_tmp[key])
            	
            else:
                d_data=read_data(path_time_dir,fd[0],delimiter=' ')
            #print(f"keys read in time dir {d}:")
            #print([key for key in d_data])
            # Endly, concatenate all data 
            if d==time_list[0] : data = copy.deepcopy(d_data)
            else:
				
# Warning it could happend that d_data[Time][0]<=data[Time][-1]
                if d_data['Time'][0] > data['Time'][-1]:
                    for key in data:
                        data[key]=np.concatenate((data[key],d_data[key]))

                elif d_data['Time'][-1] > data['Time'][-1] :
# if the greatest time in d_data is inferior or egal to the greatest
# time in data, we do nothing !
                    i=closest_rank_element(d_data['Time'],data['Time'][-1])
# i= postion of the smallest element of d_data[Time] 
# superior or egal to data[Time][-1]
# we want a strict superiority :
                    if d_data['Time'][i]==data['Time'][-1]: i+=1

                    for key in data:
# print('###### Before Concatenate :#######')
# print('read_OF_data : Time size : '+str(len(data['Time'])))
# print('read_OF_data : Cl size : '+str(len(data['Cl'])))
                        data[key]=np.concatenate((data[key],d_data[key][i:]))
# print('###### After Concatenate :#######')
# print('read_OF_data : Time size : '+str(len(data['Time'])))
# print('read_OF_data : Cl size : '+str(len(data['Cl'])))

    return (data)




def extract_period(data, period, position='last', verbose=False):
	""" take an dictionnary containing 'Time' key and the period [s] of revolution
		return a dictionary with same keys than the input dictionnary 
			but with data of only one period. 
			Last period if position='last' (defauld value)
			Period just before last one if position='penultimate'
	Args:
		data: dict{array}\n
		period: float\n
		position: str\n
		verbose: False\n
	Returns:
		rank: dict{array}\n

	A way you might use me is:\n
		data_per=extract_last_period(data,period,position='penultimate')
	"""

	if 'Time' not in data:
		raise ValueError('The key "Time" is not define in data')

	if len(data['Time'])==0:
		raise ValueError('data dictionnary is empty')

	for key in data:
		if type(data[key])==dict :
			for subkey in data[key]:
							if len(data[key][subkey])!=len(data['Time']):
									raise ValueError('The array "' + subkey +
					             		'" and "Time" are not of the same size')
		else:
			if len(data[key])!=len(data['Time']):
				raise ValueError('The array "' + key +
					             '" and "Time" are not of the same size')
		
	data_per={}
	if verbose: print(f"Period = {period}\n")

	if position=='last':
		start=closest_rank_element(data['Time'],data['Time'][-1]-period)
		for key in data:
			if type(data[key])==dict :
				data_per[key]={}
				for subkey in data[key]:
					data_per[key][subkey]=data[key][subkey][start:]
			else :
				data_per[key]=data[key][start:]

	elif position=='penultimate':
		start=closest_rank_element(data['Time'],data['Time'][-1]-2*period)
		end=closest_rank_element(data['Time'],data['Time'][-1]-period)
		for key in data:
			if type(data[key])==dict :
				data_per[key]={}
				for subkey in data[key]:
					data_per[key][subkey]=data[key][subkey][start:end]
			else :
				data_per[key]=data[key][start:end]

	else:
		raise ValueError('"position" is neighter egal to "last"'+
						 'nor to "penultimate"')

	return (data_per)


def extract_adim_periods(data):
	""" take an dictionnary containing 'TimeAdim' key.
		return a list of dictionary with same keys than the input dictionnary 
			each dictionary will contain the data of a single period 
			(ie TimeAdim in [N, N + 1[ ) 
			The last element of the list will
	Args:
		data: dict{array}\n

	Returns:
		list: list[dict]\n

	A way you might use me is:\n
		list=extract_adim_periods(data)
	"""

	if 'TimeAdim' not in data:
		raise ValueError('The key "TimeAdim" is not define in data')

	if len(data['TimeAdim'])==0:
		raise ValueError('data dictionnary is empty')

	for key in data:
		if type(data[key])==dict :
			for subkey in data[key]:
							if len(data[key][subkey])!=len(data['Time']):
									raise ValueError('The array "' + subkey +
					             	'" and "TimeAdim" are not of the same size')
		else:
			if len(data[key])!=len(data['TimeAdim']):
				raise ValueError('The array "' + key +
					             '" and "TimeAdim" are not of the same size')
		
	tad_min=data['TimeAdim'][0]
	tad_max=data['TimeAdim'][-1]

	L_int=integer_in_interval(tad_min,tad_max)
	L_index=[]
	for a in L_int :
		L_index.append(closest_rank_element(data['TimeAdim'],a))
		# L_index contains  indice of data['TimeAdim'] array of closest element
		# of "integer adim time"(ie Tad=1,tad=2,tad=3...)
	#print('tad_min='+str(tad_min))
	#print('tad_max='+str(tad_max))
	#print('integer between tad_min and tad_max:')
	#print(L_int)

	cycle=[]

	size=len(L_int)

	if tad_min<L_int[0] : #example :if tad_min=0.05, we had [0.05,1[ as a cycle
		one_cycle={}
		for key in data:
			if type(data[key])==dict :
				one_cycle[key]={}
				for subkey in data[key]:
					one_cycle[key][subkey]=data[key][subkey][:L_index[0]-1]
			else :
				one_cycle[key]=data[key][:L_index[0]-1]
		cycle.append(copy.deepcopy(one_cycle))


	for i in range (size-1): #here, one_cycle contain complete period of time
									# [N,N+1[
		one_cycle={}

		start=L_index[i]
		end=L_index[i+1]-1

		for key in data:
			if type(data[key])==dict :
				one_cycle[key]={}
				for subkey in data[key]:
					one_cycle[key][subkey]=data[key][subkey][start:end]
			else :
				one_cycle[key]=data[key][start:end]
		cycle.append(copy.deepcopy(one_cycle))

	if tad_max>L_int[-1] : #example: if tad_max=1.95, we had [1,1.95] as a cycle
		one_cycle={}
		for key in data:
			if type(data[key])==dict :
				one_cycle[key]={}
				for subkey in data[key]:
					one_cycle[key][subkey]=data[key][subkey][L_index[-1]:]
			else :
				one_cycle[key]=data[key][L_index[-1]:]
		cycle.append(copy.deepcopy(one_cycle))



	return (cycle)


def add_G(data, xA, xB, xG):
	""" take an dictionnary containing 'CofR' and 'B' 
		return a copy of dictionnary with additionnal 'G' key
	Args:
		data: dict{array}\n
		xA: float\n
		xB: float\n
		xG: float\n

	Returns:
		data_up: dict{array}\n
	"""

	if 'CofR' not in data:
		raise ValueError('The key "CofR" is not define in data')

	if 'B' not in data:
		raise ValueError('The key "B" is not define in data')

	data_up = copy.deepcopy(data)
	data_up['G']={}

	size = len(data_up['CofR']['x'])
	data_up['G']['x'] = np.zeros(size) 
	data_up['G']['y'] = np.zeros(size) 

    #weight
	wA = (xB - xG)/(xB - xA)  
	wB = (xG - xA)/(xB - xA)   

	for i, (xa_coord, xb_coord) in enumerate(zip(data_up['CofR']['x'], data_up['B']['x'])):
		data_up['G']['x'][i] = wA * xa_coord + wB * xb_coord

	for i, (ya_coord, yb_coord) in enumerate(zip(data_up['CofR']['y'], data_up['B']['y'])):
		data_up['G']['y'][i] = wA * ya_coord + wB * yb_coord


	return (data_up)

def add_dynMomentG(data, IG):
	""" take an dictionnary containing 'omega' key 
		return a copy of dictionnary with additionnal 'dynMomentG' key 
	Args:
		data: dict{array}\n
		IG: float (moment of inertia [kg.m²])\n

	Returns:
		data_up: dict{array}\n
	"""

	if 'omega' not in data:
		raise ValueError('The key "omega" is not define in data')

	data_up = copy.deepcopy(data)

	cinetic_momentG = IG * data_up['omega']
	data_up['dynMomentG'] = firstDerCentered(cinetic_momentG,  data_up['Time'])

	return (data_up)


def add_dynMomentA(data, m):
	""" take an dictionnary containing keys:
            - dynMomentG
            - CofR
            - G
            - d2G_dt2
        and m : mass of the blade [kg]
		return a copy of dictionnary with additionnal 'dynMomentA' key 
	Args:
		data: dict{array}\n
		m: float\n

	Returns:
		data_up: dict{array}\n
	"""

	if 'dynMomentG' not in data:
		raise ValueError('The key "dynMomentB" is not define in data')
	if 'CofR' not in data:
		raise ValueError('The key "CofR" is not define in data')
	if 'G' not in data:
		raise ValueError('The key "G" is not define in data')
	if 'd2G_dt2' not in data:
		raise ValueError('The key "d2G_dt2" is not define in data')

	data_up = copy.deepcopy(data)

	size = len(data_up['dynMomentG'])
	data_up['dynMomentA'] =np.zeros(size)

	for i in range(size):
		coordxG = data_up['G']['x'][i]
		coordyG = data_up['G']['y'][i]
		coordxA = data_up['CofR']['x'][i]
		coordyA = data_up['CofR']['y'][i]
		vx = m * data_up['d2G_dt2']['x'][i]
		vy = m * data_up['d2G_dt2']['y'][i]

		[useless, useless2, addz] = vect_product([coordxG-coordxA, coordyG-coordyA, 0], [vx, vy, 0])

		data_up['dynMomentA'][i] = data_up['dynMomentG'][i] + addz
        

	return (data_up)


def add_fnB(data):
	""" take an dictionnary containing keys:
            - dynMomentA : dynamique moment of the blade around (A,z) axis [kg.m²]
            - gammaB : angle between xo axe and tangent to the trajectory at B point [rad]
            - MA : moment on axe (Az) )exerced by the fluid on the blade 
            - B : position of B point
            - CofR : position of A point
		return a copy of dictionnary with additionnal 'fnB' key 
	Args:
		data: dict{array}\n
	Returns:
		data_up: dict{array}\n
	"""

	if 'dynMomentA' not in data:
		raise ValueError('The key "dynMomentA" is not define in data')
	if 'gammaB' not in data:
		raise ValueError('The key "gammaB" is not define in data')
	if 'MA' not in data:
		raise ValueError('The key "MA" is not define in data')
	if 'B' not in data:
		raise ValueError('The key "B" is not define in data')
	if 'CofR' not in data:
		raise ValueError('The key "CofR" is not define in data')


	data_up = copy.deepcopy(data)

	MA_liaisonB = data_up['dynMomentA'] - data_up['MA']['value']
                  #moment around (A, z) axis exerced by the chain
                  #at point B liaison on the blade.
    
	size = len(MA_liaisonB)
	data_up['fnB'] = np.zeros(size)

	for i in range(size):

		coordxA = data_up['CofR']['x'][i]
		coordyA = data_up['CofR']['y'][i]
		coordxB = data_up['B']['x'][i]
		coordyB = data_up['B']['y'][i]

		vecAB = np.array([coordxB - coordxA,\
						coordyB - coordyA,\
						0])
		nb = np.array([-sin(data_up['gammaB'][i]), \
						cos(data_up['gammaB'][i]), \
						0])
              #nb=vector normal to the trajectory of B point

		[useless, useless2, div] = vect_product(vecAB, nb)

		data_up['fnB'][i] = MA_liaisonB[i] / div
        

	return (data_up)

def add_ftA_fnA(data, m):
	""" take an dictionnary containing keys:
            - gamma : angle between xo axe and tangent to the trajectory at A point [rad]
            - gammaB : angle between xo axe and tangent to the trajectory at B point [rad]
            - d2G_dt2 : acceleration of gravity center G
            - fnB : norm of the force exerced by the chain on the blade at B point.
                    It direction is normal to the trajectory of B point.  [N]
            - fx : force exerced by the fluid on the blade in xo direction [N]
            - fy : force exerced by the fluid on the blade in yo direction [N]
        and m: mass of the blade [kg]
		return a copy of dictionnary with additionnal 'ftA' key 
	Args:
		data: dict{array}\n
		m: float\n
	Returns:
		data_up: dict{array}\n
	"""

	if 'gamma' not in data:
		raise ValueError('The key "gamma" is not define in data')
	if 'gammaB' not in data:
		raise ValueError('The key "gammaB" is not define in data')
	if 'd2G_dt2' not in data:
		raise ValueError('The key "d2G_dt2" is not define in data')
	if 'fnB' not in data:
		raise ValueError('The key "fnB" is not define in data')
	if 'fx' not in data:
		raise ValueError('The key "fx" is not define in data')
	if 'fy' not in data:
		raise ValueError('The key "fy" is not define in data')

	data_up = copy.deepcopy(data)

	size = len(data['gamma'])
	data_up['ftA'] = np.zeros(size)
	data_up['fnA'] = np.zeros(size)

	for i in range(size):

		gammaA = data['gamma'][i]
		gammaB = data['gammaB'][i]
		fnB = data['fnB'][i]

		ta = np.array([cos(gammaA), sin(gammaA)])
		na = np.array([-sin(gammaA), cos(gammaA)])
		accG = np.array([data['d2G_dt2']['x'][i], data['d2G_dt2']['y'][i]])
		f = np.array([data['fx'][i], data['fy'][i]])
		
		data_up['ftA'][i] = m * np.dot(accG, ta) - np.dot(f, ta) + sin(gammaB - gammaA) * fnB
		data_up['fnA'][i] = m * np.dot(accG, na) - np.dot(f, na) - cos(gammaB - gammaA) * fnB

	return (data_up)

def add_fxA_fyA_fxB_fyB(data):
	""" take an dictionnary containing keys:
            - ftA 
            - fnA 
            - fnB
            - gamma
            - gammaB
		return a copy of dictionnary with additionnal 'fxA', 'fyA','fxB','fyB', keys
	Args:
		data: dict{array}\n
	Returns:
		data_up: dict{array}\n
	"""

	data_up = copy.deepcopy(data)

	size = len(data['ftA'])
	data_up['fxA'] = np.zeros(size)
	data_up['fyA'] = np.zeros(size)
	data_up['fxB'] = np.zeros(size)
	data_up['fyB'] = np.zeros(size)

	for i in range(size):

		gammaA = data['gamma'][i]
		gammaB = data['gammaB'][i]
		fnB = data['fnB'][i]
		fnA = data['fnA'][i]
		ftA = data['ftA'][i]
		
		data_up['fxA'][i] = cos(gammaA) * ftA - sin(gammaA) * fnA
		data_up['fyA'][i] = sin(gammaA) * ftA + cos(gammaA) * fnA
		
		data_up['fxB'][i] = - sin(gammaB) * fnB  #ftB = 0
		data_up['fyB'][i] =   cos(gammaB) * fnB

	return (data_up)


def add_CpEffectif(data, Uref, Lambda, pdyn, b, d):
	""" take a dictionnary containing ftA key and :
		    - Uref: reference velocity [m/s]
		    - Lambda: TSR []
		    - pdyn: reference dynamic pressure [Pa]
		    - d: height of the area swept by the blades [m]
		    - c: chord lenght [m]
		return a copy of dictionnary with additionnal 'CpEffectif' key 
	Args:
		data: dict{array}\n
		Uref: float\n
		Lambda: float\n
		pdyn: float\n
		c: float\n
		b: float\n
	Returns:
		data_up: dict{array}\n
	"""

	if 'ftA' not in data:
		raise ValueError('The key "ftA" is not define in data')

	data_up = copy.deepcopy(data)

	S = b * d
	V = Uref * Lambda
	#ftA is the force exerced by the chain on the blade 
    #tangentially to the A point trjectory
    # 
	data_up['CpEffectif'] =  -(V / (Uref * pdyn * S)) * data_up['ftA']

	return (data_up)




def add_alpha_values(data,AlphaMean,amplitude,omega,inverse_sign=False):
	""" take an dictionnary containing time step data and the period of oscillation
		add to keys to the dictionnary : 
			'alpha': pitch angle
			'alpha_dot': first derive of pitch angle with respect to time
	Args:
		data: dict{array}\n
		AlphaMean: float\n
		amplitude: float\n
		omega: float\n
		inverse_sign: bool\n

	Returns:
		data_up: dict{array}\n

	A way you might use me is:\n
		data_per=add_alpha_values(data,AlphaMean,amplitude,omega)
	"""

	if 'T' not in data:
		raise ValueError('The key "T" is not define in data')

	if len(data['T'])==0:
		raise ValueError('data dictionnary is empty')

	for key in data:
		if len(data[key])!=len(data['T']):
			raise ValueError('The array ' + key +
			                 'and "T" are not of the same size')


	size=len(data['T'])
	alpha=np.zeros(size)
	alpha_dot=np.zeros(size)

	for i in range(size):
		alpha[i]=AlphaMean+amplitude*sin(omega*data['T'][i])
		alpha_dot[i]=amplitude*omega*cos(omega*data['T'][i])

	data_up = copy.deepcopy(data)

	if inverse_sign:
		data_up['alpha']=[ -i for i in alpha ]
		data_up['alpha_dot']=[ -i for i in alpha_dot ]
	else:
		data_up['alpha']=alpha
		data_up['alpha_dot']=alpha_dot


	return (data_up)


def split_data(data):
	""" take an dictionnary containing time step data with a 'alpha_dot' key 
		(angular speed of pitching angle)
		add to keys to the dictionnary : 
			'alpha': pitch angle
			'alpha_dot': first derive of pitch angle with respect to time
	Args:
		data: dict{array}\n
		AlphaMean: float\n
		amplitude: float\n
		omega: float\n

	Returns:
		[data_inc,data_dec]: dict{array}[2]\n

	A way you might use me is:\n
		[data_inc,data_dec]=split_data(data)
	"""

	if 'alpha_dot' not in data:
		raise ValueError('The key "alpha_dot" is not define in data')


	for key in data:
		if len(data[key])!=len(data['alpha_dot']):
			raise ValueError('The array ' + key +
				             'and "data_dot" are not of the same size')

	if len(data['alpha_dot'])==0:
		raise ValueError('data dictionnary is empty')


	size=len(data['alpha_dot'])

	data_inc={}
	data_dec={}
	for key in data:
		data_inc[key]=[]
		data_dec[key]=[]	

	for i in range(size):
		if data['alpha_dot'][i]>=0:
				for key in data:
					data_inc[key].append(data[key][i])
		if data['alpha_dot'][i]<0:
				for key in data:
					data_dec[key].append(data[key][i])


	return ([data_inc,data_dec])


def compare_data(data1,data2):
	""" The purpose is to compare data1 and data2
		data1 and data2 are dictionnaries composed of array of the same size
			example : data1={'T':[1.6,5],'Cl':[7.,3.]}
				      data2={'T':[6,7],'Cl':[4.,3.],'alpha':[6.,2.]}
		all keys of data1 must be include in data2
	Args:
		data1: dict{array}\n
		data2: dict{array}\n

	Returns:
		data_comp: dict{dict}\n

	A way you might use me is:\n
		data_comp=compare_data(data1,data2)
	"""


	for key in data1:
		if not(key in data2):
			raise ValueError('The key ' + key +
				             ' is in data1 but not in data2')

	for key in data1:
		if len(data1[key])!=len(data2[key]):
			raise ValueError('The size of ' + key +
				             ' array is different between data1 and data2 :\n'+
							 str(len(data1[key])) +' vs '+ str(len(data2[key])) )

	data_comp={}
	for key in data1:
		data_comp[key]={'max':0.,'mean':0.,'max_rel':0.,'mean_rel':0.}	

	for key in data_comp:
		diff=np.absolute(data1[key]-data2[key])
		diff_max=np.max(diff)
		diff_mean=np.mean(diff)
		mean=np.mean(np.absolute(data1[key]))
		data_comp[key]['max']=diff_max
		data_comp[key]['mean']=diff_mean
		data_comp[key]['max_rel']=100*diff_max/mean
		data_comp[key]['mean_rel']=100*diff_mean/mean
		
	return (data_comp)


def read_courant(path,log,OF_version='1912plus'): 
	""" reads log file to extrat courant number max and mean for every time step
	return a dictionnary 
	Args:
		path: str\n
        time_list: list[str]\n
        OF_version: str\n

	Returns:
		dictionnary of array :	\n
		data['T'] 
		data['CFLmax']
		data['CFLmean']

	A way you might use me is:\n
		data = read_courant(path_directory_case,log,OF_version='1912plus')

	"""
	data={}
	data['Time']=[]
	data['CFLmax']=[]
	data['CFLmean']=[]

	if OF_version!='1912plus':
		raise ValueError('OF version ' + OF_version+
                         ' : read_courant do not deal with this version')

	path_log= os.path.join(path,log)
	if not os.path.exists(path_log):
	   raise ValueError('No file ' + path_log)

	if OF_version=='1912plus':
		file = open(path_log).readlines()
		for line in file:
			if 'Time = ' in line and not('ExecutionTime ' in line):   
				line=line.replace('Time = ', '')       
				line=line.replace('\n', '')                
				data['Time'].append(float(line))
			elif 'Courant Number mean: ' in line:
				line=line.replace('\n', '') 
				line=line.split(" ")
				data['CFLmax'].append(float(line[-1])) 
				data['CFLmean'].append(float(line[-3]))

		data['CFLmax'].pop(0)
		data['CFLmean'].pop(0)

		data['CFLmax']=np.array(data['CFLmax'])
		data['CFLmean']=np.array(data['CFLmean'])

	#print('len(data[CFLmax])='+str(len(data['CFLmax'])))
	#print('len(data[CFLmean])='+str(len(data['CFLmean'])))


	return (data)

def read_execution_time(path,log,OF_version='2006plus'): 
	""" reads log file to extrat information about mesh update execution time
		and pimple loop execution time of home modifying overPimpleFoam solver
		Warning : execution time are expressed in milliseconds
					but Time is given in second
	Args:
		path: str\n
        log: str\n
        OF_version: str\n

	Returns:
		dictionnary of array :	\n
		data['T'] 
		data['MeshUp']
		data['PimpleLoop']

	A way you might use me is:\n
		data = read_execution_time('./','log',OF_version='2006plus')

	"""
	data={}
	data['Time']=[]
	data['MeshUp']=[]
	data['PimpleLoop']=[]

	if OF_version!='2006plus':
		raise ValueError('OF version ' + OF_version+
                         ' : read_execution_time do not deal with this version')

	path_log= os.path.join(path,log)
	if not os.path.exists(path_log):
	   raise ValueError('No file ' + path_log)

	if OF_version=='2006plus':
		file = open(path_log).readlines()
		for line in file:
			if 'Time = ' in line and not('ExecutionTime ' in line):   
				line=line.replace('Time = ', '')       
				line=line.replace('\n', '')                
				data['Time'].append(float(line))
			elif 'Total mesh update time :' in line:
				line=line.replace('\n', '') 
				line=line.replace('Total mesh update time :', '')  
				line=line.replace(' ', '') 
				data['MeshUp'].append(float(line))
			elif 'Total pimple loop time :' in line:
				line=line.replace('\n', '') 
				line=line.replace('Total pimple loop time :', '')  
				line=line.replace(' ', '') 
				data['PimpleLoop'].append(float(line))

		data['Time']=np.array(data['Time'])
		data['MeshUp']=np.array(data['MeshUp'])
		data['PimpleLoop']=np.array(data['PimpleLoop'])

	return (data)

def add_posParam(data,h,r,T,Nb,rank):
	""" read data containing 'Time' key and geometric parameters :
		h : length of vertical translation\n
		r : radius of rotation \n
		T : period of revolution\n
		Nb : nombre de pale\n
		rank : nombre de pale\n
	return a dictionnary with same keys than data plus the key 's'
		which is a 1D array which gives postion parameter s (0<=s<1)
		It use the Time information, geometric parameter of the trajectory,
	    the velocity of the blade and rank of blade to deduct s values. 
	The point caracterize by s=0 is taken at the start of the downward 
		vertical translation. If Nb=3, for the blade of rank=0, s(t=0)=O;
		for the blade of rank=1, s(t=0)=1/3; for the blade of rank=2, s(t=0)=2/3

	Args:
		data: {}\n
		h: float\n
		r: float\n
		w: float\n
		T: float\n
		Nb: int\n
		rank: int\n

	Returns:
		dataUp : {}\n

	A way you might use me is:\n
		data = add_posParam(data,h,r,T,Nb,rank)

	"""
	if 'Time' not in data:
			raise ValueError(' data does not have "Time" key')

	dataUp = copy.deepcopy(data)
	size=len(dataUp['Time'])

	dataUp['s']=np.zeros(size)

	for i in range(size):
		t=dataUp['Time'][i]
		s0=(rank-1)/Nb #position at t=0
		s=((t+s0*T)%T)/T #lets s3(t) the posiion of the rank 3 blade.
						 #s3(t=0)=s0=s1(t=s0*T)
		dataUp['s'][i]=s
		

	return (dataUp)



def add_sB(data, h, r):
	""" read data containing 'B' key sub-dict containing 'x' and 'y' keys
	"""
	if 'B' not in data:
			raise ValueError(' data does not have "B" key')

	dataUp = copy.deepcopy(data)
	size = len(dataUp['B']['x'])  	
	dataUp['sB'] = np.zeros(size)
	
	for i in range(size):
		x = dataUp['B']['x'][i]
		y = dataUp['B']['y'][i]
		dataUp['sB'][i] = sFromXY(h, r, x, y)
	return dataUp
	
def sFromXY(h, r, x, y):
	"""return s position parameter from location (x,y)
	"""
	l = 2*h + 2*pi*r
	if 0 >= y >= -h and x < 0:
		s = -y/l 
	elif y < -h :
		theta = acos(-x/r)
		arc = theta*r
		s = (h + arc)/l
	elif 0 >= y >= -h and x > 0:
		s = 0.5 + (h + y)/l
	elif y > 0:
		theta = acos(x/r)
		arc = theta*r
		s = 0.5 + (h + arc)/l
	else:
		raise ValueError(f"sFromXY : x, y = {x}, {y} tuple unconsistent value")
		
	return s

def add_d2G_dt2(data):
	""" read data containing 'G' and 'Time' keys
        return copy of data with additionnal key 'd2G_dt2' containing acceleration
        of G point.
	"""
	if 'G' not in data:
			raise ValueError(' data does not have "G" key')
	if 'Time' not in data:
			raise ValueError(' data does not have "Time" key')
    
	dataUp = copy.deepcopy(data)
	dataUp['d2G_dt2'] = {}
	dataUp['d2G_dt2']['x'] = secondDerCentered(dataUp['G']['x'],dataUp['Time'])
	dataUp['d2G_dt2']['y'] = secondDerCentered(dataUp['G']['y'],dataUp['Time'])
	return dataUp


def add_Ct_Cn(data):
	""" read data containing 'gamma', 'Cx', 'Cy' keys 
		gamma=angle between x direction and CofR trajectory's tangent
		gamma is in radian
	return a dictionnary with same keys than data plus the keys 'Ct'
		and 'Cn' which are  1D arrays which gives force coefficient
		respectly in tangential to the trajectory direction and
		normal to the trajectory direction.
		

	Args:
		data: {}\n

	Returns:
		dataUp : {}\n

	A way you might use me is:\n
		data = add_Ct_Cn(data)

	"""
	if 'gamma' not in data:
			raise ValueError(' data does not have "gamma" key')
	if 'Cx' not in data:
			raise ValueError(' data does not have "Cx" key')
	if 'Cy' not in data:
			raise ValueError(' data does not have "Cy" key')

	dataUp = copy.deepcopy(data)
	size=len(dataUp['s'])

	dataUp['Ct']=np.zeros(size)
	dataUp['Cn']=np.zeros(size)

	for i in range(size):
		Cx=dataUp['Cx'][i]
		Cy=dataUp['Cy'][i]
		gamma=dataUp['gamma'][i]

		(dataUp['Ct'][i],dataUp['Cn'][i])=frame_change(Cx,Cy,gamma)
		

	return (dataUp)



def read_motion(path,file_name,OF_version='1912plus'): 
	""" reads motion.dat file
	return a dictionnary 
	Args:
		path: str\n
        file_name: str\n
        OF_version: str\n

	Returns:
		dictionnary of array :	\n
		data['Time'] \n
		data['dx']\n
		data['dy']\n
		data['drot']\n

	A way you might use me is:\n
		data = read_motion('./contant','motion.dat',OF_version='1912plus')

	"""
	data={}
	data['Time']=[]
	data['x-x0']=[]
	data['y-y0']=[]
	data['theta-theta0']=[]
	data['dx']=[]
	data['dy']=[]
	data['dtheta']=[]
	data['d2x']=[]
	data['d2y']=[]
	data['d2theta']=[]

	if OF_version!='1912plus':
		raise ValueError('OF version ' + OF_version+
                         ' : read_courant do not deal with this version')

	path_motion_file= os.path.join(path,file_name)
	if not os.path.exists(path_motion_file):
	   raise ValueError('No file ' + path_motion_file)

	if OF_version=='1912plus':

		file = open(path_motion_file).readlines()
		for line in file:
			indices=[]
			line=line.replace('(', '') 
			line=line.replace(')', '') 
			line=line.replace('\n', '') 
			line=line.split(" ")
			for i in range(len(line)) : 
				if line[i]=='' : indices.append(i)
			for j in sorted(indices, reverse=True):
   				 del line[j]
			if len(line)>2:   
				line = [float(item) for item in line]
				data['Time'].append(line[0])
				data['x-x0'].append(line[1])
				data['y-y0'].append(line[2])
				data['theta-theta0'].append(line[6])

		for key in data:
			data[key]=np.array(data[key])

		print('len(data[Time])=')
		print(len(data['Time']))
		print('len(data[x-x0])=')
		print(len(data['x-x0']))

		data['dx']=firstDerCentered(data['x-x0'],data['Time'])
		data['dy']=firstDerCentered(data['y-y0'],data['Time'])
		data['dtheta']=firstDerCentered(data['theta-theta0'],data['Time'])

		data['d2x']=secondDerCentered(data['x-x0'],data['Time'])
		data['d2y']=secondDerCentered(data['y-y0'],data['Time'])
		data['d2theta']=secondDerCentered(data['theta-theta0'],data['Time'])

	return (data)



def add_omega(data):
	""" read data containing 'thetaTot' 
		(angle between x direction and chord in radian) and 'Time' keys.
	return a dictionnary with same keys than data plus the keys 'omega'
		which is  1D arrays which gives angular velocity as the time
		derivative of thetaTot.
		A second order centered scheme is used : 
			omega(t)=[thetaTot(t+dt)-thetaTot(t-dt)]/2dt
		For value at t=0, omega(0)=[thetaTot(dt)-totRot(0)]/dt
		For value at t=tmax, omega(tmax)=[thetaTot(tmax)-thetaTot(tmax-dt)]/dt
		Time step dt is supposed to be constant

	Args:
		data: {}\n

	Returns:
		dataUp : {}\n

	A way you might use me is:\n
		data = add_omega(data)

	"""
	if 'thetaTot' not in data:
			raise ValueError(' data does not have "thetaTot" key')
	if 'Time' not in data:
			raise ValueError(' data does not have "Time" key')

	dataUp = copy.deepcopy(data)
	size=len(dataUp['thetaTot'])

	dataUp['omega']=np.zeros(size)
	dataUp['omega']=firstDerCentered(dataUp['thetaTot'],dataUp['Time']) 

	return (dataUp)




def add_CmA(data,c):
	""" read data containing 'CmO', 'Cy', 'Cx', 'CofR' keys
	return a dictionnary with same keys than data plus the keys 'CmA'
		 which is a subDict containing 3 keys :\n
 			'x':float=data['CofR']['x']\n
			'y':float=data['CofR']['y']\n
			'value':float\n
		  A is the blade center of rotation point (CofR)
		CmA is the pitching moment coeffient calculated at A point.
		c is the chord lenth (needed to applay babar formula from force coeff)

	Args:
		data: {}\n
		c: float\n

	Returns:
		dataUp : {}\n

	A way you might use me is:\n
		data = add_CmA(data)

	"""
	if 'CmO' not in data:
			raise ValueError(' data does not have "CmO" key')
	if 'Cy' not in data:
			raise ValueError(' data does not have "Cy" key')
	if 'Cx' not in data:
			raise ValueError(' data does not have "Cx" key')
	if 'CofR' not in data:
			raise ValueError(' data does not have "CofR" key')

	dataUp = copy.deepcopy(data)
	size=len(dataUp['Cy'])

	dataUp['CmA']={}
	dataUp['CmA']['x']=np.zeros(size)
	dataUp['CmA']['y']=np.zeros(size)
	dataUp['CmA']['value']=np.zeros(size)

	for i in range (size):
		dataUp['CmA']['x'][i]=	dataUp['CofR']['x'][i]
		dataUp['CmA']['y'][i]=	dataUp['CofR']['y'][i]

		CmO=dataUp['CmO']['value'][i]
		xO=dataUp['CmO']['x'][i]
		yO=dataUp['CmO']['y'][i]
		xA=	dataUp['CofR']['x'][i]
		yA=	dataUp['CofR']['y'][i]
		Cx=	dataUp['Cx'][i]
		Cy=	dataUp['Cy'][i]

		#dataUp['CmA']['value'][i]=CmO+(xO-xA)*Cy/c-(yO-yA)*Cx/c
		dataUp['CmA']['value'][i]=-CmO+(xO-xA)*Cy/c-(yO-yA)*Cx/c


	return (dataUp)

def add_MA(data):
	""" read data containing 'MO', 'fy', 'fx', 'CofR' keys,
        calcul the moment exerced by the fluid on the blade at point A from the 
        value of this moment at point O (use BABAR rules)
	return a dictionnary with same keys than data plus the keys 'MA'
		 which is a subDict containing 3 keys :\n
 			'x':float=data['CofR']['x']\n
			'y':float=data['CofR']['y']\n
			'value':float\n
	
	Args:
		data: {}\n

	Returns:
		dataUp : {}\n

	A way you might use me is:\n
		data = add_MA(data)

	"""
	if 'MO' not in data:
			raise ValueError(' data does not have "MO" key')
	if 'fy' not in data:
			raise ValueError(' data does not have "fy" key')
	if 'fx' not in data:
			raise ValueError(' data does not have "fx" key')
	if 'CofR' not in data:
			raise ValueError(' data does not have "CofR" key')

	dataUp = copy.deepcopy(data)
	size=len(dataUp['MO']['value'])

	dataUp['MA']={}
	dataUp['MA']['x']=np.zeros(size)
	dataUp['MA']['y']=np.zeros(size)
	dataUp['MA']['value']=np.zeros(size)

	for i in range (size):
		xa = dataUp['CofR']['x'][i]
		ya = dataUp['CofR']['y'][i]
		MO = dataUp['MO']['value'][i]

		vecAO = np.array([-xa, -ya, 0]) #point O coordonate = (0, 0)
		vecf = np.array([dataUp['fx'][i], dataUp['fy'][i], 0]) 
								#vecf = force exerced by fluid on the blade

		dataUp['MA']['value'][i] = MO + vect_product(vecAO, vecf)[2]
		dataUp['MA']['x'][i] = xa
		dataUp['MA']['y'][i] = ya

	return (dataUp)

def add_Cp_Cpt_Cpr(data,c,h,r,Uref,Lambda):
	""" read data containing 'Ct', 'CmA','omega' keys and several arguments :
			c : chord lenght \n
			h : length of vertical translation\n
			r : radius of rotation \n
			Uref : upstream velocity magnitude\n
			Lambda : TSR=V/Uref, with V magnitude of the CofR in fixed ref. frame\n
	return a dictionnary with same keys than data plus the keys 'Cp'
		 which is a 1D array containing the power coefficient 

	Args:
		data: {}\n
		c : float\n
		h : float\n
		r : float \n
		Uref : float\n
		Lambda : float\n

	Returns:
		dataUp : {}\n

	A way you might use me is:\n
		data = add_Cp_Cpt_Cpr(data,c,h,r,Uref)

	"""
	if 'Ct' not in data:
			raise ValueError(' data does not have "Ct" key')
	if 'CmA' not in data:
			raise ValueError(' data does not have "CmA" key')
	if 'omega' not in data:
			raise ValueError(' data does not have "omega" key')

	dataUp = copy.deepcopy(data)
	size=len(dataUp['Ct'])

	dataUp['Cp']=np.zeros(size)
	dataUp['Cpt']=np.zeros(size)
	dataUp['Cpr']=np.zeros(size)

	d=h+2*r #height of the area swept by the center of rotation.
			#Attention: due to the width of the blade and its pitch angle, 
				#the total height is greater !	
 
	for i in range (size):
		
		Ct=data['Ct'][i]
		CmA=data['CmA']['value'][i]
		omega=data['omega'][i]

		dataUp['Cpt'][i]=Ct*c*Lambda/d
		dataUp['Cpr'][i]=CmA*omega*(c**2)/(Uref*d)

		dataUp['Cp'][i]=dataUp['Cpt'][i]+dataUp['Cpr'][i]

	return (dataUp)


def add_Cp(data,c,h,r,Uref,Lambda):
	""" read data containing 'Cat' key and several arguments :
			c : chord lenght \n
			h : length of vertical translation\n
			r : radius of rotation \n
			Uref : upstream velocity magnitude\n
			Lambda : TSR=V/Uref, with V magnitude of the CofR in fixed ref. frame\n
	return a dictionnary with same keys than data plus the keys 'Cp'
		 which is a 1D array containing the power coefficient 

	Args:
		data: {}\n
		c : float\n
		h : float\n
		r : float \n
		Uref : float\n
		Lambda : float\n

	Returns:
		dataUp : {}\n

	A way you might use me is:\n
		data = add_Cp(data,c,h,r,Uref)

	"""
	if 'Cat' not in data:
			raise ValueError(' data does not have "Cat" key')


	dataUp = copy.deepcopy(data)
	size=len(dataUp['Cat'])

	dataUp['Cp']=np.zeros(size)

	d=h+2*r #height of the area swept by the center of rotation.
			#Attention: due to the width of the blade and its pitch angle, 
				#the total height is greater !	
 
	for i in range (size):
		
		Cat=data['Cat'][i]
		dataUp['Cp'][i]=Cat*c*Lambda/d

	return (dataUp)

def add_CpO(data,c,h,r,Uref,Lambda):
	""" read data containing 'CmA' keys and several arguments :
			c : chord lenght \n
			h : length of vertical translation\n
			r : radius of rotation \n
			Uref : upstream velocity magnitude\n
			Lambda : TSR=V/Uref, with V magnitude of the CofR in fixed ref. frame\n
	return a dictionnary with same keys than data plus the keys 'CpO'
		 which is a 1D array containing the power coefficient 
	Warning : calcul ok only for Darrieus turbine with constant setting angle

	Args:
		data: {}\n
		c : float\n
		h : float\n
		r : float \n
		Uref : float\n
		Lambda : float\n

	Returns:
		dataUp : {}\n

	A way you might use me is:\n
		data = add_CpO(data,c,h,r,Uref)

	"""
	if 'CmA' not in data:
			raise ValueError(' data does not have "CmA" key')

	omega=Lambda*Uref/r

	dataUp = copy.deepcopy(data)
	size=len(dataUp['CmA']['value'])
	dataUp['CpO']=np.zeros(size)

	d=h+2*r #height of the area swept by the center of rotation.
			#Attention: due to the width of the blade and its pitch angle, 
				#the total height is greater !	
 
	for i in range (size):
		
		CmO=data['CmO']['value'][i]

		dataUp['CpO'][i]=-CmO*omega*(c**2)/(Uref*d)

	return (dataUp)



	
def add_totRot_cumul_postProcess(data):
	""" read data containing 'totRot' key.
	return a dictionnary with same keys than data plus a modify dotRot key.
		which is a 1D array which gives the z-component of rotation vector 
		between the chord in a position s and the s=0 in degree.
		The value of totRot is not any more inferior to 360° but has a continue value.
		During the second blade revolution, 360°<totRot<=720°
	The point caracterize by s=0 is taken at the start of the downward 
		vertical translation.

	Args:
		data: {}\n


	Returns:
		dataUp : {}\n

	A way you might use me is:\n
		data = add_totRot_cumul(data)

	"""
	if 'totRot' not in data:
			raise ValueError(' data does not have "totRot" key')

	dataUp = copy.deepcopy(data)
	size=len(dataUp['totRot'])

	rev=0 #number of revolution
	for i in range(1,size):
		if dataUp['s'][i]<dataUp['s'][i-1]: rev+=1
		dataUp['totRot'][i]=dataUp['totRot'][i]+360*rev
		

	return (dataUp)



def extract_Cp_evolution(data, T, verbose=False):
	""" take an dictionnary containing 'Time' and 'Cp' key and the period T [s] 
		of revolution
		return a 1D array containing the average value of Cp for every cycle.
		len()=number of cycle, the last array value is the average cp from
		tmax-T to tmax, the penultimate value is the average cp from tmax-2T to
		tmax-T.

	Args:
		data: dict{array}\n
		T: float\n
		verbose: booleen\n

	Returns:
		Cp_evol: array[float]\n

	A way you might use me is:\n
		Cp_evol=extract_Cp_evolution(data,period)
	"""

	if 'Time' not in data:
		raise ValueError('The key "Time" is not define in data')

	if 'Cp' not in data:
		raise ValueError('The key "Time" is not define in data')

	if len(data['Time'])==0:
		raise ValueError('data dictionnary is empty')
	
	t_max=data['Time'][-1]
	t_min=data['Time'][0]

	if t_max<T:
		raise ValueError('Can not compute average Cp on cycles because data'
			+'[Time] max is inferior to cycle period T')
	
	if len(data['Cp'])!=len(data['Time']):
		raise ValueError('The array "Cp" and "Time" are not of the same size')
	
	nb_cycle=int((t_max-t_min)//T)
	Cp_evol=np.zeros(nb_cycle)

	print(f"Number of cycles read = {nb_cycle}")
	if verbose: print(f"from tmin={t_min} to tmax={t_max}")

	for i in range(nb_cycle):
		if verbose: print(f"i={i};\tt_max={t_max};\tt_min={t_min}\n")
		if (nb_cycle-1==i):
			end=len(data['Time']-1)
		else:	
			end=closest_rank_element(data['Time'],t_max-(nb_cycle-1-i)*T)
		start=closest_rank_element(data['Time'],t_max-(nb_cycle-i)*T)

		if verbose: print(f"start={start};\tend={end}")
		cp_cycleI=data['Cp'][start:end]
		Cp_evol[i]=np.mean(cp_cycleI)
	
	return (Cp_evol)


def s_from_position(data,point,h,r):
	""" read data containing 'point' key and geometric parameters :
		h : length of vertical translation\n
		r : radius of rotation \n
		'point' is a subdict containing keys 'x' and 'y' (1D arrays which give 
		the postion point)
	return a 1D array given the position parameter s of the point.
	The point caracterize by s=0 is taken at the start of the downward 
		vertical translation.

	Args:
		data: {}\n
		point: string\n
		h: float\n
		r: float\n

	Returns:
		s : 1D array[float]\n

	A way you might use me is:\n
		s = s_from_position(data,'B',h,r)

	"""
	if point not in data:
			raise ValueError(' data does not have '+point+' key')

	size=len(data[point]['x'])
	s=np.zeros(size)
	l=2*h+2*pi*r

	for i in range(size):
		x=data[point]['x'][i]
		y=data[point]['y'][i]

		if x<0 and y<=0 and y>=-h: #down translation part
			s[i]=-y/l

		elif   y<-h: #lower bend
			'''
			x0=-r
			y0=-h
			chord=sqrt((x0-x)**2+(y0-y)**2) #distance between point (x,y) and 
											#the begining of the bend (x0,y0)
			arc=2*asin(0.5*chord/r)*r#arc lenght between point (x,y) and (x0,y0)
			'''
			x1=-x
			y1=-y+h
			
			theta=acos(x1/r)
			arc=r*theta
			s[i]=(h+arc)/l
			
		elif x>0 and y<=0 and y>=-h: #down translation part
			s[i]=0.5+(h+y)/l

		elif   y>0: #lower bend
			x1=x
			y1=y
			
			theta=acos(x1/r)
			arc=r*theta
			s[i]=0.5+(h+arc)/l	
		

	return (s)



def s_from_one_position(x, y, h, r):
	""" return s position parameter corresponding to the point (x, y)
	in a turbine of radius r and translation h.

	Args:
		x, y, h, r float\n

	Returns:
		s : 1D array[float]\n

	A way you might use me is:\n
		s = s_from_position(data,'B',h,r)

	"""
	
	l = 2 * h + 2 * pi * r


	if x<0 and y<=0 and y>=-h: #down translation part
		s=-y/l
	elif   y<-h: #lower bend
		x1=-x
		theta=acos(x1/r)
		arc=r*theta
		s=(h+arc)/l
		
	elif x>0 and y<=0 and y>=-h: #down translation part
		s=0.5+(h+y)/l

	elif   y>0: #lower bend
		x1=x
		theta=acos(x1/r)
		arc=r*theta
		s=0.5+(h+arc)/l	
		
	return (s)
	
	
def add_Cbn(data,c):
	""" read data containing 'gammaB', 'CofR', 'B', 'CmA' keys 
		gammaB=angle between x direction and B point trajectory's tangent
		gamma is in radian
	return a dictionnary with same keys than data plus the keys 'Ct'
		and 'Cn' which are  1D arrays which gives force coefficient
		respectly in tangential to the trajectory direction and
		normal to the trajectory direction.
		

	Args:
		data: {}\n

	Returns:
		dataUp : {}\n

	A way you might use me is:\n
		data = add_Ct_Cn(data)

	"""
	if 'gammaB' not in data:
			raise ValueError(' data does not have "gammaB" key')
	if 'CofR' not in data:
			raise ValueError(' data does not have "CofR" key')
	if 'B' not in data:
			raise ValueError(' data does not have "B" key')
	if 'CmA' not in data:
			raise ValueError(' data does not have "CmA" key')

	dataUp = copy.deepcopy(data)
	size=len(dataUp['s'])

	Cbn=np.zeros(size)
	Cbx=np.zeros(size)
	Cby=np.zeros(size)

	CmA=copy.deepcopy(dataUp['CmA']['value'])
	gammaB=copy.deepcopy(dataUp['gammaB'])

	for i in range(size):

		ey=np.array([0,1,0]) #en is the normal vector to trajectory at B point
		en=rot_vect_ez(ey,gammaB[i])

		vec_AB=np.zeros(3) #vector AB
		vec_AB[0]=dataUp['B']['x'][i]-dataUp['CofR']['x'][i]
		vec_AB[1]=dataUp['B']['y'][i]-dataUp['CofR']['y'][i]

		x=vect_product(vec_AB,en)[2] #ez component of vectoriel product	

		Cbn[i]=c*CmA[i]/x
		CbLocal=np.array([0,Cbn[i],0])
		CbGlobal=np.array([])
		CbGlobal=rot_vect_ez(CbLocal,-gammaB[i])
		Cbx[i]=-Cbn[i]*sin(gammaB[i])
		Cby[i]=Cbn[i]*cos(gammaB[i])

		if i==2300:
			print('at i=:'+str(i))
			print('sB:')
			print(dataUp['sB'][i])
			print('gammaB:')
			print(gammaB[i])
			print('frame_change(0,1,gammaB[i])')
			print(frame_change(0,1,gammaB[i]))
			print('en:')
			print(en)
			print('vec_AB:')
			print(vec_AB)
			print('produit vectoriel vec_AB^en:')
			print(x)
			print('Cbn:')
			print(Cbn[i])
			print('CmA:')
			print(CmA[i])
			print('Cbx:')
			print(Cbx[i])
			print('Cby:')
			print(Cby[i])

	dataUp['Cbn']=Cbn
	dataUp['Cbx']=Cbx
	dataUp['Cby']=Cby

	return (dataUp)



def add_Cat(data,h,r):
	""" read data containing 'Cx', 'Cy', 's', 'Cbx', 'Cby','gamma' keys 
		
	return a dictionnary with same keys than data plus the keys 'Cat'
		which is  1D arrays which gives force coefficient
		in tangential direction at A point.
		

	Args:
		data: {}\n

	Returns:
		dataUp : {}\n

	A way you might use me is:\n
		data = add_Cat(data)

	"""
	if 'Cx' not in data:
			raise ValueError(' data does not have "Cx" key')
	if 'Cy' not in data:
			raise ValueError(' data does not have "Cy" key')
	if 's' not in data:
			raise ValueError(' data does not have "s" key')
	if 'Cbx' not in data:
			raise ValueError(' data does not have "Cbx" key')
	if 'Cby' not in data:
			raise ValueError(' data does not have "Cby" key')
	if 'gamma' not in data:
			raise ValueError(' data does not have "gamma" key')

	dataUp = copy.deepcopy(data)
	size=len(dataUp['s'])



	Cx=copy.deepcopy(dataUp['Cx'])
	Cy=copy.deepcopy(dataUp['Cy'])
	s=copy.deepcopy(dataUp['s'])
	Cbx=copy.deepcopy(dataUp['Cbx'])
	Cby=copy.deepcopy(dataUp['Cby'])
	gamma=copy.deepcopy(dataUp['gamma'])

	Cat=np.zeros(size)
	Cax=Cx-Cbx
	Cay=Cy-Cby



	for i in range(size):
		Cat[i]=frame_change(Cax[i],Cay[i],gamma[i])[0]


	dataUp['Cat']=Cat
	dataUp['Cax']=Cax
	dataUp['Cay']=Cay

	return (dataUp)


def reorder_points(data):
	'''data is a dict containing 'faceCentres' key containing points coordonates
       of a blade profil. 
       'faceCentres' is a subdict containing 'x0' and 'y0' keys (1D arrays).
       return a deep copy of data with faceCentres reorganized so it print nicely.
	'''
	#import matplotlib.pyplot as plt
	dataUp = copy.deepcopy(data)
	xA = dataUp['pointA']['x0']
	yA = dataUp['pointA']['y0']
	
	size = len(dataUp['faceCenter']['x0'])
	x = dataUp['faceCenter']['x0']
	y = dataUp['faceCenter']['y0']
	p = dataUp['p']
	wss_x =  dataUp['wss']['x0']
	wss_y =  dataUp['wss']['y0']
	
	#apply translation to make A coincide with center of reference frame 
	x1 = x - xA
	y1 = y - yA
	'''
	x2, y2 = rot_array(x1, y1, -pi)
	# apply rotation to make the chord coincide with x axis
	fig, axs = plt.subplots()
	axs.plot(x1,y1)
	axs.plot(x2, y2)
	axs.set_xlabel('x')
	axs.set_ylabel('y')
	axs.grid()
	axs.set_aspect('equal', 'box')
	plt.show()
	plt.close(fig)
	'''
	x1, y1 = rot_array(x1, y1, -dataUp['gamma']-pi) #leading edge at xmin
	#trailing egde at xmax
	#make min(x1) = O and max (x1) = 1:
	x1 = x1 - min(x1)
	x1 = (1/max(x1)) * x1
    # split point, and data associeted, upper and down part
	data_unsorted = {'x':x1, 'y':y1, 'x0':x, 'y0':y, 'p':p, 'wss_x':wss_x, 'wss_y':wss_y}
	data_upper = {'x':[], 'y':[], 'x0':[], 'y0':[], 'p':[], 'wss_x':[], 'wss_y':[]}
	data_down = {'x':[], 'y':[], 'x0':[], 'y0':[], 'p':[], 'wss_x':[], 'wss_y':[]}
	for i in range(len(x)):
		if y1[i] >= 0 :
			for key in data_unsorted : data_upper[key].append(data_unsorted[key][i])   	
		elif y1[i] <  0 :
			for key in data_unsorted : data_down[key].append(data_unsorted[key][i]) 
		else:
			raise ValueError(f'point (x={x[i]},y={y[i]}) do not belong to any subpart ') 

	#sort eatch part:
	#    - data_trailing in y ascending order
	#    - data_upper in x ascending order
	#    - data_down in x descending order
	data_upper = buble_sort_dict(data_upper, 'x')
	data_down = buble_sort_dict(data_down, 'x')
	for l in data_down.values():
		l.reverse()

	dataUp['x_adim'] = np.array(data_upper['x']+ data_down['x'])		
	dataUp['faceCenter']['x0'] = np.array(data_upper['x0']+ data_down['x0'])
	dataUp['faceCenter']['y0'] = np.array(data_upper['y0']+ data_down['y0'])
	dataUp['p'] = np.array(data_upper['p']+ data_down['p'])
	dataUp['wss']['x0'] = np.array(data_upper['wss_x']+ data_down['wss_x'])						 
	dataUp['wss']['y0'] = np.array(data_upper['wss_y']+ data_down['wss_y'])				
	if len(dataUp['faceCenter']['x0']) != size:
		raise ValueError('dataUp[faceCenter][x0] have changed it size!')
	return dataUp
	
	
def read_xflr5_polar(path, file_name):
    """ reads file containing xflr5 polar data, example:
        xflr5 v6.50

 	Calculated polar for: NACA 2p37_515

 	1 1 Reynolds number fixed          Mach number fixed         

 	xtrf =   0.050 (top)        0.050 (bottom)
 	Mach =   0.000     Re =     0.322 e 6     Ncrit =   9.000

  	alpha     CL        CD       CDp       Cm    Top Xtr Bot Xtr   Cpmin    Chinge    XCp    
 	------- -------- --------- --------- -------- ------- ------- -------- --------- ---------
 	-16.500  -0.5059   0.17952   0.17352  -0.0201  0.0500  0.0473  -1.9738   0.0000   0.1835
 	-13.500  -0.9457   0.03995   0.03250  -0.0898  0.0500  0.0406  -5.8948   0.0000   0.1354
 	-13.000  -0.9394   0.03615   0.02845  -0.0893  0.0500  0.0419  -5.7232   0.0000   0.1363
    and return them in a dictionnaty. In this example:
                     data={'alpha':[-16.5,-13.5, -13],'CL':[...]}

    Args:
            path: str\n
            file_name: str\n
            delimiter: str\n

    Returns:
            data: dict{}\n


    A way you might use me is:\n
            data = read_xflr5_polar(path_expe,'my_polar.txt')

    """
    path_file = os.path.join(path, file_name)
    if not os.path.exists(path_file):
        raise ValueError('No file ' + path_file)

    still_header = True
    header = []
    l = []

    f = open(path_file, "r")
    for x in f:
        if  still_header == True and not '-----' in x:
            header.append(x)
        elif '-----' in x:
            still_header = False
        else:
            l.append(x)
    f.close()

    
    size = len(l)
    keys = ['alpha', 'CL', 'CD', 'CDp', 'Cm', 'TopXtr', 'BotXtr', 'Cpmin', 'Chinge', 'XCp']
    nb_key = len(keys)
    data = {}
    for key in keys:
        data[key] = []
    for i in range(size):
        l[i] = l[i].split()
        if len(l[i]) > 1: #last lines can be empty
            for j in range(nb_key):
                if will_it_float(l[i][j]):
                    data[keys[j]].append(float(l[i][j]))
                else:
                    data[keys[j]].append(l[i][j])
    for key in data:
        if isinstance(data[key][0], float):
            data[key] = np.array(data[key])

    return (data)

