#!/usr/bin/env python

import os
import numpy as np
import runpy as rp
import copy
from math import *

def closest_rank_element(tab, t):
    """ take an array of float or int sorted in ascending order and a float t
            find the smallest element of tab superior or egal to t and return it position 
            return error if not(tab[0]<=t<=tab[-1])
    Args:
            tab: array[float]\n
            t: float\n

    Returns:
            rank: int\n

    A way you might use me is:\n
            i=closest_rank_element(my_tab,0.4)
    """
    size = len(tab)
    if size < 1:
        raise ValueError('closest_rank_element : empty array')

    for i in range(size-1):
        if tab[i] > tab[i+1]:
            raise ValueError('closest_rank_element : ' +
                             'The array is not sorted in ascending order')
    if tab[0] > t:
        raise ValueError('closest_rank_element : the smallest tab element' +
                         'is higher than t')
    if tab[-1] < t:
        raise ValueError('closest_rank_element : the highest tab element' +
                         'is smaller than t')

    rank = 0
    while tab[rank] < t and rank < size-1:
        rank += 1
    #print('closest_rank_element of t='+str(t)+':'+str(rank))
    # print('tab['+str(rank-1)+']='+str(tab[rank-1]))
    # print('tab['+str(rank)+']='+str(tab[rank]))

    # print('tab['+str(rank+1)+']='+str(tab[rank+1])+'\n')
    return(rank)


def non_sorted_closest_rank_element(tab, t):
    """ take an array of float or int  and a float t
            find the closest element and return it position and the difference
            between t and it closest element in absolute value
    Args:
            tab: array[float]\n
            t: float\n

    Returns:
            rank: int\n

    A way you might use me is:\n
            i, diff =non_sorted_closest_rank_element(my_tab,0.4)
    """

    rank = 0
    diff = abs(tab[0] - t)
    for i, elmt in enumerate(tab):
        if abs(elmt - t) < diff:
            rank = i
            diff = abs(elmt - t)
    return(rank, diff)
    
def will_it_int(s):
    """ Returns True is string is a int. """
    try:
        int(s)
        return True
    except ValueError:
        return False

def will_it_float(s):
    """ Returns True is string is a float. """
    try:
        float(s)
        return True
    except ValueError:
        return False

def RadToDeg(a):
	""" convert a angle from radian to degree.
	Args:
		a : float or int\n

	Returns:
		float\n """
	return (a*180/pi)

def DegToRad(a):
	""" convert a angle from degree to radian.
	Args:
		a : float or int\n

	Returns:
		float\n """
	return (a*pi/180)


def read_simulation_properties(parameter, path = '.', file = 'simulationProperties.py'):
    """ read the simulation properties file
            return the value of parameter
    Args:
            parameter: str\n
            path: str\n

    Returns:
            value: float\n

    A way you might use me is:\n
            myRe= read_simulation_properties('Re')

    """
    path_in_simuProperties = os.path.join(path, file)
    if not os.path.exists(path_in_simuProperties):
        raise ValueError('No file ' + path_in_simuProperties)

    d = rp.run_path(path_in_simuProperties, {})

    if parameter not in d:
        raise ValueError('The variable ' + parameter +
                         'has to be defined in the file ' + path_in_simuProperties)

    return (d[parameter])
    

def read_str_simulation_properties(str_parameter): 
	""" read the simulation properties file
		return the value of str_parameter 
	Args:
		parameter: str\n

	Returns:
		value: str\n

	A way you might use me is:\n
		myRe= read_simulation_properties('Re')

	"""
	path_in_simuProperties = os.path.join('.', 'simulationProperties.py')
	if not os.path.exists(path_in_simuProperties):
		raise ValueError('No file ' + path_in_simuProperties)

	file = open(path_in_simuProperties).readlines()
	isPresent=False
	for line in file:
		if str_parameter in line and not(line[0]=='#'):  
			line=line.split("'")
			isPresent=True

	value=line[1]      

	if not(isPresent):
		    raise ValueError('The str variable ' + str_parameter +
		                     'has to be defined in the file ' + path_in_simuProperties)

	return (value)



def buble_sort_dict(data, key):
    """ data is a dictionary composed of array or list of same size
            or subdict composed of array or list. 
            example : 
                    data={'T':[1.6,5],'alpha':[7.,3.],'CofR':{'x':[4.,6.],'y':[1.,-2.]}}
            it sorts the dictionnary keys in increasing order of key
            example with key='alpha': 
                    buble_sort_dict(data)={'T':[5,1.6],'alpha':[3.,7.],'CofR':{'x':[6.,4.],'y':[-2.,1.]}}

    Args:
            data: dict{array}\n
            key: str\n

    Returns:
            data_sorted: dict{array}\n

    A way you might use me is:\n
            data= buble_sort_dict(data,'alpha')

    """
    if key not in data:
        raise ValueError('buble_sort_dict:' +
                         'The key '+key+' is not define in data')

    for keys in data:
        if type(data[keys]) == dict:
            for subkeys in data[keys]:
                if len(data[key]) != len(data[keys][subkeys]):
                    raise ValueError('The array "' + subkeys +
                                     '" and '+key+' are not of the same size')
        elif len(data[key]) != len(data[keys]):
            raise ValueError('The array "' + keys +
                             '" and '+key+' are not of the same size')

    data_sort = copy.deepcopy(data)

    size = len(data_sort[key])

    for i in range(size):
        for j in range(size-i-1):
            if data_sort[key][j] > data_sort[key][j+1]:
                for keys in data_sort:
                    if type(data[keys]) == dict:
                        for subkeys in data[keys]:
                            tmp = data_sort[keys][subkeys][j]
                            data_sort[keys][subkeys][j] = data_sort[keys][subkeys][j+1]
                            data_sort[keys][subkeys][j+1] = tmp
                    else:
                        tmp = data_sort[keys][j]
                        data_sort[keys][j] = data_sort[keys][j+1]
                        data_sort[keys][j+1] = tmp

    return (data_sort)

def write_data(data,path,file_name,decimal): 
	""" data is a dictionary composed of array or list of same size. 
		example : 
			data={'Time':[1.6,5],'alpha':[7.,3.]}
		this function store data in a file call 'file_name' into the
			'path' folder
		decimal is the number of digit of the decimal part.
	Args:
		data: dict{array}\n
		path: str\n
		file_name: str\n
		decimal: int\n

	A way you might use me is:\n
		write_data(data,'./py_postprocess_data','mydata.txt',8)
	"""
	path_txt_file=os.path.join(path, file_name)
	DATA=[]
	my_header=''

	if 'Time' in data:
		DATA.append(data['Time'])
		my_header+=' Time '
	for key in data:
		if key!='Time':
			DATA.append(data[key])
			my_header+=' '+key+' '
	my_header=my_header[:-1]

	DATA=np.array(DATA)
	DATA=DATA.T

	with open(path_txt_file, 'w') as f:
		np.savetxt(f, DATA,delimiter=' ',header=my_header, fmt='%.'+str(decimal)+'f')



def read_data(path,file_name,delimiter=' '): 
	""" reads file containing data, example (with delimiter=' '):
			#data from simulation 234
			# alpha Cd  solver
			-6.38	-0.51   pimple
			-5.30	-0.43   simple
	and return them in a dictionnaty. In this example:
			 data={'alpha':[-6.38,-0.51],'Cd':[-5.30,-0.43],'solver':['pimple','simple']}
				type(data['alpha'])->numpy.array()
				type(data['solver'])->list


	the header part is compose of one or several lines with "#" in first character
	the last line of header must be the name of fields stored

	Args:
		path: str\n
		file_name: str\n
		delimiter: str\n

	Returns:
		data: dict{}\n
	

	A way you might use me is:\n
		data = read_data(path_expe,'Cl_data.txt')

	"""
	path_file= os.path.join(path,file_name)
	if not os.path.exists(path_file):
	   raise ValueError('No file ' + path_file)

	header=''
	l=[]


	f = open(path_file, "r")
	for x in f:
		if x[0]=='#':
			header=x
		else :
			if len(x)!=1:l.append(x)
	f.close()
	

	size=len(l)

	header=header.replace('"', '')
	header=header.replace('#', '')
	header=header.replace('\n', '')
	header=header.replace(')', '')
	header=header.replace('(', '')
	if delimiter==' ': keys=header.split()
	else:
		header=header.replace(' ', '')
		keys=header.split(delimiter)
	keys = [key for key in keys if key != '']
	nb_key=len(keys)

	#print('\nkeys=\n')
	#print(keys)

	data={}
	for key in keys:
		data[key]=[]

	#print('\nl[-1]=\n')
	#print(l[-1])
	#print('\nlen(l[-1])=\n')
	#print(len(l[-1]))

	for i in range(size):
		#if i==0 : print('delimiter=\n')
		#if i==0 : print(delimiter)
		#if i==0 : print('\nl[i]=\n')
		#if i==0 : print(l[i])
		l[i]=l[i].replace(')', '')
		l[i]=l[i].replace('(', '')
		if delimiter==' ': l[i]=l[i].split()
		else:
			l[i]=l[i].split(delimiter)
		#if i==0 : print('\nl[i]=\n')
		#if i==0 : print(l[i])
		#if i==0 : print(l[i])
		#print('i='+str(i))
		#print('nb_key='+str(nb_key))
		for j in range (nb_key):
			if will_it_float(l[i][j]):
				data[keys[j]].append(float(l[i][j]))
			else:
				data[keys[j]].append(l[i][j])

	for key in data : 
		if isinstance(data[key][0], float) : data[key]=np.array(data[key])
		#print(key)


	return (data)


def concatDict(data,N,popLast=False): 
	""" data is a dictionary composed of array or list or subdictionnary of array or list. 
		example : 
			data={'T':[1.6,5],'alpha':[7.,3.]}
		this function concatenate N time array with them self in a new dict.
		if N=1, the no modification is done
		if popLast==True, the last alement of arrays is not concatenated
	Args:
		data: dict{array}\n
		N: int\n

	A way you might use me is:\n
		newData=concatDict(data,N)
	"""
	
	data_concat = copy.deepcopy(data)


	

	if N>1:
		for i in range(N-1):
			for key in data:
				if type(data[key])==dict:
					for subkey in data[key]:
						if popLast==True :		
							data_concat[key][subkey]=np.concatenate((data_concat[key][subkey][:-1],data[key][subkey]))
						else :
							data_concat[key][subkey]=np.concatenate((data_concat[key][subkey],data[key][subkey]))
				else:
					if popLast==True :		
						data_concat[key]=np.concatenate((data_concat[key][:-1],data[key]))
					else :
						data_concat[key]=np.concatenate((data_concat[key],data[key]))


	return (data_concat)
				
			

def frame_change(x0,y0,gamma):
	""" (x0,y0) are component of a vector in R0 reference frame.
		this function give component of this vector in R1 reference
		frame. gamma is the angle in radian between x0 and x1. 
	Args:
		x0: float\n
		y0: float\n
		gamma: float\n

	Returns:
		x1: float\n
		y1: float\n

	A way you might use me is:\n
		(x1,y1)=frame_change(x0,y0,gamma)
	"""
	x1=cos(gamma)*x0  + sin(gamma)*y0
	y1=-sin(gamma)*x0 + cos(gamma)*y0
	

	return(x1,y1)


def firstDerCentered(phi,x):
	""" phi and x are 1D array of the same size.
		Return a 1D array containing the first derivative of phi 
		with respect to x. It use a second order centered scheme
	Args:
		phi: array\n
		x: array\n

	Returns:
		derPhi: array\n

	A way you might use me is:\n
		derPhi=firstDerCentered(Rot,Time)
	"""
	if len(phi)!=len(x) :
			raise ValueError('phi and x do not have the same lenght')

	if len(phi)<2 :
			raise ValueError('you can not calculated the derivate of'+
								' an array with less than two elements')

	size=len(phi)
	derPhi=np.zeros(size)
	
	derPhi[0]=(phi[1]-phi[0])/(x[1]-x[0])	
	
	for i in range(1,size-1):
			derPhi[i]=(phi[i+1]-phi[i-1])/(x[i+1]-x[i-1])

	imax=size-1
	derPhi[imax]=(phi[imax]-phi[imax-1])/(x[imax]-x[imax-1])


	return (derPhi)


def secondDerCentered(phi,x):
	""" phi and x are 1D array of the same size.
		Return a 1D array containing the second derivative of phi 
		with respect to x. It use a second order centered scheme
	Args:
		phi: array\n
		x: array\n

	Returns:
		derPhi: array\n

	A way you might use me is:\n
		derPhi=secondDerCentered(Rot,Time)
	"""
	if len(phi)!=len(x) :
			raise ValueError('phi and x do not have the same lenght')

	if len(phi)<2 :
			raise ValueError('you can not calculated the derivate of'+
								' an array with less than two elements')

	size=len(phi)
	derPhi=np.zeros(size)
	
	for i in range(1,size-1):
			derPhi[i]=(phi[i+1]-2*phi[i]+phi[i-1])/((x[i]-x[i-1])**2)

	derPhi[0]=derPhi[1]
	derPhi[-1]=derPhi[-2]


	return (derPhi)



def read_OF_scalar_data(path, file_name = None, OF_version='1912plus'):
    """ reads foam file containing scalar data with the format :
            "
            3
            (
            1.45
            4.5
            3.6
            )
            "
    return a 1D array and an integer which is the lenght of data
    Args:
            path: str\n
    file_name: str\n
    OF_version: str\n

    Returns:
            l: int\n
            x: array\n


    A way you might use me is:\n
            (lenght,myArray) = read_OF_scalar_data('./myPath','myfile',OF_version='1912plus')

    """
    x = []

    if OF_version != '1912plus':
        raise ValueError('OF version ' + OF_version +
                         ' : read_OF_scalar_data do not deal with this version')

    if file_name == None:
        path_motion_file = path
    else:
        path_motion_file = os.path.join(path, file_name)
    if not os.path.exists(path_motion_file):
        raise ValueError('No file ' + path_motion_file)

    if OF_version == '1912plus':

        file = open(path_motion_file).readlines()
        for line in file:
            x.append(line)

        x.pop(0)  # first line is just : ''
        l = int(x.pop(0))  # second line is data lenght
        x.pop(0)  # third line is opening parenthesis
        x.pop(-1)  # last line is closing parenthesis
        x = [float(elmt) for elmt in x]
        x = np.array(x)

        if l != len(x):
            raise ValueError('l='+str(l)+' is not egal to len(x)='+str(len(x)))

    return (l, x)



def read_OF_vector_data(path, file_name = None, OF_version='1912plus'):
    """ reads foam file containing vector data with the format :
            "
            3
            (
            (1.45 7.5 -9.3)
            (4.5  6.54 6.1)
            (3.6  3.2 7.1)
            )
            "
    return a integer which is the lenght of data and three 1D array
            containing rdata espectly from the first, second and trird column.
    Args:
            path: str\n
    file_name: str\n
    OF_version: str\n

    Returns:
            l: int\n
            x: array\n
            y: array\n
            z: array\n


    A way you might use me is:\n
            (lenght,myx,myy,myz) = read_OF_vector_data('./myPath','myfile',OF_version='1912plus')

    """
    data = []

    if OF_version != '1912plus':
        raise ValueError('OF version ' + OF_version +
                         ' : read_courant do not deal with this version')

    if file_name == None:
        path_motion_file = path
    else:
        path_motion_file = os.path.join(path, file_name)
    if not os.path.exists(path_motion_file):
        raise ValueError('No file ' + path_motion_file)

    if OF_version == '1912plus':

        file = open(path_motion_file).readlines()
        for line in file:
            data.append(line)

        data.pop(0)  # first line is just : ''
        l = int(data.pop(0))  # second line is data lenght
        data.pop(0)  # third line is opening parenthesis
        data.pop(-1)  # last line is closing parenthesis

        size = len(data)
        x = np.zeros(size)
        y = np.zeros(size)
        z = np.zeros(size)

        for i in range(size):
            data[i] = data[i].replace('(', '')
            data[i] = data[i].replace(')', '')
            data[i] = data[i].replace('\n', '')
            [x[i], y[i], z[i]] = data[i].split(' ')
            #data[i]=data[i].split(' ')

        x = [float(elmt) for elmt in x]
        y = [float(elmt) for elmt in y]
        z = [float(elmt) for elmt in z]

        if l != len(x):
            raise ValueError('l='+str(l)+' is not egal to len(x)='+str(len(x)))
        if l != len(y):
            raise ValueError('l='+str(l)+' is not egal to len(y)='+str(len(y)))
        if l != len(y):
            raise ValueError('l='+str(l)+' is not egal to len(z)='+str(len(z)))

    return (l, x, y, z)


def find_3orderPolCoeff(x0,f0,df0,x1,f1,df1):
	""" Calculate coefficient of a third order polynomial function : 
		f(x)=ax³+bx²+cx+d
		read several parameters :
		x0,x1 : boundary points position\n
		f0,f1 : value of the fonction at this points \n
		df0,df1  : value of the derivative at this points \n

	Returns coefficient (a,b,c,d).

	Args:
		x0 : float\n
		f0 : float\n
		df0 : float\n
		x1 : float\n
		f1 : float\n
		df1 : float\n


	Returns:
		(a,b,c,d) : floats\n

	A way you might use me is:\n
		(a,b,c,d)= find_3orderPolCoeff(x0,f0,df0,x1,f1,df1)

	"""
	#the system of equation is :

	#f0=ax0³+bx0²+cx0+d
	#f1=ax1³+bx1²+cx1+d
	#df0=3ax0²+2bx0+c
	#df1=3ax1²+2bx1+c
	
	#we put it in a matricial form : Ax=b, with x=[a,b,c,d]⁻1 and b=[f0,f1,df0,df1]⁻1
	A=np.array([[x0**3,x0**2,x0,1],[x1**3,x1**2,x1,1],[3*x0**2,2*x0,1,0],[3*x1**2,2*x1,1,0]])
	
	b=np.array([f0,f1,df0,df1])

	#and we solve the system x=A⁻1b with a np function :
	x = np.linalg.solve(A, b)
	a=x[0]
	b=x[1]
	c=x[2]
	d=x[3]

	return (a,b,c,d)


def solve_2ndOrderPolynom(a,b,c):
	""" Calculate zeros of a second order polynomial function : 
		f(x)=ax²+bx+c
	Returns error message if there is no solution or (x1,x2) the zeros
		of the polynom

	Args:
		a : float\n
		b : float\n
		c : float\n

	Returns:
		(x1,x2) : floats\n

	A way you might use me is:\n
		(x1,x2)= solve_2ndOrderPolynom(1,2,3.4)

	"""
	delta=b**2-4*a*c

	if delta<0:
	   raise ValueError('The second order polynome have no solution')

	x1=(-b-sqrt(delta))/(2*a)
	x2=(-b+sqrt(delta))/(2*a)

	return (x1,x2)



def how_many_intersec_2circles(x0,y0,r0,x1,y1,r1):
	""" Calculate the number of intersection points between 2 circles.
		first circle of radius r0 and center (x0,y0)
		second circle of radius r1 and center (x1,y1)
		for description of method check :
		https://stackoverflow.com/questions/3349125/circle-circle-intersection-points
	Returns the number of intersection (O,1,2 or 'inf')

	Args:
		x0 : float\n
		y0 : float\n
		r0 : float\n
		x1 : float\n
		y1 : float\n
		r1 : float\n

	Returns:
		nb : int or char\n

	A way you might use me is:\n
		nb_intersection= how_many_intersec_2circles(x0,y0,r0,x1,y1,r1)

	"""
	d=sqrt((x1-x0)**2+(y1-y0)**2)

	if d>r0+r1:
		#there are no solutions, the circles are separate. 
		nb=0

	elif d<abs(r0-r1):
		#then there are no solutions because one circle is contained within the other.
		nb=0

	elif d==0 and r0==r1:
		#then the circles are coincident and there are an infinite number of solutions.
		nb='inf'

	elif d==r0+r1:
		#the circle have just one commun point
		nb=1

	else:
		nb=2

	return (nb)


def coord_intersec_2circles(x0,y0,r0,x1,y1,r1):
	""" Calculate coordonate of intersection points between 2 circles.
		first circle of radius r0 and center (x0,y0)
		second circle of radius r1 and center (x1,y1)
		for description of method check :
		https://stackoverflow.com/questions/3349125/circle-circle-intersection-points
		https://www.ipgirl.com/65697/points-dintersection-du-cercle.html
	Returns error message if there is no intersection or an infinity.

	Args:
		x0 : float\n
		y0 : float\n
		r0 : float\n
		x1 : float\n
		y1 : float\n
		r1 : float\n

	Returns:
		(x3,y3,x4,y4) : floats\n

	A way you might use me is:\n
		(x3,y3,x4,y4)= coord_intersec_2circles(x0,y0,r0,x1,y1,r1)

	"""

	d=sqrt((x1-x0)**2+(y1-y0)**2)
	if d>r0+r1:
	   raise ValueError('there are no solutions, the circles are separate. ')

	elif d<abs(r0-r1):
	   raise ValueError('there are no solutions because one circle is contained within the other.')


	elif d==0 and r0==r1:
	   raise ValueError('the circles are coincident and there are an infinite number of solutions.')


	a = (r0**2-r1**2+d**2)/(2*d) 
	h=sqrt(r0**2-a**2)

	#P2 = P0 + a ( P1 - P0 ) / d 
	x2=x0 + a*(x1-x0)/d
	y2=y0 + a*(y1-y0)/d

	x3 = x2 + h*(y1-y0)/d
	y3 = y2 - h*(x1-x0)/d

	x4 = x2 - h*(y1-y0)/d
	y4 = y2 + h*(x1-x0)/d



	return (x3,y3,x4,y4)


def coord_intersec_circle_line(x0,y0,r0,x1,y1,x2,y2):
	""" Calculate coordonate of intersection points between a line and a cercle.
		Circle of radius r0 and center in P0(x0,y0)
		Line contain points P1(x1,y1) and P2(x2,y2)
	Returns the coordonates of intersection points.

	Args:
		x0 : float\n
		y0 : float\n
		r0 : float\n
		x1 : float\n
		y1 : float\n
		x2 : float\n
		y2 : float\n

	Returns:
		(x3,y3,x4,y4) : floats\n

	A way you might use me is:\n
		(x3,y3,x4,y4)= coord_intersec_circle_line(x0,y0,r0,x1,y1,x2,y2)

	"""

	if x1==x2 and y1==y2:
	   raise ValueError('P1 and P2 are coincident so they do not define a line ')

	#general form line equation y=ax+b.
	#circle equation (x-x0)²+(y-y0)²=r0². Intersection points respect both equation.

	if x1==x2:
		#line equation x=x1
		if r0<abs(x0-x1):
	   		raise ValueError('there are no solutions, the circle and the line are separate. ')
		else :
			x3=x1
			y3=y0+sqrt(r0**2-(x1-x0)**2)
			x4=x1
			y4=y0-sqrt(r0**2-(x1-x0)**2)

	elif y1==y2:
		#line equation y=y1
		if r0<abs(y0-y1):
	   		raise ValueError('there are no solutions, the circle and the line are separate. ')
		else :
			y3=y1
			x3=x0+sqrt(r0**2-(y1-y0)**2)
			y4=y1
			x4=x0-sqrt(r0**2-(y1-y0)**2)

	else:
		#line equation y=ax+b
	
		a=(y1-y2)/(x1-x2)
		b=(y1*x2-y2*x1)/(x2-x1)

		#we replace y in circle equation by ax+b (using line equation)
		#the 2nd order polynome Ax²+Bx+C give the x3 and x4 coordonate 
		A=(1+a**2)
		B=-2*x0+2*a*(b-y0)
		C=x0**2+(b-y0)**2-r0**2

		(x3,x4)= solve_2ndOrderPolynom(A,B,C)
		y3=a*x3+b
		y4=a*x4+b

	return (x3,y3,x4,y4)


def make_it_continuous(f,a,b):
	""" Transform list of scalar to make it continus.
		f is 1D array containing scalar values.
		a and b are scalar with a<b. 
		a and b are respectly the lower and upper bound of f.
		f is a continus function except for rank n0 where
			lim(n->n0-) f(n)=a et lim(n->n0+) f(n)=b
		 or lim(n->n0-) f(n)=b et lim(n->n0+) f(n)=a

	Returns g 1D array of continus values by translating each of
		them by Nfx*(b-a), whit Nfx the adapted integer fot he fx values. 

	Args:
		f : 1DArray(float)\n
		a : float\n
		b : float\n

	Returns:
		g : 1DArray(float)\n

	A way you might use me is:\n
		g=make_it_continuous(f,a,b)

	"""

	if a>=b:
	   	raise ValueError('scalar a is egal or superior to b, should strictly inferior.')

	size=len(f)
	g=np.zeros(size)

	nb=0 #count the discontinuity nb+=-1 when lim(n->n0-) f(n)=a et lim(n->n0+) f(n)=b
		 #count the discontinuity nb+=1 when lim(n->n0-) f(n)=b et lim(n->n0+) f(n)=a
	
	eps=(b-a)/50 #array values will be considered discontinus in n when 
						# fn<a+eps and fn+1>b-eps	or 	fn>b-eps and fn+1<a+eps

	g[0]=f[0]
	for i in range(size-1):
		if f[i]>b-eps and f[i+1]<a+eps:
			nb+=1
		elif f[i]<a+eps and f[i+1]>b-eps:
			nb+=-1
			
		g[i+1]=f[i+1]+nb*(b-a)
									

	return (g)



def shift_array(tab,x):
	""" Transform 1D array of scalar tab sort in ascending order.

		N is the first element of tab superior or egal to x.

		newTab=tab[N:-1]+tab[:N]+tab[N+1]


	Args:
		tab : array[float]\n
		x : float\n

	Returns:
		newTab : array[float]\n

	A way you might use me is:\n
		newTab=shift_data(tab,0.6)

	"""

	N=closest_rank_element(tab,x)

	a=tab[:N] 
	b=[tab[N]]
	c=tab[N:-1]



	newTab=np.concatenate((c,a,b))
		

	return (newTab)



def shift_data(data,mykey,x):
	""" Transform dictionnary composed of 1D array and subdict of 1D array.
		All arrays have the same size 

		data contains key mykey=1DArray of scalar sorted in ascending order.
		N is the first element of mykey superior or egal to x.

		All array of data and it sub dict will be re-ordered as follow :
			newTab=oldTab[N:-1]+oldTab[:N]+oldTab[N+1]

	Returns dictionnary containing same key as data.

	Args:
		data : {}\n
		mykey : char\n
		x : float\n

	Returns:
		dataBld : {}\n

	A way you might use me is:\n
		dataBld=shift_data(data,'s',0.5)

	"""
	dataBld= copy.deepcopy(data)

	N=closest_rank_element(dataBld[mykey],x)

	for key in dataBld:
		if type(dataBld[key])==dict:
			for subkey in dataBld[key]:

				a=dataBld[key][subkey][N:-1]
				b=[dataBld[key][subkey][N]]
				c=dataBld[key][subkey][:N] 

				dataBld[key][subkey]=np.concatenate((c,a,b))
		else:
			a=dataBld[key][N:-1]
			b=[dataBld[key][N]]
			c=dataBld[key][:N] 

			dataBld[key]=np.concatenate((c,a,b))
		

	return (dataBld)

def integer_in_interval(x_min,x_max):
	""" this fonction return 1D array of integer in the interval 
			[x_min,x_max]. We assume x_min<=x_max
		If there is no integer between x_min and x_max (example 2.3 and 2.5), 
		it return empty array.
			
	Args:
		x_min : float\n
		x_max : float\n

	Returns:
		I : array[int]\n\n

	A way you might use me is:\n
		myintegers=integer_in_interval(2.54,9.3)

	"""
	if x_min>x_max:
	   	raise ValueError('x_max should be superior or egal to w_min')

	I=[]

	int_max=floor(x_max) 	 #floor(5.2)==5 , floor(-8.2)==-9
	if (floor(x_min)<x_min):
		int_min=floor(x_min)+1
	else:
		int_min=floor(x_min)

	
	#there will be at least one element in I if x_min<=int_min<=int_max<=x_max
	
	while int_min<=int_max:
		I.append(int_min)
		int_min+=1
	
	I=np.array(I)

	return(I)


def vect_product(u,v):
	""" u=[ux,uy,uz] and v=[vx,vy,vz] are coordonates of 3D vector expressed
			in a direct and "orthonormé" base of R^3. 
		it return the vector egal to the vector product of u^v
	
	Args:
		u : array[float]\n
		v : array[float]\n

	Returns:
		w : array[float]\n

	A way you might use me is:\n
		[wx,wy,wz]=vect_product([1.,2.5,56],[8,6,3.4])
	"""

	sizeU=len(u)
	sizeV=len(v)

	if sizeU!=3:
	   	raise ValueError('u is not vector of dimension 3')
	if sizeV!=3:
	   	raise ValueError('V is not vector of dimension 3')

	wx=u[1]*v[2]-u[2]*v[1]
	wy=u[2]*v[0]-u[0]*v[2]
	wz=u[0]*v[1]-u[1]*v[0]

	w=np.array([wx,wy,wz])
	return(w)

def scalar_product(u,v):
	""" u and v are coordonates of  vector of dimension n expressed
			in a direct and "orthonormé" base of R^n. 
		it return float egal to the scalar product u.v
		
	Args:
		u : array[float]\n
		v : array[float]\n

	Returns:
		x : float\n

	A way you might use me is:\n
		x=vect_product([1.,2.5,56,2,5],[8,6,3.4,5,9])x
	"""

	sizeU=len(u)
	sizeV=len(v)

	if sizeU!=sizeV:
	   	raise ValueError('u and v not of the same size')
	
	x=0.0
	for i in range(sizeU):
		x=u[i]*v[i]

	return(x)


def rot_vect_ez(u,theta):
	""" u=[ux,uy,uz] are component of a 3D vector in R0 reference frame.
		this function give component of this vector after a rotation of theta 
		angle of z-axis.
	Args:
		u: array[float]\n
		theta: float\n

	Returns:
		v: array[float]\n

	A way you might use me is:\n
		v=rot_vect_ez(u,gamma)
	"""
	v=np.zeros(3)
	v[0]=cos(theta)*u[0]  - sin(theta)*u[1]
	v[1]=sin(theta)*u[0] + cos(theta)*u[1]
	v[2]=u[2]

	return(v)


def rot_array(x, y, theta):
	""" x, y are array containing coordonate of point in R0 reference frame.
		This function give new coordonates array after o rotation of theta [rad] 
		arround (O, z) axis.
	Args:
		x: array[float]\n
		y: array[float]\n
		theta: float\n

	Returns:
		x_new: array[float]\n
		y_new: array[float]\n

	A way you might use me is:\n
		x1, y1 = rot_array(x0, y0, gamma)
	"""
	if len(x) != len(y):
		raise ValueError('x and y not of the same size')
    
	
	x_new = [x_old*cos(theta) - y_old*sin(theta) for x_old, y_old in zip(list(x), list(y))]
	y_new = [x_old*sin(theta) + y_old*cos(theta) for x_old, y_old in zip(list(x), list(y))]
	
	return np.array(x_new), np.array(y_new)

def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def signed_angle_between(v1, v2):
    """ Returns the signed angle in radians between vectors 'v1' and 'v2'.
        Both vector are in the (x,y) plane
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    angle = np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
    
    if v1_u[0] * v2_u[1] - v1_u[1] * v2_u[0] < 0:
        return -angle
    else:
        return angle





def lin_interpolation(x, f, xo, verbose = False):
        """ x is 1D array of float in ascending order
            f is 1D array corresponding to evalution of a fonction for x values
            xo is the float we want to determine fo by linear interpolation
            x[0]<= xo <= x[-1]
        Args:
	        x: array[float]\n
	        f: array[]\n
	        x: float\n

        Returns:
	        fo: float\n

        A way you might use me is:\n
	        fo = lin_interpolation(x, f, xo)
        """
	        
        i = closest_rank_element(x, xo) # rank of the closest superior or egal element, 
                                        #if not x[0]<= xo <= x[-1]  will retrn an error
        if i == 0 and x[i] == xo:
                # xo is egal to the fisrt element of x, so we just return f[0]:
                return f[i]
        elif i ==0 and x[i] != xo:
                raise ValueError('xo is inferior to the smallest element of x: interpolation impossible')
        else:
                x1, f1, x2, f2 = x[i-1], f[i-1], x[i], f[i]

                if verbose == True: 
                        print(f"x1 = {x1} \t x2 = {x2} \t f1 = {f1} \t f2 = {f2} \t f0={f1 + ((xo - x1) / (x2 - x1)) * (f2 - f1)}")

                return f1 + ((xo - x1) / (x2 - x1)) * (f2 - f1)

def read_vector_probes(path, file_name):
    """ reads file containing multi points probes data of a vectoriel field,
     example:
      # Probe 0 (0 -0.8 0)
	  # Probe 1 (-0.3061 -0.7391 0)
	  # Probe 2 (-0.5657 -0.5657 0)
	  #    Probe      0    1    2    3    
      #     Time
         12.33     (7.120 -2.35 0)   (5.46 -2.14 0)     (5.435 -1.45 1.67)   
         12.43     (7.523 -1.85 0)   (4.25  5.14 0)     (5.6 -1.45 1.67)           
         12.53     (8.180 -5.022 0)  (5.46  5.14 0)    (5.623 -1.45 1.67)  
         
    and return them in a dictionnaty. In this example:
       data = {'Probe 0': {
                           'location': (0, -0.8, 0);
                           'Time': np.array([12.33, 12.43, 12.53]);
                           'x': np.array([1.120, 7.523, 8.180]);
                           'y': np.array([-2.35, -1.85, -5.022]);
                           'z': np.array([0, 0, 0])
                           }
                           
               'Probe 1': {
                           'location': (-0.3061, -0.7391, 0);
                           'Time': np.array([12.33, 12.43, 12.53]);
                           'x': np.array([5.46, 4.25,  5.46]);
                           'y': np.array([-2.14, 5.14, 5.14]);
                           'z': np.array([0,        0,    0])
                           }
                etc.
                }

    Args:
            path: str\n
            file_name: str\n
    Returns:
            data: dict{}\n


    A way you might use me is:\n
            data = read_vector_probes(my_path, 'U')

    """
    path_file = os.path.join(path, file_name)
    if not os.path.exists(path_file):
        raise ValueError('No file ' + path_file)

    data = {}

    f = open(path_file, "r")
    for x in f:
        if x[0] == '#' and '(' in x:
            # first lines who defined a probe name and location:
            x = x.replace('#', '')
            x = x.replace('\n', '')
            x = x.replace(')', '')
            name, location = x.split('(')
            
            # remove useless space from probe name:
            name = name.replace(' ', '')
            
            location = location.split(' ')
            data[name] = {'location': [float(a) for a in location]}
        elif '#' not in x:
            # this line contain data at a given Time:
            l = x.split('(')
            
            Time = l.pop(0)
            Time = float(Time.replace(' ',''))
            
            if len(data) != len(l):
                raise ValueError('number of probes stored in data do not '+
                                 'match probes readen the line.')
                                 
            for v, key in zip(l, data):
                v = v.split(')')[0] 
                v = v.split(' ')
                x, y, z = (float(i) for i in v)
                if 'x' in data[key] and 'Time' in data[key]:
                    #keys 'x', 'y', 'z' and 'Time' already exist:
                    data[key]['Time'].append(Time)
                    data[key]['x'].append(x)
                    data[key]['y'].append(y)
                    data[key]['z'].append(z)
                else:
                    data[key]['Time'] = [Time]
                    data[key]['x'] = [x]
                    data[key]['y'] = [y]
                    data[key]['z'] = [z]
            
       
    f.close()
    for key in data:
        for subkey in data[key]:
            data[key][subkey] = np.array(data[key][subkey])
    
    return (data)

def read_multi_file_vector_probes(path, file_name, time_list):
    """ reads files containing multi points probes data of a vectoriel field,
     example:
      # Probe 0 (0 -0.8 0)
	  # Probe 1 (-0.3061 -0.7391 0)
	  # Probe 2 (-0.5657 -0.5657 0)
	  #    Probe      0    1    2    3    
      #     Time
         12.33     (7.120 -2.35 0)   (5.46 -2.14 0)     (5.435 -1.45 1.67)   
         12.43     (7.523 -1.85 0)   (4.25  5.14 0)     (5.6 -1.45 1.67)           
         12.53     (8.180 -5.022 0)  (5.46  5.14 0)    (5.623 -1.45 1.67)  
         
    and return them in a dictionnaty. In this example:
       data = {'Probe 0': {
                           'location': (0, -0.8, 0);
                           'Time': np.array([12.33, 12.43, 12.53]);
                           'x': np.array([1.120, 7.523, 8.180]);
                           'y': np.array([-2.35, -1.85, -5.022]);
                           'z': np.array([0, 0, 0])
                           }
                           
               'Probe 1': {
                           'location': (-0.3061, -0.7391, 0);
                           'Time': np.array([12.33, 12.43, 12.53]);
                           'x': np.array([5.46, 4.25,  5.46]);
                           'y': np.array([-2.14, 5.14, 5.14]);
                           'z': np.array([0,        0,    0])
                           }
                etc.
                }

    Args:
            path: str\n
            file_name: str\n
            time_list: list[str]\n
    Returns:
            data: dict{}\n


    A way you might use me is:\n
            data = read_multi_file_vector_probes(my_path, 'U', my_time_list)

    """
    data = {}

    for time in time_list:
        path_folder = os.path.join(path, time)

        if not os.path.exists(os.path.join(path_folder, file_name)):
            raise ValueError('No file ' + os.path.join(path_folder, file_name))
        
          
        tmp_data = read_vector_probes(path_folder, file_name) 

        # concatenate tmp_data into data: 
        if time == time_list[0]: 
            data = copy.deepcopy(tmp_data)
        else:
            # Warning: it could happend that time write in files overlaps,
            # i.e. we check if the latest time writen in data is indeed smaller
            # than the smallest time in tmp_data:
            for key in tmp_data:
                if tmp_data[key]['Time'][0] > data[key]['Time'][-1]:
                    for field in tmp_data[key]:
                        if field != 'location':
                            data[key][field] = np.concatenate((data[key][field], 
                                                           tmp_data[key][field]))

                elif tmp_data[key]['Time'][-1] > data[key]['Time'][-1] :
                    # if the greatest time in d_data is inferior or egal to the 
                     # greatest time in data, we do nothing !
                    i = closest_rank_element(tmp_data[key]['Time'],
                                             data[key]['Time'][-1])
                    # i = postion of the smallest element of tmp_data[key]['Time'] 
                    # superior or egal to data[key]['Time'][-1]
                    # we want a strict superiority :
                    if tmp_data[key]['Time'][i] == data[key]['Time'][-1]: i+=1

                
                for field in tmp_data[key]:
                    if field != 'location':
                        data[key][field] = np.concatenate((data[key][field], 
                                                    tmp_data[key][field][i:]))

    return (data)



