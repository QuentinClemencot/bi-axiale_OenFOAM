#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate  files:
       -'./0_orig/include/initialConditions'
       -'./constant/transportProperties'

    This template use './simulationProperties.py' file as input.
'''

import os
from math import sqrt 
from runpy import run_path
import sys
sys.path.append(
    "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axial_turbine/py/func")
from preprocess_function_def import read_simulation_properties

path_OF_header = "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/\
bi_axial_turbine/templates/OF_header.txt" #header for openfoam file

#################### Read usefull parameters #############################
c = read_simulation_properties('c')
Cmu = read_simulation_properties('Cmu')
I = read_simulation_properties('I')
Lturb = read_simulation_properties('Lturb')
nu = read_simulation_properties('nu')
Pref = read_simulation_properties('Pref')
Re = read_simulation_properties('Re')
Uref = read_simulation_properties('Uref')
Lambda = read_simulation_properties('Lambda')


if (nu == 0 and Re == 0) or (nu != 0 and Re != 0):
    raise ValueError('nu and Re cannot both be equal to zeros '+
                     'or both different from zeros')

####################### Calculating of nu ################################
if nu == 0 and Re != 0:
    W = Uref * sqrt(1 + Lambda**2) #relative speed of the blade with respect to 
								   #the flow in verticale translation parts.
    nu = W * c/Re  
    nu = round(nu, 15)

####################### Calculating k and omega ###########################

# Here is the source of information:
# https://www.openfoam.com/documentation/guides/latest/doc/guide-turbulence-ras-k-omega-sst.html


turbulentKE = 1.5 * (Uref * I)**2
turbulentKE = round(turbulentKE, 12)

turbulentOmega = (turbulentKE**0.5) / ((Cmu**0.25) * Lturb)
turbulentOmega = round(turbulentOmega, 12)

####################### writing ./constant/transportProperties ############

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()


txt += f'''

FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      transportProperties;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * // 


transportModel  Newtonian;

nu               [ 0 2 -1 0 0 0 0 ] {nu};

// ************************************************************************* //
'''



path = os.path.join('./constant', 'transportProperties')
print('Save transportProperties in path ' + path)

with open(path, 'w') as f:
    f.write(txt)


####################### writing ./0.orig/include/initialConditions ########

with open(path_OF_header, 'r', encoding = 'utf8') as f:
	txt = f.read()


txt += f'''


flowVelocity         ({Uref} 0 0);
pressure             {Pref};
turbulentKE          {turbulentKE};
turbulentOmega       {turbulentOmega};

// ************************************************************************* //
'''



path = os.path.join('0_orig/include', 'initialConditions')
print('Save initialConditions in path ' + path)

with open(path, 'w') as f:
    f.write(txt)

