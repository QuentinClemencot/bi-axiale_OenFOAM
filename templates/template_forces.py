#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file './system/forceCoeffs'. 
It contains drag and lift directions, blade dimension and fluid volumic mass 
in order to compute force coefficients acting on eatch blades.

This template use './simulationProperties.py' file as input '''

import os
import sys
sys.path.append(
    "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axial_turbine/py/func")
from preprocess_function_def import read_simulation_properties

path_OF_header = "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/\
bi_axial_turbine/templates/OF_header.txt" #header for openfoam file
#################### Read usefull parameters #############################

Aref = read_simulation_properties('Aref')
c = read_simulation_properties('c')
Nb = read_simulation_properties('Nb')
rho = read_simulation_properties('rho')
Uref = read_simulation_properties('Uref')


#################### Define Output directory #############################
path_out = os.path.join('.', 'system')


#################### Create forceCoeffs file #############################
with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()



for i in range(1, Nb + 1):

    txt += f'''\n
forcesBlade{i}
{{
    type        forces;

    libs ( "libforces.so" );

    writeControl   timeStep;
    writeInterval   1;

    log         yes;

    patches     ( blade{i}Wall);

    rho         rhoInf; 
    rhoInf      {rho}; 
    CofR        (0 0 0);
    pitchAxis   (0 0 1);

}}
'''

txt += '''

// ************************************************************************* //'''


path = os.path.join('system', 'forces')
print('Save forces in path ' + path)

with open(path, 'w') as f:
    f.write(txt)



