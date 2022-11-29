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


#################### Read usefull parameters #############################

Aref = read_simulation_properties('Aref')
c = read_simulation_properties('c')
Nb = read_simulation_properties('Nb')
rho = read_simulation_properties('rho')
Uref = read_simulation_properties('Uref')


#################### Define Output directory #############################
path_out = os.path.join('.', 'system')


#################### Create forceCoeffs file #############################
txt = '''/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2012                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
'''

for i in range(1, Nb + 1):

    txt += f'''\n
forceCoeffsBlade{i}
{{
    type        forceCoeffs;

    libs ( "libforces.so" );


    writeControl   timeStep;
    writeInterval   1;

    log         yes;

    patches     ( blade{i}Wall);
    rho         rhoInf;      // Indicates incompressible
    rhoInf      {rho};       // Redundant for incompressible
    liftDir     (0 1 0);
    dragDir     (1 0 0);
    CofR        (0 0 0);     // Center of rotation
    pitchAxis   (0 0 1);
    magUInf     {Uref};      // Reference flow speed
    lRef        {c};         // Reference length
    Aref        {Aref};      // Reference area
}}
'''

txt += '''

// ************************************************************************* //'''


path = os.path.join('system', 'forceCoeffs')
print('Save forceCoeffs in path ' + path)

with open(path, 'w') as f:
    f.write(txt)



