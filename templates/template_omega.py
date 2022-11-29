#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file
	'./0_orig/omega'
    This template use './simulationProperties.py' file as input.
'''

import os
import sys
sys.path.append(
    "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axial_turbine/py/func")
from preprocess_function_def import read_simulation_properties

path_OF_header = "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/\
bi_axial_turbine/templates/OF_header.txt" #header for openfoam file

#################### Read usefull parameters #############################
Nb = read_simulation_properties('Nb')

#################### Define Output directory #############################
path_out = os.path.join('.', '0_orig')


#################### Create omega file ###################################

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()


txt += '''
FoamFile
{
    version     2.0;
    format      ascii;
    class       volScalarField;
    object      omega;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include        "include/initialConditions"

dimensions      [ 0 0 -1 0 0 0 0 ];

internalField   uniform $turbulentOmega;

boundaryField
{
    #includeEtc "caseDicts/setConstraintTypes"

'''

for i in range(1, Nb + 1):

    txt += f'''\n
    blade{i}Overset
    {{
        type            overset;
        value           $internalField;
    }}'''


for i in range(1, Nb + 1):

    txt += f'''\n
    blade{i}Wall
    {{
        type            omegaWallFunction;
        value           $internalField;
    }}'''

txt += '''

    inlet
    {
        type            fixedValue;
        value           $internalField;
    }

    outlet
    {
        type            inletOutlet;
        inletValue      $internalField;
        value           $internalField;
    }

    top
    {
        type            symmetryPlane;
    }

    bottom
    {
        type            symmetryPlane;
    }

	frontAndBack
    {
        type            empty;
    }
}

// ************************************************************************* //'''
path = os.path.join(path_out, 'omega')
print('Save omega in path ' + path)

with open(path, 'w') as f:
    f.write(txt)
