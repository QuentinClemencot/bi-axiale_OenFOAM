#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file
	'./0_orig/pointDisplacement'
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


#################### Create pointDisplacement file #######################
with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()


txt += '''
FoamFile
{
    version     2.0;
    format      ascii;
    class       pointVectorField;
    object      pointDisplacement;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dimensions      [0 1 0 0 0 0 0];

internalField   uniform (0 0 0);

boundaryField
{

'''

for i in range(1, Nb+1):

	txt += f'''\n
	blade{i}Overset
    {{
        patchType       overset;
        type            zeroGradient;
    }}'''


for i in range(1, Nb+1):

	txt += f'''\n
	blade{i}Wall
    {{
        type            uniformFixedValue;
        uniformValue    (0 0 0);
    }}'''

txt += '''

    inlet
    {
        type            uniformFixedValue;
        uniformValue    (0 0 0);
    }

    outlet
    {
        type            uniformFixedValue;
        uniformValue    (0 0 0);
    }

    top
    {
        type            uniformFixedValue;
        uniformValue    (0 0 0);
    }

    bottom
    {
        type            uniformFixedValue;
        uniformValue    (0 0 0);
    }

	frontAndBack
    {
        type            empty;
    }
}

// ************************************************************************* //'''
path = os.path.join(path_out, 'pointDisplacement')
print('Save pointDisplacement in path ' + path)

with open(path, 'w') as f:
    f.write(txt)




