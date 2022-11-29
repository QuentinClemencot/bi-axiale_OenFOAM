#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file
	'./system/setFieldsDict'
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
path_out = os.path.join('.', 'system')


#################### Create topoSetDict file #############################

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += '''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      setFieldsDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


defaultFieldValues
(
    volScalarFieldValue zoneID 123
);

regions
(
    // Set cell values
    // (does zerogradient on boundaries)
    cellToCell
    {
        set c0;

        fieldValues
        (
            volScalarFieldValue zoneID 0
        );
    }'''

for i in range(1, Nb + 1):
    txt += f'''

    cellToCell
    {{
        set set_blade{i};

        fieldValues
        (
            volScalarFieldValue zoneID {i}
        );
    }}'''

txt+='''
);

// ************************************************************************* //'''
path = os.path.join(path_out, 'setFieldsDict')
print('Save setFieldsDict in path ' + path)

with open(path, 'w') as f:
    f.write(txt)




