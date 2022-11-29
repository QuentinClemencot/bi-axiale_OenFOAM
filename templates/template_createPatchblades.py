#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file
	'./system/createPatchDict'
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
path_out = os.path.join('../bladeMesh', 'system')


#################### Create createPatchDict file #########################

for i in range(1, Nb+1):
	with open(path_OF_header, 'r', encoding='utf8') as f:
		txt = f.read()

	txt += f'''

FoamFile
{{
	version     2.0;
	format      ascii;
	class       dictionary;
	object      createPatchDict;
}}

pointSync false;

// Patches to create.
patches
(
    {{
        name blade{i}Overset;
        patchInfo
        {{   
            type overset;
        }}
        constructFrom patches;
        patches (oversetBlade);
    }}

    {{ 
        name blade{i}Wall;
        patchInfo
        {{   
            type wall;
        }}
        constructFrom patches;
        patches (blade);
    }}

    {{
        name blade{i}frontAndBack;
        patchInfo
        {{   
            type empty;
        }}
        constructFrom patches;
        patches (frontAndBack);
    }}
	   
);

// *************************************************************************
'''
	path = os.path.join(path_out, f'createPatchDictblade{i}')
	print(f'Save createPatchDictblade{i} in path ' + path)

	with open(path, 'w') as f:
		f.write(txt)




