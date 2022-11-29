#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file
	'./system/surfaces'
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
Ninterval = read_simulation_properties('writeInterval_surface')



#################### Define Output directory #############################
path_out = os.path.join('.', 'system')


#################### writing ./system/forceCoeffs ########################

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += f'''

surfaces
{{
    type            surfaces;
    libs            (sampling);
    writeControl    timeStep;
    writeInterval   {Ninterval};

    surfaceFormat   foam;
    fields          (p wallShearStress);

    // interpolationScheme cellPoint;  //<- default

    surfaces
    {{

        blade1
        {{
            type            patch;
            patches         ("blade1Wall");
        }}

    }}
}} '''

txt += '''

// ************************************************************************* //'''


path = os.path.join('system', 'surfaces')
print('Save surfaces in path ' + path)

with open(path, 'w') as f:
    f.write(txt)



