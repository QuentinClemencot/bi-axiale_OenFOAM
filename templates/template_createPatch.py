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
path_out = os.path.join('.', 'system')


#################### Create createPatchDict file #########################

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += '''
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      createPatchDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Do a synchronisation of coupled points after creation of any patches.
// Note: this does not work with points that are on multiple coupled patches
//       with transformations (i.e. cyclics).
pointSync false;

// Patches to create.
patches
(
    {   
        // Name of new patch
        name frontAndBack;

        // Dictionary to construct new patch from
        patchInfo
        {   
            type empty;
        }

        // How to construct: either from 'patches' or 'set'
        constructFrom patches;

        // If constructFrom = patches : names of patches. Wildcards allowed.
        patches (
'''
for i in range(1, Nb + 1):
    txt += f'''				blade{i}frontAndBack
 '''

txt += '''
				);

        // If constructFrom = set : name of faceSet
       // set f0;
    }
   
);


// ************************************************************************* //
	'''
path = os.path.join(path_out, 'createPatchDict')
print('Save createPatchDict in path ' + path)

with open(path, 'w') as f:
    f.write(txt)




