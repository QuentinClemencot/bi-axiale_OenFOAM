#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file
	'./system/topoSetDict'
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
    object      topoSetDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

actions
(
    {
        name    c0;
        type    cellSet;
        action  new;
        source  regionToCell;
        insidePoints ((-1.99 0 0));
    }'''


for i in range(1, Nb+1):

	txt += f'''

	{{
		name    blade{i}FaceSet;
		type    faceSet;
		action  new;
		source  patchToFace;
		sourceInfo
		{{
		  patch blade{i}frontAndBack;
		}}
	}}



    {{
        name    set_blade{i};
        type    cellSet;
        action  new;
        source  faceToCell;
        sourceInfo
        {{
            set blade{i}FaceSet;             // Name of faceSet

       //option neighbour; // cell with neighbour in faceSet
        //option owner;     //  ,,       owner
            option any;         // cell with any face in faceSet
            //option all;       // cell with all faces in faceSet
        }} 
    }}


	{{
		name    blade{i}FaceSet;
		type    faceSet;
		action  remove;
	}}'''

txt += '''

);

// ************************************************************************* //'''

path = os.path.join(path_out, 'topoSetDict')
print('Save topoSetDict in path ' + path)

with open(path, 'w') as f:
	f.write(txt)




