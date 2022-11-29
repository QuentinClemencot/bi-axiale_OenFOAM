#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file
	'./system/extrudeMeshDict'
    This template use './simulationProperties.py' file as input.
'''

import os
import sys
sys.path.append(
    "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axial_turbine/py/func")
from preprocess_function_def import read_simulation_properties
from decimal import Decimal

path_OF_header = "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/\
bi_axial_turbine/templates/OF_header.txt" #header for openfoam file

#################### Read usefull parameters #############################
b = read_simulation_properties('b')

b = Decimal(str(b))
#in Decimal format to avoid writing of "0.40000000001" and have "0.4" instead.


#################### Define Output directory #############################
path_out = os.path.join('.', 'system')



#################### Create blockMeshDict file ############################

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      extrudeMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// What to extrude:
//      patch   : from patch of another case ('sourceCase')
//      mesh    : as above but with original case included
//      surface : from externally read surface

constructFrom patch;
sourceCase ".";
sourcePatches (front);

// If construct from patch: patch to use for back (can be same as sourcePatch)
exposedPatchName back;

// Flip surface normals before usage. Valid only for extrude from surface or
// patch.
flipNormals true;

//- Linear extrusion in point-normal direction
extrudeModel        linearNormal;

nLayers             1;

expansionRatio      1.0;

linearNormalCoeffs
{{
    thickness       {b};
}}

// Do front and back need to be merged? Usually only makes sense for 360
// degree wedges.
mergeFaces false;   //true;

// Merge small edges. Fraction of bounding box.
mergeTol 0;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //'''


path = os.path.join(path_out, 'extrudeMeshDict')
print('Save extrudeMeshDict in path ' + path)

with open(path, 'w') as f:
    f.write(txt)



