#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the files  
	- '../bladeMesh/Allrun.pre' 
	- '../bladeMesh/Allrun.merge'
	- '../bladeMesh/system/createPatchDictblade*'  
    This template use './simulationProperties.py' file as input.
'''

import os
import sys
sys.path.append(
    "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axial_turbine/py/func")
from preprocess_function_def import read_simulation_properties
from general_function_def import read_data

from decimal import Decimal
from math import cos, sin

path_py = "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axial_turbine/py/func"

path_OF_header = "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/\
bi_axial_turbine/templates/OF_header.txt" #header for openfoam file

##################### Read usefull parameters ##############
b = read_simulation_properties('b')
c = read_simulation_properties('c')
Nb = read_simulation_properties('Nb')
xa = read_simulation_properties('xa')

b = Decimal(str(b)) #convert in decimal to avoid for example 
                    #to write "0.40000000005" and write "0.4" instead
xa = Decimal(str(xa))
c = Decimal(str(c))


##################### Define Output directory ##############
path_out = os.path.join('..', 'bladeMesh')
counter = 0
for file in os.listdir(path_out):
    if file.endswith(".msh"):
        counter += 1
        mesh_file = file

print('Mesh file for blade overset :' + mesh_file)

if counter == 0:
    raise ValueError('no .msh file in ' + path_out)

if counter > 1:
    raise ValueError('more than one .msh file in ' + path_out)



##################### Read initial blades position ##############
data = read_data('./constant','blade_init_pos.dat')

for key in data :
    if Nb != len(data[key]):
        raise ValueError('Number of blade Nb not egal to informations'
                          + 'in ./constant/blade_init_pos.dat')


#################### Create createPatchDict file ################

for i in range(1, Nb + 1):
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

// ************************************************************************* //'''
    path = os.path.join('../bladeMesh/system', f'createPatchDictblade{i}')
    print(f'Save createPatchDictblade{i} in path {path}')

    with open(path, 'w') as f:
        f.write(txt)


#################### Create Allrun.pre file #####################

txt = f'''
#!/bin/sh
cd ${{0%/*}} || exit 1                        # Run from this directory
. $WM_PROJECT_DIR/bin/tools/RunFunctions    # Tutorial run functions

# Remove old log.* files 
rm -f log.* 
rm -rf constant/polyMesh

# Generate mesh from fluent 3d mesh 
runApplication fluentMeshToFoam -2D {b} {mesh_file}

# Make front and back type empty
cp ./system/createPatchDict0 ./system/createPatchDict
runApplication createPatch -overwrite
mv log.createPatch log.createPatch0

# transform the mesh to be at corect initiale angle and correct chord lenght (default c=1m)

# This translation places point A at the origin of the coordinate system
runApplication transformPoints -translate "(-{xa} 0 0)"  
mv log.transformPoints log.transformPoints0translate

# A rotation of pi is necessary to have theta = 0rad.
# To avoid bugs, we perform two successive rotations of pi/2
runApplication transformPoints -rotate "((1 0 0)(0 1 0))"
mv log.transformPoints log.transformPoints0rotate1
runApplication transformPoints -rotate "((1 0 0)(0 1 0))"
mv log.transformPoints log.transformPoints0rotate2

runApplication transformPoints -scale "({c} {c} 1)"
mv log.transformPoints log.transformPoints0scale
'''

for i in range(1, Nb + 1):

    dirName = f'blade{i}Mesh'
    bladeName = f'blade{i}'
    x0 = data['x'][i-1]
    y0 = data['y'][i-1]
    theta0 = data['theta'][i-1]
    #print('for blade'+str(i+1)+' :')
    #print('theta0='+str(theta0))
    Xrot = round(cos(theta0), 10)
    Yrot = round(sin(theta0), 10)
	
    txt += f'''\n
# Create {bladeName} overset mesh:
rm -rf ../{dirName}
mkdir ../{dirName}
cp -r constant ../{dirName}/
cp -r system ../{dirName}/

runApplication transformPoints -case ../{dirName}/ -rotate "((1 0 0)({Xrot} {Yrot} 0))" 
mv log.transformPoints log.transformPointsRotate{bladeName}

runApplication transformPoints -case ../{dirName}/ -translate "({x0} {y0} 0)" 
mv log.transformPoints log.transformPointsTranslate{bladeName}

cp ./system/createPatchDict{bladeName} ../{dirName}/system/createPatchDict
runApplication createPatch -case ../{dirName} -overwrite
mv log.createPatch	log.createPatch{bladeName}
'''

path = os.path.join('../bladeMesh', 'Allrun.pre')
print('Save Allrun.pre in path ' + path)

with open(path, 'w') as f:
    f.write(txt)


#################### Create Allrun.merge file ###########################

txt = '''

#!/bin/sh
cd ${0%/*} || exit 1                        # Run from this directory
. $WM_PROJECT_DIR/bin/tools/RunFunctions    # Tutorial run functions


# merge overset blade together :
cp -rf  ../blade1Mesh/constant/polyMesh/ ./constant/

'''

for i in range(2, Nb + 1):
    txt += f'''

python {path_py}/mergeMesh.py . ../blade{i}Mesh/
'''

txt += f'''

python {path_py}/mergeMesh.py . ../bladeAndBackground/

'''

txt += '''

#------------------------------------------------------------------------------'''

path = os.path.join(path_out, 'Allrun.merge')
print('Save Allrun.merge in path ' + path)

with open(path, 'w') as f:
    f.write(txt)




