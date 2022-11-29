#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file
	'./system/blockMeshDict'
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
c = read_simulation_properties('c')
cellmin = read_simulation_properties('cellmin')
h = read_simulation_properties('h')
r = read_simulation_properties('r')

b = Decimal(str(b))
d = Decimal(str(h)) + Decimal('2')*Decimal(str(r))     
#width of the area covered by the center of rotation of the foil
#in Decimal format to avoid writing of "0.40000000001" and have "0.4" instead.


#################### Define Output directory #############################
path_out = os.path.join('.', 'system')


#################### Calcul of background domain size ###################

H = Decimal('60')*d #weight of the domain
ymax = Decimal('30')*d  # y-limit of the domain
ymin = Decimal('-30')*d

L = Decimal('55')*d  #lenght of the domain
xmin = Decimal('-15')*d # x-limit of the domain
xmax = Decimal('40')*d

zmin = Decimal('-0.5')*b # x-limit of the domain
zmax = Decimal('0.5')*b

#parameters to be change here :
N = 9 #number of raffinement box for the background mesh 


cellmax=cellmin*(2**N) #size of biggest cell of background 

Nx=float(L)//cellmax #Number of cell along x direction 
Nx=int(Nx)
Ny=float(H)//cellmax #Number of cell along y direction
Ny=int(Ny)


#################### Create blockMeshDict file ############################

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      blockMeshDict;
}}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

scale   1;

vertices
(
    ({xmin} {ymin} {zmin})
    ({xmax} {ymin} {zmin})
    ({xmax} {ymax} {zmin})
    ({xmin} {ymax} {zmin})
    ({xmin} {ymin} {zmax})
    ({xmax} {ymin} {zmax})
    ({xmax} {ymax} {zmax})
    ({xmin} {ymax} {zmax})
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({Nx} {Ny} 1) simpleGrading (1 1 1)
);
'''

txt += '''


edges
(
);

boundary
(

   //this is an empty patch to ensure that the overset patch will be the first 
   //patch on final constant/polyMesh/boundary file
    oversetAirfoil_1 
    {
        type overset;
        faces ();
    }

    
    top
    {
        type symmetryPlane;
        faces
        (
            (2 3 7 6)
        );
    }

    bottom
    {
        type symmetryPlane;
        faces
        (
            (0 1 5 4)
        );
    }

    inlet
    {
        type patch;
        faces
        (
            (0 4 7 3)
        );
    }

    outlet
    {
        type patch;
        faces
        (
            (1 2 6 5)
        );
    }

    front
    {
        type empty;
        faces
        (
            (0 3 2 1)
        );
    }
    
    back
    {
        type empty;
        faces
        (
            (4 5 6 7)
        );
    }
);


// ************************************************************************* //'''


path = os.path.join(path_out, 'blockMeshDict')
print('Save blockMeshDict in path ' + path)

with open(path, 'w') as f:
    f.write(txt)



