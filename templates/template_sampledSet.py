#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file
	'./system/sampledSet'
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
r = read_simulation_properties('r')
h = read_simulation_properties('h')

r = Decimal(str(r))
h = Decimal(str(h))
#################### Define Output directory #############################
path_out = os.path.join('.', 'system')



#####################################################################################
#####################       writing ./system/sampledSet   #########################
#####################################################################################


with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += f'''

sampledSet
{{
    type        sets;

    libs ( "libsampling.so" );


    writeControl   timeStep;
    writeInterval   10;

    interpolationScheme   cellPoint;

	setFormat   raw;

    fields (U p);

    sets
    (
        x_minus5r
             type    uniform;
             axis    xyz;
             start   ( {Decimal('-5')*r} {-h-r} 0.0);
             end     ( {Decimal('-5')*r} {r} 0.0);
             nPoints 60;'
         }}

        x_minus2r
        {{
             type    uniform;
             axis    xyz;
             start   ( {Decimal('-2')*r} {-h-r} 0.0);
             end     ( {Decimal('-2')*r} {r} 0.0);
             nPoints 60;
         }}

        x_minus1p2r
        {{
             type    uniform;
             axis    xyz;
             start   ( {Decimal('-1.2')*r} {-h-r} 0.0);
             end     ( {Decimal('-1.2')*r} {r} 0.0);
             nPoints 60;
         }}

        x_minus0p8r
        {{
             type    uniform;
             axis    xyz;
             start   ( {Decimal('-0.8')*r} {-h} 00.0);
             end     ( {Decimal('-0.8')*r} 0 0.0);
             nPoints 40;
        }}

        x_0
        {{
             type    uniform;
             axis    xyz;
             start   ( 0 {-h} 0.0);
             end     ( 0 0 0.0);
             nPoints 40;
        }}

        x_0p8r
        {{
             type    uniform;
             axis    xyz;
             start   ( {Decimal('0.8')*r} {-h} 0.0);
             end     ( {Decimal('0.8')*r} 0 0.0);
             nPoints 40;
        }}

        x_1p2r
        {{
             type    uniform;
             axis    xyz;
             start   ( {Decimal('1.2')*r} {-h-r} 0.0);
             end     ( {Decimal('1.2')*r} {r} 0.0);
             nPoints 60;
         }}

         x_2r
        {{
             type    uniform;
             axis    xyz;
             start   ( {Decimal('2')*r} {-h-r} 0.0);
             end     ( {Decimal('2')*r} {r} 0.0);
             nPoints 60;
        }}

        x_5r
        {{
             type    uniform;
             axis    xyz;
             start   ( {Decimal('5')*r} {-h-r} 0.0);
             end     ( {Decimal('5')*r} {r} 0.0);
             nPoints 60;
         }}
    );
}}


// ************************************************************************* //
'''


path = os.path.join('system', 'sampledSet')
print('Save sampledSet in path ' + path)

with open(path, 'w') as f:
    f.write(txt)



