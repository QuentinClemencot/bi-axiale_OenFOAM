#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file
	'./system/probes'
    This template use './simulationProperties.py' file as input.
'''

import os
import sys
sys.path.append(
    "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axial_turbine/py/func")
from preprocess_function_def import read_simulation_properties
from math import pi
import numpy as np
path_OF_header = "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/\
bi_axial_turbine/templates/OF_header.txt" #header for openfoam file
#################### Read usefull parameters #############################
r = read_simulation_properties('r')
h = read_simulation_properties('h')
c = read_simulation_properties('c')

t = 0.12 # blade tickness
#################### Define Output directory #############################
path_out = os.path.join('.', 'system')


#################### Compute probes location #############################

# probes will be uniformly distributed over the trajectory along 2 lines 
# which follow blade trajectory:
# an 'innner line', at r - 2 * c * t in the turns
# (we want to be as close as possible from blade but not in contact with it)
Rin = r - 2 * c * t

# an 'outer line', at r + 2 * c * t in the turns
Rout = r + 2 * c * t

# compute number all probes
# we want at least 10 probes in eatch turn
# and probes at every transition between turn and translation
l = 2 * pi * r + 2 * h

Nb_turn = 10
Nb_total = Nb_turn * l / (pi * r)
Nb_trans = int(Nb_total * h / l)

# probes in first translation:
x_t1_in = np.linspace(start = -Rin, stop = -Rin, num = Nb_trans, endpoint = False)
x_t1_out = np.linspace(start = -Rout, stop = -Rout, num = Nb_trans, endpoint = False)

y_t1_in = np.linspace(start = 0, stop = -h, num = Nb_trans, endpoint = False)
y_t1_out = np.linspace(start = 0, stop = -h, num = Nb_trans, endpoint = False)


# probes in first turn:
theta = np.linspace(start = 0, stop = pi, num = Nb_turn, endpoint = False)
x_r1_in = - Rin * np.cos(theta)
x_r1_out = - Rout * np.cos(theta)

y_r1_in = -h - Rin * np.sin(theta)
y_r1_out = -h - Rout * np.sin(theta)

# probes in second translation:
x_t2_in = np.linspace(start = Rin, stop = Rin, num = Nb_trans, endpoint = False)
x_t2_out = np.linspace(start = Rout, stop = Rout, num = Nb_trans, endpoint = False)

y_t2_in = np.linspace(start = -h, stop = 0, num = Nb_trans, endpoint = False)
y_t2_out = np.linspace(start = -h, stop = 0, num = Nb_trans, endpoint = False)

# probes in second turn:
theta = np.linspace(start = 0, stop = pi, num = Nb_turn, endpoint = False)
x_r2_in =  Rin * np.cos(theta)
x_r2_out =  Rout * np.cos(theta)

y_r2_in = Rin * np.sin(theta)
y_r2_out = Rout * np.sin(theta)


# concatenate arrays:
x_in = np.concatenate((x_t1_in, x_r1_in, x_t2_in, x_r2_in))
y_in = np.concatenate((y_t1_in, y_r1_in, y_t2_in, y_r2_in))

x_out = np.concatenate((x_t1_out, x_r1_out, x_t2_out, x_r2_out))
y_out = np.concatenate((y_t1_out, y_r1_out, y_t2_out, y_r2_out))

#################### Compute y=cste probes location #############################

mid_y = -h/2 #line at yhe middle of the turbine
d = h + 2 * r
mid_x = np.concatenate((np.linspace(start = -5*d, stop = -d, num = 11, endpoint = False),
                        np.linspace(start = -d, stop = d, num = 14, endpoint = False),
                        np.linspace(start = d, stop = 5*d, num = 11, endpoint = False)))


#####################################################################################
#####################       writing ./system/sampledSet   #########################
#####################################################################################


with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += '''

probesIn
{

    type        probes;
    libs        ("libsampling.so");

    writeControl    timeStep;
    writeInterval   1;

    fields (U k p);
    setFormat   raw;

    probeLocations
    ( 
'''

for x, y in zip(x_in, y_in):
    txt +=f'''   ({x:.6f}   {y:.6f} 0)
'''

txt +='''
);

}

probesOut
{

    type        probes;
    libs        ("libsampling.so");

    writeControl    timeStep;
    writeInterval   1;

    fields (U k p);
    setFormat   raw;

    probeLocations
    ( 
'''

for x, y in zip(x_out, y_out):
    txt +=f'''   ({x:.6f}   {y:.6f} 0)
'''

txt += '''
);

}



probesMidLine
{

    type        probes;
    libs        ("libsampling.so");

    writeControl    timeStep;
    writeInterval   1;

    fields (U k p);
    setFormat   raw;

    probeLocations
    ( 
'''

for x in mid_x:
    txt +=f'''   ({x:.6f}   {mid_y:.6f} 0)
'''

txt += '''
);

}
// ************************************************************************* //
'''


path = os.path.join('system', 'probes')
print('Save probes in path ' + path)

with open(path, 'w') as f:
    f.write(txt)



