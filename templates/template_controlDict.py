#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the files:
    -'./system/controlDict_regular'
    -'./system/controlDict_start'
    -'./system/controlDict_probes'

    This template use './simulationProperties.py' file as input.
'''


import sys
sys.path.append("/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/bi_axe_multi_blade/py/func")


from general_function_def import read_simulation_properties
import os
from math import pi
from runpy import run_path

path_OF_header = "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/\
bi_axial_turbine/templates/OF_header.txt" #header for openfoam file

#################### Read usefull parameters #############################
deltaT_probes = read_simulation_properties('deltaT_probes')
deltaT_regular = read_simulation_properties('deltaT_regular')
deltaT_start = read_simulation_properties('deltaT_start')

Nb = read_simulation_properties('Nb')

nRevRegular = read_simulation_properties('nRevRegular')
nRevProbes = read_simulation_properties('nRevProbes')

nStep_per_revProbes = read_simulation_properties('nStep_per_revProbes')
nStep_per_revRegular = read_simulation_properties('nStep_per_revRegular')

nStepStart = read_simulation_properties('nStepStart')

T = read_simulation_properties('T')
writeInterval_timeDir_reg = read_simulation_properties('writeInterval_timeDir_reg')
timeDir_perCycle_probes = read_simulation_properties('timeDir_perCycle_probes')

#################### Define Output directory #############################
path_out = os.path.join('.', 'system')


#################### Calculation of endTime ##############################

endTimeReg = T * (nRevRegular)
endTimeReg = round(endTimeReg, 6)
endTimeProb = T * (nRevRegular + nRevProbes)
endTimeProb = round(endTimeProb, 6)

#################### writeInterval_timeDir_probes ###########################
writeInterval_timeDir_probes = nStep_per_revProbes // timeDir_perCycle_probes


# 3 differents controlDict file are use : 

#  controlDict_start: for the first iterations
#  controlDict_regular: it is used until the permanent state is reached
#  controlDict_probes: the final one, used for surfaces and lines probes on few periods


#################### Writing controlDict_start ###########################

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


libs            (overset fvMotionSolvers myoverset myfieldFunctionObjects);


application     overPimpleDyMFoam;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         {deltaT_start*(nStepStart+1)};

deltaT          {deltaT_start};

writeControl    timeStep;

writeInterval   {nStepStart};

purgeWrite      2;

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   7;

runTimeModifiable true;


functions
{{
	#include	"forces"
	#includeFunc  solverInfo
	
    yPlusExtremum
    {{
        type            yPlusExtremum;
        libs            (myfieldFunctionObjects);
        patches         (blade1Wall);
        writeControl    timeStep;
        executeControl  timeStep;
        writeInterval   5;
        executeInterval 5;
    }}

	
}}

// ************************************************************************* //
'''



path = os.path.join('./system', 'controlDict_start')
print('Save controlDict_start in path ' + path)

with open(path, 'w') as f:
    f.write(txt)

path = os.path.join('./system', 'controlDict')
print('Save controlDict in path ' + path)

with open(path, 'w') as f:
    f.write(txt)
    
#################### Writing controlDict_regular #########################

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


libs            (overset fvMotionSolvers myoverset myfieldFunctionObjects);


application     overPimpleDyMFoam;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         {endTimeReg};

deltaT          {deltaT_regular}; //{nStep_per_revRegular} step per revolution

writeControl    timeStep;

writeInterval   {writeInterval_timeDir_reg};

purgeWrite      2;

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   11;

runTimeModifiable true;



functions
{{
	#include	"forces"
	#includeFunc  solverInfo


    yPlus
    {{
        type            yPlus;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        patches         (blade1Wall);
    }}

    yPlusExtremum
    {{
        type            yPlusExtremum;
        libs            (myfieldFunctionObjects);
        patches         (blade1Wall);
        writeControl    timeStep;
        executeControl  timeStep;
        writeInterval   5;
        executeInterval 5;
    }}
	
}}


// ************************************************************************* //'''


path = os.path.join('./system', 'controlDict_regular')
print('Save controlDict_regular in path ' + path)

with open(path, 'w') as f:
    f.write(txt)



#################### Writing controlDict_probes ##########################

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    location    "system";
    object      controlDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

libs            (overset fvMotionSolvers myoverset myfieldFunctionObjects);


application     overPimpleDyMFoam;

startFrom       latestTime;

startTime       0;

stopAt          endTime;

endTime         {endTimeProb};

deltaT          {deltaT_probes}; // {nStep_per_revProbes} step per revolution

writeControl    timeStep;

writeInterval   {writeInterval_timeDir_probes};

purgeWrite      {timeDir_perCycle_probes};

writeFormat     ascii;

writePrecision  8;

writeCompression off;

timeFormat      general;

timePrecision   11;

runTimeModifiable true;


functions
{{
    #include	"forces"
    #include	"surfaces"
    #include	"probes"
    #includeFunc  solverInfo


    vorticity1
    {{
        // Mandatory entries (unmodifiable)
        type        vorticity;
        libs        (fieldFunctionObjects);

        // Optional (inherited) entries
        executeControl  timeStep;
        executeInterval {writeInterval_timeDir_probes};
        writeControl    timeStep;
        writeInterval   {writeInterval_timeDir_probes};
    }}

    yPlus
    {{
        type            yPlus;
        libs            (fieldFunctionObjects);
        writeControl    writeTime;
        patches         (blade1Wall);
    }}
    
    yPlusExtremum
    {{
        type            yPlusExtremum;
        libs            (myfieldFunctionObjects);
        patches         (blade1Wall);
        writeControl    timeStep;
        executeControl  timeStep;
        writeInterval   1;
        executeInterval 1;
    }}

    wallShearStressBlades
    {{
        type        wallShearStress;
        libs        ("libfieldFunctionObjects.so");
        writeControl    writeTime;
        executeControl  timeStep;
        executeInterval 1;
        patches     (
'''
for i in range(1, Nb+1):
	txt += f'''
                     blade{i}Wall
           '''  


txt += '''
                    );
    }

	
}


// ************************************************************************* //'''

path = os.path.join('./system', 'controlDict_probes')
print('Save controlDict_probes in path ' + path)

with open(path, 'w') as f:
    f.write(txt)


