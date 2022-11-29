#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate the file
	- './constant/globalProcIds' 
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
from fluidfoam import readmesh, readscalar
import numpy as np


path_OF_header = "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/\
bi_axial_turbine/templates/OF_header.txt" #header for openfoam file
##################### Read usefull parameters ##############
(xa, ya) = read_simulation_properties('pointA')
(xb, yb) = read_simulation_properties('pointB')
(xc, yc) = read_simulation_properties('pointC')
(xd, yd) = read_simulation_properties('pointD')

# check conditions on points:
if xa != xc:
	raise ValueError(f"xa:{xa}, should be egal to xc:{xc}")
if xb != xd:
	raise ValueError(f"xb:{xb}, should be egal to xd:{xd}")
if ya != yb:
	raise ValueError(f"ya:{ya}, should be egal to yb:{yb}")
if yc != yd:
	raise ValueError(f"yc:{yc}, should be egal to yc:{yd}")

nzone1 = read_simulation_properties('nzone1')
nzone2 = read_simulation_properties('nzone2')
nzone3 = read_simulation_properties('nzone3')
nzone4 = read_simulation_properties('nzone4')
nzone5 = read_simulation_properties('nzone5')

nproc = read_simulation_properties('nproc') # toal number of proc

meshPlane = read_simulation_properties('meshPlane')

Nb = read_simulation_properties('Nb') # Nb number of overset

tolerance = read_simulation_properties('tolerance') # Tolerance to non-homogeneity

##################### Define Output directory ##############
path_out_globalProcIds = os.path.join('.', 'constant')
path_out_decomposeParDict = os.path.join('.', 'system')

##################### define function ############################

def count_point_in_2Dbox(x, y, xmin, xmax, ymin, ymax):
    """
    return the number of points inside the box delimited by
    xmin, xmax, ymin, ymax).
    
    agr:
    x,y : 1D array of float (same lenght)
    xmin, xmax, ymin, ymax: floats
    
    return:
    N: int, number of points inside box
    """
    
    N = 0
    for (my_x, my_y) in zip(x, y):
        if (xmin <= my_x <= xmax) and (ymin <= my_y <= ymax):
            N += 1
    return N
    
    
def count_point_interval(x, xmin, xmax):
    """
    return the number of points inside the interval delimited by
    xmin, xmax).
    
    agr:
    x: 1D array of float 
    xmin, xmax: floats
    
    return:
    N: int, number of points inside interval
    """
    
    N = 0
    for my_x in x:
        if xmin <= my_x <= xmax:
            N += 1
    return N
    
def create_subset(x, y,  xmin, xmax, ymin, ymax):
    """
        return sub set of points (x, y) that are inside  box
    """
    xin, yin = [], []
    for (my_x, my_y) in zip(x, y):
        if (xmin <= my_x <= xmax) and (ymin <= my_y <= ymax):
            xin.append(my_x)
            yin.append(my_y)
    return (np.array(xin), np.array(yin))

def is_inside(x, y,  box):
    """
        return True if the point (x, y) is in the box
        delimited by box = (xmin, xmax, ymin, ymax), False otherwise
    """
    return ((box[0] <= x <= box[1]) and (box[2] <= y <= box[3]))
    
    
def init_col_lin_position(ncol, nlin, xmin, xmax, ymin, ymax):
    """
    return an initial position of line and column  separation coordonates
    if ncol, nlin = 1, 1: return ([ymin, ymax], [xmin, xmax])
    arg:
    ncol, nlin: int (must be >= 1)
    xmin, xmax, ymin, ymax: floats
    
    return (coord_lin, coord_col) len(coord_lin) = ncol + 1
    """
    
    width_lin = (ymax - ymin) / nlin
    width_col = (xmax - xmin) / ncol
    
    coord_lin, coord_col = [], []
    
    for i in range(nlin + 1):
        coord_lin.append(ymin + i * width_lin)
        
    for j in range(ncol + 1):
        coord_col.append(xmin + j * width_col)
    
    
    return (coord_lin, coord_col)
        
    

        
def compute_col_lin_position(x, y, xmin, xmax, ymin, ymax, ncol, nlin, ite = 10, tol = 20, verbose = False):
    """
    return list of zone_box with a homogene 
    repartition of number of point per zone
    if ncol, nlin = 1, 1: return ([ymin, ymax], [xmin, xmax])
    arg:
    x, y: 1D array of float (same lenght, can be outside box)
    xmin, xmax, ymin, ymax: floats, limit of the box
    ncol, nlin: int, number of columns and lines (lust be >= 1 )
    ite: int, number of iteration to ajust the homogeneity of point repartition
    tol: float, tolerance in %, if max or min number of cell is tol% from the
           average, we consider the repartition homogeneous, so we stop the
           iteration
           
    -----------------
    |      |        | 
    | box2 |        | 
    |------|  box4  |
    |      |        |
    |      |--------|
    | box1 |        |
    |      |  box3  |
    ----------------
    
    return:
     zone_box: list of lenght ncol, nlin 4-uplet [(xmin1, xmax1, ymin1, ymax1), (xmax2,..), ..]
     max_non_homo: float, give the % of max non-homo |max-ave/ave|
    """
    
    # in all this function, i will be indice of line and j indice of column
    zone_box = []
    
    # first, we work only on cell inside box
    xin, yin = create_subset(x, y,  xmin, xmax, ymin, ymax)
    
    # number of zone:
    nzone = ncol * nlin
    
    # average number of cell per zone: 
    ave = len(xin) / nzone
    
    # average number of cell per column: 
    ave_col = len(xin) / ncol
    
    # init  lin and column separation coordonates
    coord_col = np.linspace(xmin, xmax, num = ncol + 1)
    # init number of points per column:
    npoints_col = -1 * np.ones(ncol)
    for j in range(ncol):
            npoints_col[j] = count_point_interval(x = xin,
                                                  xmin = coord_col[j],
                                                  xmax = coord_col[j + 1])
        
    # compute homogeneity of column
    max_non_homo_col = 0
    for j in range(ncol):
        if (npoints_col[j] - ave_col) / ave_col * 100 > max_non_homo_col:
            max_non_homo_col = (npoints_col[j] - ave_col) / ave_col * 100
            
    if verbose:
        print(f"initial col x:{[round(x, 5) for x in coord_col]}")
        for j in range(ncol):
              print(f"initial npoint of column {j}: {npoints_col[j]}")
        print(f"initial max non-homogeneity: {max_non_homo_col:.4f}%")
        print(f"average number of cell per column: {ave_col:.1f}")
    
    # do a loop to decrease the non-homogenity of point distribution by
    # moving internal columns:
    ite_i = 0
    while  max_non_homo_col > tol / 2 and ite_i < ite:
        # for every sub zone, if npoint > ave: we ask the zone to be smaller, 
        # i.e. the zone limit to be reduce.
        # if npoint < ave: we ask the zone to be bigger
        # the changing of coordonate will be the average demand of all
        # zones sharing the column or line separation
        new_col_weight = np.zeros(ncol)
        
       
        for j in range(ncol):
            xmin_j, xmax_j = coord_col[j], coord_col[j + 1]
            dx = xmax_j - xmin_j
            w = 0.5
             #if 0: no change, 1: fast change
            new = dx * ave_col / npoints_col[j]
            new_col_weight[j] = (1 - w) * dx + w * new
                       
        # sum of new_col_weight must be egal to xmin - xmax
        new_col_weight =((xmax - xmin) / np.sum(new_col_weight)) *  new_col_weight
        new_cum_col_weight = np.cumsum(new_col_weight)
        
        # update coord_col:
        for j in range(ncol - 1):
            coord_col[j + 1] = xmin + new_cum_col_weight[j]
            
        # update number of points per column:
        for j in range(ncol):
            npoints_col[j] = count_point_in_2Dbox(x = xin,
                                                  y = yin,
                                                  xmin = coord_col[j],
                                                  xmax = coord_col[j + 1],
                                                  ymin = -10**5,
                                                  ymax = 10**5)
            
        # update homogeneity:
        max_non_homo_col = 0
        for j in range(ncol):
            if (npoints_col[j] - ave_col) / ave_col * 100 > max_non_homo_col:
                max_non_homo_col = (npoints_col[j] - ave_col) / ave_col * 100
                
        # increment iteration counter:
        ite_i += 1
        if verbose:
            print(f"iteration {ite_i}:")
            print(f"col x:{[round(x, 5) for x in coord_col]}")
            print("npoint:")
            print(f"{[npoints_col[j] for j in range(ncol)]}")
            print(f"max non-homogeneity: {max_non_homo_col:.4f}%")
            
            
    # now we do the same or every colomn:
    for j in range(ncol):
        xmin_j, xmax_j = coord_col[j], coord_col[j + 1]
        
        # first, we work only on cell inside this column
        xin_j, yin_j = create_subset(x, y,  xmin = xmin_j, xmax = xmax_j, ymin = ymin, ymax = ymax)
        
        # average number of cell per line: 
        ave_lin = len(xin_j) / nlin
        
        # init  lin separation coordonates
        coord_lin = np.linspace(ymin, ymax, num = nlin + 1)
        
        # init number of points per line:
        npoints_lin = -1 * np.ones(nlin)
        for i in range(nlin):
                npoints_lin[i] = count_point_interval(x = yin_j,
                                                      xmin = coord_lin[i],
                                                      xmax = coord_lin[i + 1])
            
        # compute homogeneity of line
        max_non_homo_lin = 0
        for i in range(nlin):
            if (npoints_lin[i] - ave_lin) / ave_lin * 100 > max_non_homo_lin:
                max_non_homo_lin = (npoints_lin[i] - ave_lin) / ave_lin * 100
                
        if verbose:
            print(f"For column {j}:")
            print(f"initial lin y:{[round(y, 5) for y in coord_lin]}")
            for i in range(nlin):
                  print(f"initial npoint of lin {i}: {npoints_lin[i]}")
            print(f"initial max non-homogeneity: {max_non_homo_lin:.4f}%")
            print(f"average number of cell per line: {ave_lin:.1f}")
        
        # do a loop to decrease the non-homogenity of point distribution by
        # moving internal columns:
        ite_i = 0
        while  max_non_homo_lin > tol / 2 and ite_i < ite:
            # for every sub zone, if npoint > ave: we ask the zone to be smaller, 
            # i.e. the zone limit to be reduce.
            # if npoint < ave: we ask the zone to be bigger
            new_lin_height = np.zeros(nlin)
            
           
            for i in range(nlin):
                ymin_i, ymax_i = coord_lin[i], coord_lin[i + 1]
                dy = ymax_i - ymin_i
                w = 0.1
                 #if 0: no change, 1: fast change
                new = dy * ave_lin / npoints_lin[i]
                new_lin_height[i] = (1 - w) * dy + w * new
                           
            # sum of new_col_weight must be egal to xmin - xmax
            #print(f"new_col_weight before sizing: {new_col_weight}")
            new_lin_height =((ymax - ymin) / np.sum(new_lin_height)) *  new_lin_height
            #print(f"new_col_weight after sizing: {new_col_weight}")
            new_cum_lin_height = np.cumsum(new_lin_height)
            #print(f"new_col_cum: {new_cum_col_weight}")
            
            # update coord_col:
            for i in range(nlin - 1):
                coord_lin[i + 1] = ymin + new_cum_lin_height[i]
                
            # update number of points per column:
            for i in range(nlin):
                npoints_lin[i] = count_point_interval(x = yin_j,
                                                      xmin = coord_lin[i],
                                                      xmax = coord_lin[i + 1])
                
            # update homogeneity:
            max_non_homo_lin = 0
            for i in range(nlin):
                if (npoints_lin[i] - ave_lin) / ave_lin * 100 > max_non_homo_lin:
                    max_non_homo_lin = (npoints_lin[i] - ave_lin) / ave_lin * 100
                    
            # increment iteration counter:
            ite_i += 1
            if verbose:
                print(f"For column {j}:")
                print(f"iteration {ite_i}:")
                print(f"lin y:{[round(y, 5) for y in coord_lin]}")
                print("npoint:")
                print(f"{[npoints_lin[i] for i in range(nlin)]}")
                print(f"max non-homogeneity: {max_non_homo_lin:.4f}%")
        
        for i in range(nlin):
            zone_box.append((xmin_j, xmax_j, coord_lin[i], coord_lin[i + 1]))  
            
    if verbose:
        print(f"Zone size: xmin:{xmin:.4f}, xmax:{xmax:.4f}, ymin:{ymin:.4f}, ymax:{ymax:.4f}")
        for box in zone_box:
            print(f"add box xmin:{box[0]:.4f}, xmax:{box[1]:.4f}, ymin:{box[2]:.4f}, ymax:{box[3]:.4f}")
    return (zone_box)
##################### Read mesh ############################
sol = '.'

x_all, y_all, z_all = readmesh(sol)



xmin, xmax = np.min(x_all), np.max(x_all)
ymin, ymax = np.min(y_all), np.max(y_all)
zmin, zmax = np.min(z_all), np.max(z_all)
# not the true limit but the limit of cell center

### we assume the mesh is 2D plane, so only to coordonates are usefull:
if meshPlane == 'xy':
    x0_all = x_all
    x1_all = y_all
    x0min, x0max = xmin, xmax
    x1min, x1max = ymin, ymax
elif meshPlane == 'xz':
    x0_all = x_all
    x1_all = z_all
    x0min, x0max = xmin, xmax
    x1min, x1max = zmin, zmax
elif meshPlane == 'yz':
    x0_all = y_all
    x1_all = z_all
    x0min, x0max = ymin, ymax
    x1min, x1max = zmin, zmax
else:
    raise ValueError(f"meshPlane: {meshPlane}, should be either egal to \
'xy', 'xz' or 'yz'!")

### limit of every zone:
limit_zone = [((x0min, x0max), (ya, x1max)), 
            ((x0min, xa),    (yc, ya)),
            ((xa, xb),       (yc, ya)),
            ((xb, x0max),    (yc, ya)),
            ((x0min, x0max), (x1min, yc))]


### number of lines and colum for eatch zone
ncol_nlin_per_zone = [nzone1,
                        nzone2,
                        nzone3,
                        nzone4,
                        nzone5]
				  
### number of proc per zone
nproc_per_zone = [np.prod(nzone1),
                np.prod(nzone2),
                np.prod(nzone3),
                np.prod(nzone4),
                np.prod(nzone5)]


cumul_nproc = np.cumsum(nproc_per_zone)

computed_nproc = cumul_nproc[-1]

if computed_nproc != nproc:
    raise ValueError(f"total number of proc computed: {computed_nproc} \
is not egal to specified number of proc in \
simulationProperties.py file: {nproc}, must be!")

# Compute the line and column coordonate  of zone 3 to ensure a certain level of
# homogeneity in cell proc distribution:
zone3_procBox = compute_col_lin_position(x = x0_all,
                                         y = x1_all,
                                         xmin = xa,
                                         xmax = xb,
                                         ymin = yc,
                                         ymax = ya,
                                         ncol = nzone3[0],
                                         nlin = nzone3[1],
                                         ite = 100,
                                         tol = tolerance,
                                         verbose = False)


globalProcIds = -1 * np.ones(len(x_all)) # will contain procId for every cells

for i, (x, y) in enumerate(zip(x0_all, x1_all)):
    # first; determine meshPart:
    if y > ya:
        meshPart = 1
    elif y < yc:
        meshPart = 5
    elif yc <= y <= ya and x < xa:
        meshPart = 2
    elif yc <= y <= ya and xa <= x <= xb:
        meshPart = 3
    elif yc <= y <= ya and x > xb:
        meshPart = 4
    else:
        raise ValueError(f"The point ({x}, {y}) do not belong to any zone !")

    if meshPart == 3:
        my_gloabl_proc = -1
        for j, my_box in enumerate(zone3_procBox):
            if is_inside(x, y, my_box):
                my_gloabl_proc = j + cumul_nproc[1]
                continue
        if my_gloabl_proc == -1:
            print(f"List of proc box of zone 3:")
            for my_box in zone3_procBox:
                print(my_box)
            raise ValueError(f"The point ({x}, {y}) do not belong to any proc of zone 3 !")
    
    else:
        # number of 'line' and 'column' of proc for this zone
        n_col, n_lin = ncol_nlin_per_zone[meshPart - 1] #meshPart start from 1, not 0.

        # find the 'line' and 'column' local proc the cell belong to:
        my_zone_limit = limit_zone[meshPart - 1]

        x_len = my_zone_limit[0][1] - my_zone_limit[0][0]
        y_len = my_zone_limit[1][1] - my_zone_limit[1][0]
        dx = x - my_zone_limit[0][0] # x coordonate of cell relatively to the zone
						             # bottom left coorner
        dy = y - my_zone_limit[1][0]

        # find the local proc number (in [0, nproc_per_zone[meshPart - 1] -1)
        my_col = int(dx / x_len *  n_col) # goes from 0 to n_col - 1
        my_lin = int(dy / y_len *  n_lin) # goes from 0 to n_lin - 1
        # bricolage: sometime my_col == n_col: not good (prcision issues):
        if my_col == n_col:
            my_col = n_col -1
        if my_lin == n_lin:
            my_lin = n_lin -1
            
        my_local_proc = my_col + my_lin * n_col

        # find global proc:
        if meshPart == 1:
            my_gloabl_proc = my_local_proc
        else:
            my_gloabl_proc = my_local_proc + cumul_nproc[meshPart - 2]

        if my_gloabl_proc == -1:
            raise ValueError(f"The point:({x}, {y}) of mesh part:{meshPart} have a  my_gloabl_proc: {my_gloabl_proc}!")

    globalProcIds[i] = int(my_gloabl_proc)
	
print(f"Write globalProcIds, procIds goes from {np.min(globalProcIds)} to \
{np.max(globalProcIds)}.")
									 
##################### Create globalProcIds ###############################    

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()
	
txt += """
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       labelList;
    location    "constant";
    object      globalProcIds;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


"""

txt += f"""

{len(globalProcIds)}
(
"""  

for proci in globalProcIds:
    txt += f"{proci}\n" 
    
txt += """)


// ************************************************************************* //"""

               
##################### Create decomposeParDict ############################### 
with open(path_OF_header, 'r', encoding='utf8') as f:
	txt2 = f.read()
	

txt2 += """
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      decomposeParDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
"""

txt2 += f"""
numberOfSubdomains {nproc}; // your total number of processors
"""

txt2 += """
method manual;
manualCoeffs
{
    dataFile    "globalProcIds";
}

// ************************************************************************* //"""

##################### Write Files ###############################                       
path = os.path.join(path_out_globalProcIds, 'globalProcIds')
print('Save globalProcIds in path ' + path)
with open(path, 'w') as f:
    f.write(txt)
    
path2 = os.path.join(path_out_decomposeParDict, 'decomposeParDict')
print('Save decomposeParDict in path ' + path2)
with open(path2, 'w') as f:
    f.write(txt2)   
    
    
