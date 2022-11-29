#!/usr/bin/env python
# coding: utf-8

# import usefull libs:
import sys
from fluidfoam import readfield, OpenFoamFile
import os 
import copy
import numpy as np


# path to OpenFOAM file header:
path_OF_header = "/fsnet/project/nrj/2019/19TRIBINE/simulation/reference_cases/\
bi_axial_turbine/templates/OF_header.txt" 



def readMesh(path, d):
    '''
    add to dictionnary d keys 'faces', 'points', 'owner', 'neighbour', 'boundary'
    '''
    print(f"path: {path}")
    d_up = copy.deepcopy(d)
    d_up['faces'] = OpenFoamFile(path = path, name = 'faces')
    d_up['points'] = OpenFoamFile(path = path, name = 'points')
    d_up['owner'] = OpenFoamFile(path = path, name = 'owner')
    d_up['neighbour'] = OpenFoamFile(path = path, name = 'neighbour')
    d_up['boundary'] = OpenFoamFile(path = path, name = 'boundary')

    return d_up
    
def mergeMesh(master, slave):
    '''
    will merge the 2 meshes
    in the new boundary file, boundaries of master will appear first
    '''

    
    # start by concatenate point of master and slaves
    new_points = np.concatenate((master['points'].values, slave['points'].values))
    new_points_x = new_points[::3]
    new_points_y = new_points[1::3]
    new_points_z = new_points[2::3]
        
    # then correct the slave faces indices (add + nb_master_points to every points id)
    nb_master_points = len(master['points'].values_x)
    tmp_slave_faces = {}
   
    for i in range(slave['faces'].nfaces):
        # tmp_face is a dict with keys 'nptns' 'id_pts'
        tmp_face = slave['faces'].faces[i]
        tmp_face['id_pts'] = np.array([a + nb_master_points for a in tmp_face['id_pts']])
        tmp_slave_faces[i] = tmp_face

    # reorder faces to have:
    # -first: internal master faces 
    # -then: internal slave faces
    # -thend: boudary master faces
    # -finally: boundary slave faces
    nb_int_face_slave  = len(slave['neighbour'].values)
    nb_int_face_master = len(master['neighbour'].values)


    new_faces = {}
    new_owner = np.zeros(master['owner'].nb_faces + slave['owner'].nb_faces)
    
    # add internal master faces
    for i in range(nb_int_face_master):
        # i in [0, nb_int_face_master[
        new_faces[i] = master['faces'].faces[i]['id_pts']
        new_owner[i] = master['owner'].values[i]
        
    # add internal slave faces   
    for i in range(nb_int_face_slave):
        # i in [0, nb_int_face_slave[
        new_faces[i + nb_int_face_master] = slave['faces'].faces[i]['id_pts']
        new_owner[i + nb_int_face_master] = slave['owner'].values[i] + master['owner'].nb_cell
        
    # add boundary master faces   
    for i in range(nb_int_face_master, master['faces'].nfaces):
        # i in [nb_int_face_master, nb_face_master[
        new_faces[i + nb_int_face_slave] = master['faces'].faces[i]['id_pts']
        new_owner[i + nb_int_face_slave] = master['owner'].values[i] 
        
    # add boundary slave faces   
    for i in range(nb_int_face_slave, slave['faces'].nfaces):
        # i in [nb_int_face_slave, nb_face_slave[
        new_faces[i + master['faces'].nfaces] = slave['faces'].faces[i]['id_pts']
        new_owner[i + master['faces'].nfaces] = slave['owner'].values[i] + master['owner'].nb_cell
        
        
    # merge the neighbour
    new_neighbour = np.concatenate((master['neighbour'].values, slave['neighbour'].values + master['owner'].nb_cell))
    new_neighbour = np.concatenate((master['neighbour'].values, slave['neighbour'].values + master['owner'].nb_cell))
    
    # merge the boundary:
    new_boundary = copy.deepcopy(master["boundary"].boundaryface)
    # change indice of faces (+nb_int_face_slave)
    for key in new_boundary:
        tmp_startFace = int(new_boundary[key][b'startFace'].decode('UTF-8')) + nb_int_face_slave
        tmp_startFace = str(tmp_startFace).encode('UTF-8')
        new_boundary[key][b'startFace'] = tmp_startFace
    for key in slave['boundary'].boundaryface:
        new_boundary[key] = copy.deepcopy(slave['boundary'].boundaryface[key])
        # change indice of faces (+master['faces'].nfaces)
        tmp_startFace = int(new_boundary[key][b'startFace'].decode('UTF-8')) + master['faces'].nfaces
        tmp_startFace = str(tmp_startFace).encode('UTF-8')
        new_boundary[key][b'startFace'] = tmp_startFace
    
    
    # create the merge dcitionnary
    merge = {'points': {'x': new_points_x, 'y': new_points_y, 'z': new_points_z}, 
             'faces': new_faces,
             'owner': new_owner,
             'neighbour': new_neighbour,
             'boundary': new_boundary}
             
    return merge


def write_merge_mesh(path, d, header):
    """ write files 'points', 'faces', 'owner', 'neighbour', 'boundary'
    in 'path' directory
    """
    # write points:
    with open(header, 'r', encoding='utf8') as f:
	    txt = f.read()
    txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       vectorField;
    location    "constant/polyMesh";
    object      points;

}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


{len(d['points']['x'])}
(
'''
    for x, y, z in zip(d['points']['x'], d['points']['y'], d['points']['z']):
        txt += f"({x} {y} {z})\n"

    txt += ''')
    
// ************************************************************************* //
'''
    save_path = os.path.join(path, 'points')
    print('Save points in path ' + save_path)

    with open(save_path, 'w') as f:
        f.write(txt)
    
    # write faces:
    with open(header, 'r', encoding='utf8') as f:
	    txt = f.read()
    txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       faceList;
    location    "constant/polyMesh";
    object      faces;


}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


{len(d['faces'])}
('''
    for i in range(len(d['faces'])):
        nb_pts =  len(d['faces'][i])
        txt += f"\n{nb_pts}({d['faces'][i][0]} "
        for j in range(1, len(d['faces'][i]) - 1):
            txt += f"{d['faces'][i][j]} "
        txt += f"{d['faces'][i][-1]})"

    txt += '''
)
   
// ************************************************************************* //
'''
    save_path = os.path.join(path, 'faces')
    print('Save faces in path ' + save_path)

    with open(save_path, 'w') as f:
        f.write(txt)
    
    
    # write neighbour:
    with open(header, 'r', encoding='utf8') as f:
	    txt = f.read()
    txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    note        "nPoints:{len(d['points']['x'])}  nCells:{int(max([c for c in d['owner']]) + 1)}  nFaces:{len(d['faces'])}  nInternalFaces:{len(d['neighbour'])}";
    class       labelList;
    location    "constant/polyMesh";
    object      neighbour;


}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


{len(d['neighbour'])}
(
'''
    for cell_id in d['neighbour']:
        txt += f"{int(cell_id)}\n"

    txt += '''
)
    
// ************************************************************************* //
'''
    save_path = os.path.join(path, 'neighbour')
    print('Save neighbour in path ' + save_path)

    with open(save_path, 'w') as f:
        f.write(txt)
 
 
    # write owner:
    with open(header, 'r', encoding='utf8') as f:
	    txt = f.read()
    txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    note        "nPoints:{len(d['points']['x'])}  nCells:{int(max([c for c in d['owner']]) + 1)}  nFaces:{len(d['faces'])}  nInternalFaces:{len(d['neighbour'])}";
    class       labelList;
    location    "constant/polyMesh";
    object      owner;


}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


{len(d['owner'])}
(
'''
    for cell_id in d['owner']:
        txt += f"{int(cell_id)}\n"

    txt += '''
)
    
// ************************************************************************* //
'''
    save_path = os.path.join(path, 'owner')
    print('Save owner in path ' + save_path)

    with open(save_path, 'w') as f:
        f.write(txt)       
        
    
    # write boundary:
    with open(header, 'r', encoding='utf8') as f:
	    txt = f.read()
    txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       polyBoundaryMesh;
    location    "constant/polyMesh";
    object      boundary;

}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


{len(d['boundary'])}
(
'''
    for bound in d['boundary']:
        txt += f'''    {bound.decode('UTF-8')}
    {{
'''
        for key in d['boundary'][bound]:
            txt += f"        {key.decode('UTF-8')}        {d['boundary'][bound][key].decode('UTF-8')};\n"    
        txt += '''
    }
    
'''
     
    txt +='''
)
    
// ************************************************************************* //
'''
    save_path = os.path.join(path, 'boundary')
    print('Save boundary in path ' + save_path)

    with open(save_path, 'w') as f:
        f.write(txt)   
        
        
    
if __name__ == "__main__":
    # read argument pass when coalling this file:
    wd_master = sys.argv[1]
    wd_slave = sys.argv[2]

    # create dictionnary with the mesh of master and slave
    master = {"path": os.path.join(wd_master, 'constant/polyMesh')}
    slave = {"path": os.path.join(wd_slave, 'constant/polyMesh')}
    master = readMesh(master['path'], master)
    slave = readMesh(slave['path'], slave)
    
    # create a merge dictionnary: the new mesh after merging
    merge = mergeMesh(master, slave)
    
    verbose = False
    if verbose:
        print('\n nb points:')
        print(f"master: {len(master['points'].values_x)}, slave: {len(slave['points'].values_x)}, merge: {len(merge['points']['x'])}")
        
        print('\n nb faces:')
        print(f"master: {master['faces'].nfaces}, slave: {slave['faces'].nfaces}, merge: {len(merge['faces'])}")

        print('\n nb owner:')
        print(f"master: {master['owner'].nb_faces}, slave: {slave['owner'].nb_faces}, merge: {len(merge['owner'])}")
        
        print('\n nb neighbour:')
        print(f"master: {master['neighbour'].nb_faces}, slave: {slave['neighbour'].nb_faces}, merge: {len(merge['neighbour'])}")
        
        print('\n nb cell:')
        print(f"master: {master['owner'].nb_cell}, slave: {slave['owner'].nb_cell}, merge: {max([c for c in merge['owner']]) + 1}")
        
        print('\n boundary:')
        print(f"master: {[key for key in master['boundary'].boundaryface]}, \
    slave: {[key for key in slave['boundary'].boundaryface]},  merge: {[key for key in merge['boundary']]}")
    

        
        
    # write new mesh files:
    write_merge_mesh(path = master['path'], d = merge, header = path_OF_header)

