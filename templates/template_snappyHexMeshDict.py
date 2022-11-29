#!/usr/bin/env python
# coding: utf-8

'''The purpose of this template is to generate files:
	    - './system/snappyHexMeshDict'
	    - './system/snappyHexMeshDict0'
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
c = read_simulation_properties('c')
h = read_simulation_properties('h')
r = read_simulation_properties('r')

d = 2 * r + h #width of the area covered by the blade

c = Decimal(str(c)) #here we convert float parameters in Decimal type to avoid
h = Decimal(str(h)) #writing of '0.40000000000000000006' and write '0.4' instead 
r = Decimal(str(r))
d = Decimal(str(d))
####################  writing ./system/snappyHexMeshDict  ################
#the first 7 levels of division 

with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Which of the steps to run
castellatedMesh true;
snap            false;
addLayers       false;


// Geometry. Definition of all surfaces. All surfaces are of class
// searchableSurface.
// Surfaces are used
// - to specify refinement for any mesh cell intersecting it
// - to specify refinement for any mesh cell inside/outside/near
// - to 'snap' the mesh boundary to the surface
geometry
{{


    refinementBox4_up
    {{
        type cylinder;
    	point1          (0 0 -1);
   		point2          (0 0 1);
    	radius          {Decimal('2.5')*r};
    }}

    refinementBox4_center
    {{
        type box;
        min  ({Decimal('-2.5')*r} {-h} -1);
        max  ({Decimal('2.5')*r}  {0}  1);
    }}

    refinementBox4_down
    {{
        type cylinder;
    	point1          (0 {-h}  -1);
   		point2          (0 {-h}  1);
    	radius          {Decimal('2.5')*r};
    }}

    refinementBox5
    {{
        type box;
        min  ({Decimal('-3')*r} {Decimal('-3')*r-h} -1);
        max  ({Decimal('31')*r} {Decimal('3')*r}  1);
    }}

    refinementBox6
    {{
        type box;
        min  ({Decimal('-6')*r} {Decimal('-6')*r-h} -1);
        max  ({Decimal('34')*r} {Decimal('6')*r}  1);
    }}

    refinementBox7
    {{
        type box;
        min  ({Decimal('-10')*r} {Decimal('-10')*r-h} -1);
        max  ({Decimal('40')*r}  {Decimal('10')*r}  1);
    }}

    refinementBox8
    {{
        type box;
        min  ({Decimal('-16')*r} {Decimal('-16')*r-h} -1);
        max  ({Decimal('46')*r}  {Decimal('16')*r}  1);
    }}

    refinementBox9
    {{
        type box;
        min  ({Decimal('-26')*r} {Decimal('-26')*r-h} -1);
        max  ({Decimal('56')*r}  {Decimal('26')*r}  1);
    }}

    refinementBox10
    {{
        type box;
        min  ({Decimal('-40')*r} {Decimal('-40')*r-h} -1);
        max  ({Decimal('70')*r}  {Decimal('40')*r}  1);
    }}
    

}}


// Settings for the castellatedMesh generation.
castellatedMeshControls
{{

    // Refinement parameters
    // ~~~~~~~~~~~~~~~~~~~~~

    // If local number of cells is >= maxLocalCells on any processor
    // switches from from refinement followed by balancing
    // (current method) to (weighted) balancing before refinement.
    maxLocalCells 1000000;

    // Overall cell limit (approximately). Refinement will stop immediately
    // upon reaching this number so a refinement level might not complete.
    // Note that this is the number of cells before removing the part which
    // is not 'visible' from the keepPoint. The final number of cells might
    // actually be a lot less.
    maxGlobalCells 30000000;

    // The surface refinement loop might spend lots of iterations refining just a
    // few cells. This setting will cause refinement to stop if <= minimumRefine
    // are selected for refinement. Note: it will at least do one iteration
    // (unless the number of cells to refine is 0)
    minRefinementCells 10;

    // Allow a certain level of imbalance during refining
    // (since balancing is quite expensive)
    // Expressed as fraction of perfect balance (= overall number of cells /
    // nProcs). 0=balance always.
    maxLoadUnbalance 0.10;


    // Number of buffer layers between different levels.
    // 1 means normal 2:1 refinement restriction, larger means slower
    // refinement.
    nCellsBetweenLevels 4;



    // Explicit feature edge refinement
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Specifies a level for any cell intersected by its edges.
    // This is a featureEdgeMesh, read from constant/triSurface for now.
    features
    (
       // {{
          //  file "motorBike.eMesh";
        //    level 6;
      //  }}
    );



    // Surface based refinement
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    // Specifies two levels for every surface. The first is the minimum level,
    // every cell intersecting a surface gets refined up to the minimum level.
    // The second level is the maximum level. Cells that 'see' multiple
    // intersections where the intersections make an
    // angle > resolveFeatureAngle get refined up to the maximum level.

    refinementSurfaces
    {{
     //   motorBike
     //   {{
            // Surface-wise min and max refinement level
           // level (5 6);

            // Optional specification of patch type (default is wall). No
            // constraint types (cyclic, symmetry) etc. are allowed.
          //  patchInfo
            //{{
           //     type wall;
          //      inGroups (motorBikeGroup);
          //  }}
       // }}
    }}

    // Resolve sharp angles
    resolveFeatureAngle 30;


    // Region-wise refinement
    // ~~~~~~~~~~~~~~~~~~~~~~

    // Specifies refinement level for cells in relation to a surface. One of
    // three modes
    // - distance. 'levels' specifies per distance to the surface the
    //   wanted refinement level. The distances need to be specified in
    //   descending order.
    // - inside. 'levels' is only one entry and only the level is used. All
    //   cells inside the surface get refined up to the level. The surface
    //   needs to be closed for this to be possible.
    // - outside. Same but cells outside.

    refinementRegions
    {{

        refinementBox4_up
        {{
            mode inside;
            levels ((1E15 7));
        }}

        refinementBox4_center
        {{
            mode inside;
            levels ((1E15 7));
        }}

        refinementBox4_down
        {{
            mode inside;
            levels ((1E15 7));
        }}

        refinementBox5
        {{
            mode inside;
            levels ((1E15 6));
        }}

        refinementBox6
        {{
            mode inside;
            levels ((1E15 5));
        }}

        refinementBox7
        {{
            mode inside;
            levels ((1E15 4));
        }}

        refinementBox8
        {{
            mode inside;
            levels ((1E15 3));
        }}

        refinementBox9
        {{
            mode inside;
            levels ((1E15 2));
        }}

        refinementBox10
        {{
            mode inside;
            levels ((1E15 1));
        }}
        
    }}


    // Mesh selection
    // ~~~~~~~~~~~~~~

    // After refinement patches get added for all refinementSurfaces and
    // all cells intersecting the surfaces get put into these patches. The
    // section reachable from the locationInMesh is kept.
    // NOTE: This point should never be on a face, always inside a cell, even
    // after refinement.


    locationInMesh ({Decimal('-14.99')*d} 0.0001 0.00001);


    // Whether any faceZones (as specified in the refinementSurfaces)
    // are only on the boundary of corresponding cellZones or also allow
    // free-standing zone faces. Not used if there are no faceZones.
    allowFreeStandingZoneFaces true;
}}



// Settings for the snapping.
snapControls
{{
    //- Number of patch smoothing iterations before finding correspondence
    //  to surface
    nSmoothPatch 3;

    //- Relative distance for points to be attracted by surface feature point
    //  or edge. True distance is this factor times local
    //  maximum edge length.
    tolerance 2.0;

    //- Number of mesh displacement relaxation iterations.
    nSolveIter 30;

    //- Maximum number of snapping relaxation iterations. Should stop
    //  before upon reaching a correct mesh.
    nRelaxIter 5;

    // Feature snapping

        //- Number of feature edge snapping iterations.
        //  Leave out altogether to disable.
        nFeatureSnapIter 10;

        //- Detect (geometric only) features by sampling the surface
        //  (default=false).
        implicitFeatureSnap false;

        //- Use castellatedMeshControls::features (default = true)
        explicitFeatureSnap true;

        //- Detect points on multiple surfaces (only for explicitFeatureSnap)
        multiRegionFeatureSnap false;
}}



// Settings for the layer addition.
//addLayersControls
addLayersControls
{{
    // Are the thickness parameters below relative to the undistorted
    // size of the refined cell outside layer (true) or absolute sizes (false).
    relativeSizes true;

    // Per final patch (so not geometry!) the layer information
    layers
    {{
       // "(lowerWall|motorBike).*"
        //{{
          //  nSurfaceLayers 1;
        //}}
    }}

    // Expansion factor for layer mesh
    expansionRatio 1.0;

    // Wanted thickness of final added cell layer. If multiple layers
    // is the thickness of the layer furthest away from the wall.
    // Relative to undistorted size of cell outside layer.
    // See relativeSizes parameter.
    finalLayerThickness 0.3;

    // Minimum thickness of cell layer. If for any reason layer
    // cannot be above minThickness do not add layer.
    // Relative to undistorted size of cell outside layer.
    minThickness 0.1;

    // If points get not extruded do nGrow layers of connected faces that are
    // also not grown. This helps convergence of the layer addition process
    // close to features.
    // Note: changed(corrected) w.r.t 1.7.x! (didn't do anything in 1.7.x)
    nGrow 0;

    // Advanced settings

    // When not to extrude surface. 0 is flat surface, 90 is when two faces
    // are perpendicular
    featureAngle 60;

    // At non-patched sides allow mesh to slip if extrusion direction makes
    // angle larger than slipFeatureAngle.
    slipFeatureAngle 30;

    // Maximum number of snapping relaxation iterations. Should stop
    // before upon reaching a correct mesh.
    nRelaxIter 3;

    // Number of smoothing iterations of surface normals
    nSmoothSurfaceNormals 1;

    // Number of smoothing iterations of interior mesh movement direction
    nSmoothNormals 3;

    // Smooth layer thickness over surface patches
    nSmoothThickness 10;

    // Stop layer growth on highly warped cells
    maxFaceThicknessRatio 0.5;

    // Reduce layer growth where ratio thickness to medial
    // distance is large
    maxThicknessToMedialRatio 0.3;

    // Angle used to pick up medial axis points
    // Note: changed(corrected) w.r.t 1.7.x! 90 degrees corresponds to 130
    // in 1.7.x.
    minMedialAxisAngle 90;


    // Create buffer region for new layer terminations
    nBufferCellsNoExtrude 0;


    // Overall max number of layer addition iterations. The mesher will exit
    // if it reaches this number of iterations; possibly with an illegal
    // mesh.
    nLayerIter 50;
}}





// Generic mesh quality settings. At any undoable phase these determine
// where to undo.
meshQualityControls
{{
	// Include defaults parameters from master dictionary
	#includeEtc "caseDicts/meshQualityDict"

	//- minFaceWeight (0 -> 0.5)
	minFaceWeight 0.02;


    // Advanced

    //- Number of error distribution iterations
    nSmoothScale 4;
    //- Amount to scale back displacement at error points
    errorReduction 0.75;
}}


// Advanced

// Write flags
writeFlags
(
    scalarLevels
    layerSets
    layerFields     // write volScalarField for layer coverage
);


// Merge tolerance. Is fraction of overall bounding box of initial mesh.
// Note: the write tolerance needs to be higher than this.
mergeTolerance 1e-6;


// ************************************************************************* //'''

path = os.path.join('system', 'snappyHexMeshDict')
print('Save snappyHexMeshDict in path ' + path)

with open(path, 'w') as f:
    f.write(txt)


####################  writing ./system/snappyHexMeshDict1  ###############
#3 additional levels of divisions 


with open(path_OF_header, 'r', encoding='utf8') as f:
	txt = f.read()

txt += f'''
FoamFile
{{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      snappyHexMeshDict;
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Which of the steps to run
castellatedMesh true;
snap            false;
addLayers       false;


// Geometry. Definition of all surfaces. All surfaces are of class
// searchableSurface.
// Surfaces are used
// - to specify refinement for any mesh cell intersecting it
// - to specify refinement for any mesh cell inside/outside/near
// - to 'snap' the mesh boundary to the surface
geometry
{{

   refinementBox2_up
    {{
        type cylinder;
    	point1          (0 0 -1);
   		point2          (0 0 1);
    	radius          {Decimal('1.95')*r};
    }}

    refinementBox2_center
    {{
        type box;
        min  ({Decimal('-1.95')*r} {-h} -1);
        max  ({Decimal('1.95')*r}  {0}  1);
    }}

    refinementBox2_down
    {{
        type cylinder;
    	point1          (0 {-h}  -1);
   		point2          (0 {-h}  1);
    	radius          {Decimal('1.95')*r};
    }}

    refinementBox3_up
    {{
        type cylinder;
    	point1          (0 0 -1);
   		point2          (0 0 1);
    	radius          {Decimal('2.05')*r};
    }}

    refinementBox3_center
    {{
        type box;
        min  ({Decimal('-2.05')*r} {-h} -1);
        max  ({Decimal('2.05')*r}  {0}  1);
    }}

    refinementBox3_down
    {{
        type cylinder;
    	point1          (0 {-h}  -1);
   		point2          (0 {-h}  1);
    	radius          {Decimal('2.05')*r};
    }}


}}


// Settings for the castellatedMesh generation.
castellatedMeshControls
{{

    // Refinement parameters
    // ~~~~~~~~~~~~~~~~~~~~~

    // If local number of cells is >= maxLocalCells on any processor
    // switches from from refinement followed by balancing
    // (current method) to (weighted) balancing before refinement.
    maxLocalCells 1000000;

    // Overall cell limit (approximately). Refinement will stop immediately
    // upon reaching this number so a refinement level might not complete.
    // Note that this is the number of cells before removing the part which
    // is not 'visible' from the keepPoint. The final number of cells might
    // actually be a lot less.
    maxGlobalCells 30000000;

    // The surface refinement loop might spend lots of iterations refining just a
    // few cells. This setting will cause refinement to stop if <= minimumRefine
    // are selected for refinement. Note: it will at least do one iteration
    // (unless the number of cells to refine is 0)
    minRefinementCells 10;

    // Allow a certain level of imbalance during refining
    // (since balancing is quite expensive)
    // Expressed as fraction of perfect balance (= overall number of cells /
    // nProcs). 0=balance always.
    maxLoadUnbalance 0.10;


    // Number of buffer layers between different levels.
    // 1 means normal 2:1 refinement restriction, larger means slower
    // refinement.
    nCellsBetweenLevels 4;



    // Explicit feature edge refinement
    // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

    // Specifies a level for any cell intersected by its edges.
    // This is a featureEdgeMesh, read from constant/triSurface for now.
    features
    (
       // {{
          //  file "motorBike.eMesh";
        //    level 6;
      //  }}
    );



    // Surface based refinement
    // ~~~~~~~~~~~~~~~~~~~~~~~~

    // Specifies two levels for every surface. The first is the minimum level,
    // every cell intersecting a surface gets refined up to the minimum level.
    // The second level is the maximum level. Cells that 'see' multiple
    // intersections where the intersections make an
    // angle > resolveFeatureAngle get refined up to the maximum level.

    refinementSurfaces
    {{
     //   motorBike
     //   {{
            // Surface-wise min and max refinement level
           // level (5 6);

            // Optional specification of patch type (default is wall). No
            // constraint types (cyclic, symmetry) etc. are allowed.
          //  patchInfo
            //{{
           //     type wall;
          //      inGroups (motorBikeGroup);
          //  }}
       // }}
    }}

    // Resolve sharp angles
    resolveFeatureAngle 30;


    // Region-wise refinement
    // ~~~~~~~~~~~~~~~~~~~~~~

    // Specifies refinement level for cells in relation to a surface. One of
    // three modes
    // - distance. 'levels' specifies per distance to the surface the
    //   wanted refinement level. The distances need to be specified in
    //   descending order.
    // - inside. 'levels' is only one entry and only the level is used. All
    //   cells inside the surface get refined up to the level. The surface
    //   needs to be closed for this to be possible.
    // - outside. Same but cells outside.

    refinementRegions
    {{


        
        refinementBox2_up
        {{
            mode inside;
            levels ((1E15 2));
        }}

        refinementBox2_center
        {{
            mode inside;
            levels ((1E15 2));
        }}

        refinementBox2_down
        {{
            mode inside;
            levels ((1E15 2));
        }}
        
        refinementBox3_up
        {{
            mode inside;
            levels ((1E15 1));
        }}

        refinementBox3_center
        {{
            mode inside;
            levels ((1E15 1));
        }}

        refinementBox3_down
        {{
            mode inside;
            levels ((1E15 1));
        }}

        
    }}


    // Mesh selection
    // ~~~~~~~~~~~~~~

    // After refinement patches get added for all refinementSurfaces and
    // all cells intersecting the surfaces get put into these patches. The
    // section reachable from the locationInMesh is kept.
    // NOTE: This point should never be on a face, always inside a cell, even
    // after refinement.

    locationInMesh (0 0.0001 0.00001);



    // Whether any faceZones (as specified in the refinementSurfaces)
    // are only on the boundary of corresponding cellZones or also allow
    // free-standing zone faces. Not used if there are no faceZones.
    allowFreeStandingZoneFaces true;
}}



// Settings for the snapping.
snapControls
{{
    //- Number of patch smoothing iterations before finding correspondence
    //  to surface
    nSmoothPatch 3;

    //- Relative distance for points to be attracted by surface feature point
    //  or edge. True distance is this factor times local
    //  maximum edge length.
    tolerance 2.0;

    //- Number of mesh displacement relaxation iterations.
    nSolveIter 30;

    //- Maximum number of snapping relaxation iterations. Should stop
    //  before upon reaching a correct mesh.
    nRelaxIter 5;

    // Feature snapping

        //- Number of feature edge snapping iterations.
        //  Leave out altogether to disable.
        nFeatureSnapIter 10;

        //- Detect (geometric only) features by sampling the surface
        //  (default=false).
        implicitFeatureSnap false;

        //- Use castellatedMeshControls::features (default = true)
        explicitFeatureSnap true;

        //- Detect points on multiple surfaces (only for explicitFeatureSnap)
        multiRegionFeatureSnap false;
}}



// Settings for the layer addition.
//addLayersControls
addLayersControls
{{
    // Are the thickness parameters below relative to the undistorted
    // size of the refined cell outside layer (true) or absolute sizes (false).
    relativeSizes true;

    // Per final patch (so not geometry!) the layer information
    layers
    {{
       // "(lowerWall|motorBike).*"
        //{{
          //  nSurfaceLayers 1;
        //}}
    }}

    // Expansion factor for layer mesh
    expansionRatio 1.0;

    // Wanted thickness of final added cell layer. If multiple layers
    // is the thickness of the layer furthest away from the wall.
    // Relative to undistorted size of cell outside layer.
    // See relativeSizes parameter.
    finalLayerThickness 0.3;

    // Minimum thickness of cell layer. If for any reason layer
    // cannot be above minThickness do not add layer.
    // Relative to undistorted size of cell outside layer.
    minThickness 0.1;

    // If points get not extruded do nGrow layers of connected faces that are
    // also not grown. This helps convergence of the layer addition process
    // close to features.
    // Note: changed(corrected) w.r.t 1.7.x! (didn't do anything in 1.7.x)
    nGrow 0;

    // Advanced settings

    // When not to extrude surface. 0 is flat surface, 90 is when two faces
    // are perpendicular
    featureAngle 60;

    // At non-patched sides allow mesh to slip if extrusion direction makes
    // angle larger than slipFeatureAngle.
    slipFeatureAngle 30;

    // Maximum number of snapping relaxation iterations. Should stop
    // before upon reaching a correct mesh.
    nRelaxIter 3;

    // Number of smoothing iterations of surface normals
    nSmoothSurfaceNormals 1;

    // Number of smoothing iterations of interior mesh movement direction
    nSmoothNormals 3;

    // Smooth layer thickness over surface patches
    nSmoothThickness 10;

    // Stop layer growth on highly warped cells
    maxFaceThicknessRatio 0.5;

    // Reduce layer growth where ratio thickness to medial
    // distance is large
    maxThicknessToMedialRatio 0.3;

    // Angle used to pick up medial axis points
    // Note: changed(corrected) w.r.t 1.7.x! 90 degrees corresponds to 130
    // in 1.7.x.
    minMedialAxisAngle 90;


    // Create buffer region for new layer terminations
    nBufferCellsNoExtrude 0;


    // Overall max number of layer addition iterations. The mesher will exit
    // if it reaches this number of iterations; possibly with an illegal
    // mesh.
    nLayerIter 50;
}}





// Generic mesh quality settings. At any undoable phase these determine
// where to undo.
meshQualityControls
{{
	// Include defaults parameters from master dictionary
	#includeEtc "caseDicts/meshQualityDict"

	//- minFaceWeight (0 -> 0.5)
	minFaceWeight 0.02;


    // Advanced

    //- Number of error distribution iterations
    nSmoothScale 4;
    //- Amount to scale back displacement at error points
    errorReduction 0.75;
}}


// Advanced

// Write flags
writeFlags
(
    scalarLevels
    layerSets
    layerFields     // write volScalarField for layer coverage
);


// Merge tolerance. Is fraction of overall bounding box of initial mesh.
// Note: the write tolerance needs to be higher than this.
mergeTolerance 1e-6;


// ************************************************************************* //'''



path = os.path.join('system', 'snappyHexMeshDict1')
print('Save snappyHexMeshDict1 in path ' + path)

with open(path, 'w') as f:
    f.write(txt)
