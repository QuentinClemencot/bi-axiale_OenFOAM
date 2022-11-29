""" This file contain the main parameter to generate the fields views"""

#### list of parameters:

# scene view settings:
Ambient = 0.8
Diffuse = 0.5

# scale fo field:
Scale_vorticity = [-80, 80]
Scale_p = [-20, 20]
Scale_Ux = [-1, 3]

# field to read:
isVectorField_vorticity = True
Field_vorticity = ('vorticity', 'Z') 

isVectorField_p = False
Field_p = ('p')

isVectorField_Ux = True
Field_Ux = ('U', 'X') 

# image definition
pixel_width = 800

# image demention
rotor_height = 1.2 # times d = h + 2*r
rotor_width = 1.2 # times d = h + 2*r
zoom_height = 3 # times c
zoom_width = 3 # times c
