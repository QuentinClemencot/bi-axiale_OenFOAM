#!/usr/bin/env python
# coding: utf-8


from scipy import integrate
from scipy import integrate
'''
The purpose of this script is :
- The calculation of section surface **section**
- The calculation of the position of G (center of gravity for a 2D NACA airfoil)
- The calculation of the moment of inertia of the blade at G points. 

It plot result

It python version is 3.8.2

'''

###################### Define parameters #######################################
verification = False #if ==True, verification step will be launch 
if verification: c = 1 #[m] chord lenght

 
####################### Define functions #######################################      
def y_OO12(x):
    """ return thickness of the blade at position x (x=0 corresponds to the 
        leading edge, x=1 corresponds to the trailing edge).
        The top part of the airfoil is given by yt(x) = y(x).
        The bottom part of the airfoil is given by yt(x) = -y(x).
    """
    A =  0.2969 #NACA 4 digits constants
    B = -0.1260
    C = -0.3516
    D =  0.2843
    E = -0.1015
    
    t = 0.12
    c = 1
    
    return 5*t*(A*(x/c)**0.5 
               + B*(x/c)
               + C*(x/c)**2
               + D*(x/c)**3
               + E*(x/c)**4)

def section_OO12(c):
    """ return profil section [m²] for NACA0012
        airfoil.
    """
    return integrate.dblquad(lambda y, x: 1, 0, c, lambda x: -c*y_OO12(x/c), \
                                                   lambda x:  c*y_OO12(x/c))[0]


def xG_position_OO12(c, verbose=False):
    """ return position of center of gravity for NACA0012
        airfoil.
    """
    section = integrate.dblquad(lambda y, x: 1, 0, c, lambda x: -c*y_OO12(x/c), \
                                                      lambda x:  c*y_OO12(x/c))
                                                      
    int_xdS = integrate.dblquad(lambda y, x: x, 0, c, lambda x: -c*y_OO12(x/c), \
                                                      lambda x:  c*y_OO12(x/c))

    xG = int_xdS[0]/section[0]

    if verbose:
        print(f"Section of the profil is {section[0]:6f}m^2 +- {section[1]}.\n")
        print(f"The distance between the leading edge and the center of gravity is"+
              f"{xG:6f}m, or {xG/c:6f}c. \n")
    
    return xG
    
    

def area_moment_inertia_OO12(xG, c):
    """ return area moment of inertia around axe (G, ey) with G: center of gravity and
        ey: vector orhtogonal to the study plane 
    Args:
        xG: float (x-coordonate of center of gravity, x=0 being the leading edge
                   x=c being the trailing edge [m])\n
        c: float  (chord lenght [m])\n
    Returns:
        J: float (inertia momentum [m⁴])\n
    """
    if xG > c:
        raise ValueError('Warning: xG > c, center of gravity outside blade !')
    else:
        return integrate.dblquad(lambda y, x: (x-xG)**2 + y**2,\
                                     0, c,\
                                     lambda x: -c*y_OO12(x/c),\
                                     lambda x:  c*y_OO12(x/c))


####################### Main ###################################################  
if verification:
	xG = xG_position_OO12(c) #x-coordonate of center of gravity, x=0 being the leading edge
		                #x=c being the trailing edge [m]

	J = area_moment_inertia_OO12(xG, c) #tuple with moment of inertia and 
		                                   #precision on integration

	print(f"Area moment of inertia around (G,ey) axe is J={J[0]:6f} m⁴ +- {J[1]}.\n")

	#verification : for a rectangle of lenght b (x-direction) and height h (y-direction)
	#with centroid at the origin, J = (bh/12) * (b**2 + h**2)
	#let test it with our method for b=1 and h=0.12, exact answer is J=0.010144m⁴:
	L, l = 1, 0.12
	J_rec = integrate.dblquad(lambda y, x: (x-0.5)**2 + (y-0.06)**2,\
		                                 0, L,\
		                                 lambda x: 0,\
		                                 lambda x: l)
	print(f"Area moment of inertia of rectangle with L=1m and l=0.12m  is J={J_rec[0]:6f} m⁴ +- {J_rec[1]}.\n \
	Exact solution is J={l*L*(l**2 + L**2)/12}m⁴")
       
       
       
       
                                          
                                     
