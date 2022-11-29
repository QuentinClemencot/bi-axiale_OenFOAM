""" This file contain the main carateristics of the simulationInput
 It is use to fill template files with correct values"""

from math import pi




'''Blade caracteristics'''
Nb=2 #number of blade
c=0.04  #chord lenght [m]
b=0.01 #wingspan of blade (for 2D case=cell thickness in z-direction) [m]
Aref=c*b #reference surface for the force coefficient calculation Aref=chord*width [m²]
xa=1/8  #position of the leading axis of rotation (point A) on the blade  (xa=0: leading edges, xa=1: trailing edge) []

'''Flow caracteristics'''
rho=1000 #density of fluid [kg/m³] 
Pref=0 # Kinematic pressure of reference [Pa/(kg/m³)] 
Uref=1    # Velocity impose in inlet [m/s]

#Here we either use the Re or nu parameters. The parameter not use is set =0 and
#will not be consider during pre processing.
Re=1 * 10**6 #Reynolds number  []
nu=0 #cinematic viscosity [m²/s] 


'''Turbulence caracteristics'''
I=0.01  #Turbulence intensity  []
Lturb=0.01 * c  #Turbulence reference scale [m]
Cmu=0.09 # modele constant https://www.openfoam.com/documentation/guides/latest/doc/guide-turbulence-ras-k-omega-sst.html

'''CofR trajectory caracteristics'''
r=0.1 #radius of rotation [m]
h=0.4 #vertical translation distance of the foil [m]
Lambda=5 #dimensionless foil velocity (Lambda=V/Uref with V=velocity of foil) []


'''blade orientation caracteristics'''
lawType=3 #=0 if the blade has a translation mouvement (fixed theta=-beta0-pi/2)
		  #=1 if the blade orientation is dict by setting angle : fixed beta=beta0
		  #=2 if the blade orientation is dict by setting angle : beta=-beta0 descending translation
				#beta=beta0  ascending translation, a third degree polynomial law dict the
				# transition between both part
		  #=3 if the blade orientation is dict by setting angle : beta=0 in translation part
				#2 points follow the trajectory  : the CofR (A) and
				# a point close to the trailing edge (B). 
					
beta0=0 #setting angle depend on lawType. 
		 #Requisite for lawType 0,1,2

xb=7/8 #Near trailing egde point position (B point).xb>xa (xb=0 : leading edges, xb=1 : trailing edge) 
 	   #Requisite for lawType 3

'''Mesh size caracteristics'''
cellmin=0.0016 #[m] dimension of the smallest cells of the background (approximatively the overset biggest cells)

'''Time caracteristics'''
nStepStart=50# number of time step of the "starting state" []

nRevRegular=18 # number of revolutions made (to reach permanent state) []
# between 2 time step, the overset must move of ~cellmin and of ~1° in the turn
# so we compute the 2 times step to respect eatch conditions and we take the smaller
l = 2 * (pi * r + h) # lenght of blade trajectory [m]
dtheta = 1 # target azimuthal variation in turns between 2 times step [°]
dl = min(dtheta * 2 * pi * r / 360, cellmin) # distance travel by blade between 2 time step [m]
nStep_per_revRegular=int(l/dl) # integer, number of time step per period []

nRevProbes=2 # number of extra revolutions made to collect probes []
nStep_per_revProbes= 2 * nStep_per_revRegular # number of time step per period []

T=round(2*(h+pi*r)/(Uref*Lambda),10) #effective period (round to avoid numerical precision issues)[s]
deltaT_regular=round(2*(h+pi*r)/(Uref*Lambda*nStep_per_revRegular),10) #effective time step (round to avoid numerical precision issues)[s]
deltaT_probes=round(2*(h+pi*r)/(Uref*Lambda*nStep_per_revProbes),10)
deltaT_start=deltaT_regular/2

'''Sampling caracteristics'''
writeInterval_surface=10 #write p and wss along blades surface every 'writeInterval_surface' time step.

writeInterval_timeDir_reg=600 #write time directory every "writeInterval_timeDir" time step
writeInterval_timeDir_probes=120 #write time directory every "writeInterval_timeDir" time step





'''Decomposition of the domain:
 ___________________________________
 |                                 |
 |               zone1             |
 |                                 |
 |-------A-----B-------------------|
 |       |     |                   |
 | zone2 |  3  |     zone4         |
 |       |     |                   |
 |-------C-----D-------------------|
 |                                 |
 |          zoone5                 |
 |_________________________________|

 oversets are in zone 3.'''

pointA = (-r - 2 * c, r + 2 * c) # coordonates of point A [m]
pointB = ( r + 2 * c, r + 2 * c) #x coordonates of point B [m]
pointC = (-r - 2 * c, -r -h -2 * c) # coordonates of point A [m]
pointD = ( r + 2 * c, -r -h -2 * c) #x coordonates of point B [m]

meshPlane = 'xy' # plane containing the 2D mesh: 'xy', 'xz' or 'yz'

nzone1 = (1,  1) #proc decomposition of zone 1 (method simple) (2D)
nzone2 = (1,  1) #proc decomposition of zone 2 (method simple)
nzone3 = (2,  3) #proc decomposition of zone 3 (method simple)
nzone4 = (1,  1) #proc decomposition of zone 4 (method simple)
nzone5 = (1,  1) #proc decomposition of zone 5 (method simple)

tolerance = 2 #[%] tolerance of non-homogeneity of cell distribution for proc of zone 3
nproc = 10 	# total number of proc, it will be recompute and compared with
			# this user estimation for verification.

