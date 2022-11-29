#!/usr/bin/env python

from general_function_def import *
import os
import numpy as np
from math import *
import copy




def add_CofR(data,h,r):
    """ read data containing 's' key and geometric parameters :
        h : length of vertical translation\n
        r : radius of rotation \n
    return a dictionnary with same keys than data plus the key 'CofR'
        which is a sub-dictionnary containing 2 keys 'x' and 'y'.
        x' and 'y' are 1D array which give the postion of the center of
        rotation (CofR) of the foil.
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation.

    Args:
        data: {}\n
        h: float\n
        r: float\n

    Returns:
        dataUp : {}\n

    A way you might use me is:\n
        data = add_CofR(data,h,r)

    """
    if 's' not in data:
            raise ValueError(' data does not have "s" key')

    dataUp = copy.deepcopy(data)
    size=len(dataUp['s'])

    dataUp['CofR']={'x':np.zeros(size),'y':np.zeros(size)}
    dataUp['dCofR_dt']={'x':np.zeros(size),'y':np.zeros(size)}
    dataUp['d2CofR_dt2']={'x':np.zeros(size),'y':np.zeros(size)}

    for i in range(size):
        x=x_coord(dataUp['s'][i],h,r)
        

        dataUp['CofR']['x'][i]=x_coord(dataUp['s'][i],h,r)
        dataUp['CofR']['y'][i]=y_coord(dataUp['s'][i],h,r)
        

    return (dataUp)


def x_coord(s,h,r):
    """ read several parameters :
        s : position parameter 0<=s<=1\n
        h : length of vertical translation\n
        r : radius of rotation \n
    y-axe is the vertical ascending direction
    x-axe is  pointing in downstream direction 
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation. it coordonate are (x,y)=(-r,0).
    Returns the x-coordonate of the point at position s.

    Args:
        s : float\n
        h : float\n
        r : float\n

    Returns:
        x : float\n

    A way you might use me is:\n
        x=x_coord(s,h,r)

    """
    if s<0 or s>1 :
        raise ValueError(' s='+str(s)+' parameter is not 0<=s<=1')


    #distance traveled by the center of rotation in one period :
    l=2*(h+pi*r)

    if s<h/l :
        #downward vertical translation part :
        x=-r

    elif s>=h/l and s<0.5 :
        #first rotation of pi :
        alpha=(l*s-h)/r  #alpha between 0 and pi
        x=-r*cos(alpha)


    elif s>=0.5 and s<0.5+h/l :
        #vertical upward translation :
        x=r

    elif s>=0.5+h/l :
        #second rotation of pi  :
        alpha=(l*(s-0.5)-h)/r  #alpha between 0 and pi/2
        x=r*cos(alpha)


    else : 
        raise ValueError('s='+str(s)+' does not belong to any defined trajectory sections')

    return (x)


def y_coord(s,h,r):
    """ read several parameters :
        s : position parameter 0<=s<=1\n
        h : length of vertical translation\n
        r : radius of rotation \n
    y-axe is the vertical ascending direction
    x-axe is  pointing in downstream direction 
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation. it coordonate are (x,y)=(-r,0).
    Returns the y-coordonate of the point at position s.

    Args:
        s : float\n
        h : float\n
        r : float\n

    Returns:
        y : float\n

    A way you might use me is:\n
        y=y_coord(s,h,r)

    """
    if s<0 or s>1 :
        raise ValueError(' s parameter is not 0<=s<=1')


    #distance traveled by the center of rotation in one period :
    l=2*(h+pi*r)

    if s<h/l :
        #downward vertical translation part :
        y=-s*l

    elif s>=h/l and s<0.5 :
        #first rotation of pi :
        alpha=(l*s-h)/r  #alpha between 0 and pi
        y=-h-r*sin(alpha)


    elif s>=0.5 and s<0.5+h/l :
        #vertical upward translation :
        y=-(0.5+h/l-s)*l

    elif s>=0.5+h/l :
        #second rotation of pi  :
        alpha=(l*(s-0.5)-h)/r  #alpha between 0 and pi/2
        y=r*sin(alpha)

    else : 
        raise ValueError('s='+str(s)+' does not belong to any defined trajectory sections')

    return (y)



def add_gamma(data,h,r):
    """ read data containing 's' key and geometric parameters :
        h : length of vertical translation\n
        r : radius of rotation \n
    return a dictionnary with same keys than data plus the key 'gamma'
        which is a 1D array which gives the angle between the tangent to the 
        trajectory and the x-axe direction. gamma is in [-pi,+pi]
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation.

    Args:
        data: {}\n
        h: float\n
        r: float\n

    Returns:
        dataUp : {}\n

    A way you might use me is:\n
        data = add_gamma(data,h,r)

    """
    if 's' not in data:
            raise ValueError(' data does not have "s" key')

    dataUp = copy.deepcopy(data)
    size=len(dataUp['s'])

    dataUp['gamma']=np.zeros(size)

    for i in range(size):
        dataUp['gamma'][i]=gammaAngle(dataUp['s'][i],h,r)
        

    return (dataUp)

def add_gammaB(data,h,r):
    """ read data containing 'sB' key and geometric parameters :
        h : length of vertical translation\n
        r : radius of rotation \n
    return a dictionnary with same keys than data plus the key 'Bgamma'
        which is a 1D array which gives the angle between the tangent to the 
        trajectory and the x-axe direction for B point. gamma is in [-pi,+pi]
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation.

    Args:
        data: {}\n
        h: float\n
        r: float\n

    Returns:
        dataUp : {}\n

    A way you might use me is:\n
        data = add_gammaB(data,h,r)

    """
    if 'sB' not in data:
            raise ValueError(' data does not have "s" key')

    dataUp = copy.deepcopy(data)
    size=len(dataUp['s'])

    dataUp['gammaB']=np.zeros(size)

    for i in range(size):
        dataUp['gammaB'][i]=gammaAngle(dataUp['sB'][i],h,r)
        

    return (dataUp)
    
'''
def add_gammaB(data,h,r):
    """ read data containing 'sB' key and geometric parameters :
        h : length of vertical translation\n
        r : radius of rotation \n
    return a dictionnary with same keys than data plus the key 'gammaB'
        which is a 1D array which gives the angle between the tangent to the 
        trajectory and the x-axe direction of B point. gamma is in [-pi,+pi]
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation.

    Args:
        data: {}\n
        h: float\n
        r: float\n

    Returns:
        dataUp : {}\n

    A way you might use me is:\n
        data = add_gammaB(data,h,r)

    """
    if 's' not in data:
            raise ValueError(' data does not have "sB" key')

    dataUp = copy.deepcopy(data)
    size=len(dataUp['sB'])

    dataUp['gammaB']=np.zeros(size)

    for i in range(size):
        dataUp['gammaB'][i]=gammaAngle(dataUp['sB'][i],h,r)
        
    return (dataUp)
'''

def gammaAngle(s,h,r):
    """ read several parameters :
        s : position parameter 0<=s<=1\n
        h : length of vertical translation\n
        r : radius of rotation \n
    y-axe is the vertical ascending direction
    x-axe is  pointing in downstream direction 
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation. it coordonate are (x,y)=(-r,0).
    Returns the gamma angle of the point at position s. gamma is in [-pi,+pi].

    Args:
        s : float\n
        h : float\n
        r : float\n

    Returns:
        gamma : float\n

    A way you might use me is:\n
        gamma=gammaAngle(s,h,r)

    """
    if s<0 or s>1 :
        raise ValueError(' s parameter is not 0<=s<=1')


    #distance traveled by the center of rotation in one period :
    l=2*(h+pi*r)

    if s<h/l :
        #downward vertical translation part :
        gamma=-pi/2

    elif s>=h/l and s<0.5 :
        #first rotation of pi :
        alpha=(l*s-h)/r  #alpha between 0 and pi
        gamma=-pi/2+alpha


    elif s>=0.5 and s<0.5+h/l :
        #vertical upward translation :
        gamma=pi/2

    elif s>=0.5+h/l :
        #second rotation of pi  :
        alpha=(l*(s-0.5)-h)/r  #alpha between 0 and pi/2
        gamma=pi/2+alpha
        #-pi<gamma<pi
        if gamma>pi : gamma=gamma-2*pi

    else : 
        raise ValueError('s='+str(s)+' does not belong to any defined trajectory sections')
    
    return (gamma)



def add_theta(beta0, c, xl, xc, data, h, lawType, r, verbose=False):
    """ read data containing 's', 'gamma' keys and geometric parameters :
        h : length of vertical translation\n
        r : radius of rotation \n
        beta0 : basic setting angle \n
        lawType : integer give the setting law \n
        c : chord lenght \n
        xl : position of B point (0:leading point, 1:trailing point) \n
        xc : position of CofR (point A) (0:leading point, 1:trailing point) \n
        verbose: Booleen to active prints
    return a dictionnary with same keys than data plus the key 'theta'
        which is a 1D array which gives the angle between x direction and the chord.
        theta is express in radian and -pi/2<theta<pi/2
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation.

    Args:
        data: {}\n
        h: float\n
        r: float\n
        beta0: float\n
        lawType: int\n
        c: float\n
        xl: float\n
        xc: float\n
        verbose: booleen\n

    Returns:
        dataUp : {}\n


    """
    if 's' not in data:
            raise ValueError(' data does not have "s" key')
    max_s = max(data['s'])
    min_s = min(data['s'])
    if verbose: print(f"max(data[s])={max_s}\nmin(data[s])={min_s}")
    if max_s > 1:
            raise ValueError('max( data[s])>1')
    if min_s > 1:
            raise ValueError('min( data[s])<1')
    if 'gamma' not in data:
            raise ValueError(' data does not have "gamma" key')

    beta0 = DegToRad(beta0) #beta0 is in degre, we turn it in radian

    dataUp = copy.deepcopy(data)
    max_s = max(dataUp['s'])
    min_s = min(dataUp['s'])
    if verbose: print(f"max(data[s])={max_s}\nmin(data[s])={min_s}")
    size = len(dataUp['s'])

    dataUp['theta'] = np.zeros(size)


    if lawType == 0: 
        dataUp['theta'] = thetaLawType0(beta0, size)

    if lawType == 1: 
        dataUp['theta'] = thetaLawType1(dataUp['gamma'], beta0, size)

    if lawType == 2: 
        dataUp['theta'] = thetaLawType2(h, r, dataUp['gamma'], dataUp['s'], beta0, size)

    if lawType == 3:
        if verbose: print(f"lawType = {lawType}") 
        dataUp['B'] = {}
        (dataUp['theta'],dataUp['B']['x'],dataUp['B']['y'])=thetaLawType3(c,dataUp['CofR']['x'],dataUp['CofR']['y'],dataUp['gamma'],h,r,dataUp['s'],size,xl,xc)

    if verbose: print(f"'B' in dataUp: {'B' in dataUp}")
    return (dataUp)




def thetaLawType0(beta0,size):
    """ the blade has a translation mouvement ie  theta is constant
        theta(s)=-beta0-pi/2.

    Returns 1D array of size "size" containing values of theta angle.

    Args:
        beta0: float\n
        size: int\n

    Returns:
        beta : float\n

    A way you might use me is:\n
        theta=thetaLawType0(beta0,size)

    """
    return (np.ones(size)*(-beta0-pi/2))




def thetaLawType1(gamma,beta0,size):
    """ the blade orientation is dict by setting angle : fixed beta=beta0

    Returns 1D array containing values of theta angle.

    Args:
        gamma: 1D array of float\n
        beta0: float\n
        size: int\n

    Returns:
        theta: 1D array of float\n

    A way you might use me is:\n
        theta=thetaLawType1(gamma,size)

    """
    theta=np.zeros(size)
    for i in range(size):
        theta[i]=gamma[i]+beta0

    return (theta)




def thetaLawType2(h,r,gamma,s,beta0,size):
    """the blade orientation is dict by setting angle : 
            beta=-beta0 descending translation
            beta=beta0  ascending translation     
            a smooth transition is necessary to avoid discontinuity of moment and forces  
            a third degree polynomial law dict the transition between both part

    Returns 1D array containing values of theta angle.

    Args:
        h: float\n
        r: float\n
        gamma: 1D array of float\n
        s: 1D array of float\n
        beta0: float\n
        size: int\n

    Returns:
        theta: 1D array of float\n

    A way you might use me is:\n
        theta=thetaLawType2(h,r,gamma,beta0,size)

    """

    #distance traveled by the center of rotation in one period :
    l=2*(h+pi*r)

    #a smooth transition is necessary in turn to avoid discontinuity of moment and forces
    #theta is defined in turn by third order polynome : theta(s)=as³+bs²+cs+d
    #(a,b,c,d) is calculated in order to respect boundary condition at the begining and
    #the end of turn on theta and it derivative.

    #First we resolve linear system to obtain (a,b,c,d) coeff for bottom turn :  
    s0=h/l               #value of s at begining of bottom turn
    theta_0=-pi/2-beta0  #value of theta at begining of bottom turn
    dtheta_0=0           #value of d theta/ds at begining of bottom turn
    s1=0.5               #value of s at end of bottom turn
    theta_1=pi/2+beta0   #value of theta at end of bottom turn
    dtheta_1=0           #value of dtheta/dt at end of bottom turn

    (a_dwn,b_dwn,c_dwn,d_dwn)=find_3orderPolCoeff(s0,theta_0,dtheta_0,s1,theta_1,dtheta_1)
    #print('a='+str(a_dwn)+'b='+str(b_dwn)+'c='+str(c_dwn)+'d='+str(d_dwn))

    #Then we resolve linear system to obtain (a,b,c,d) coeff for top turn :  
    s0=h/l+0.5               #value of s at begining of top turn
    theta_0=pi/2+beta0  #value of theta at begining of top turn
    dtheta_0=0           #value of d theta/ds at begining of top turn
    s1=1               #value of s at end of top turn
    theta_1=1.5*pi-beta0   #value of theta at end of top turn
    dtheta_1=0           #value of dtheta/dt at end of top turn

    (a_top,b_top,c_top,d_top)=find_3orderPolCoeff(s0,theta_0,dtheta_0,s1,theta_1,dtheta_1)

    theta=np.zeros(size)

    for i in range(size):
        if s[i]<h/l :
            #downward vertical translation part :
            theta[i]=gamma[i]-beta0

        elif s[i]>=h/l and s[i]<0.5 :
            #first rotation of pi.
            theta[i]=a_dwn*s[i]**3+b_dwn*s[i]**2+c_dwn*s[i]+d_dwn

        elif s[i]>=0.5 and s[i]<0.5+h/l :
            #vertical upward translation :
            theta[i]=gamma[i]+beta0

        elif s[i]>=0.5+h/l :
            #second rotation of pi  :
            theta[i]=a_top*s[i]**3+b_top*s[i]**2+c_top*s[i]+d_top

    return (theta)


def thetaLawType3(c,x,y,gamma,h,r,s,size,xl,xc):
    """the blade orientation is dict by setting angle : beta=0 in translation part.
            a smooth transition between translation and rotation part is necessary to avoid
            discontinuity of moment and forces.
            2 points follow the trajectory of the CofR  : the CofR it self (A) and
            a point close to the trailing edge (B).

    Returns 1D array containing values of theta angle.

    Args:
        c: float\n
        r: float\n
        gamma: 1D array of float\n
        s: 1D array of float\n
        beta0: float\n
        size: int\n

    Returns:
        theta: 1D array of float\n

    A way you might use me is:\n
        theta=thetaLawType2(h,r,gamma,beta0,size)

    """

    #distance traveled by the center of rotation in one period :
    l=2*(h+pi*r)

    cbis=c*(xl-xc) #distance AB
    arc=2*asin(0.5*cbis/r)*r #arc lenght between A and B when both point are 
                                #on turn part

    theta=np.zeros(size)
    yBtab=np.zeros(size)
    xBtab=np.zeros(size)
    #print('h/l='+str(h/l))
    #print('(h+arc)/l='+str((h+arc)/l))

    for i in range(size):
        if s[i]>1:
            raise ValueError('s[i]='+str(s[i])+'>1...!')


        if h == 0:
            #there is no translation part (Darrieus type)
            sB = s[i] - arc/l if s[i] - arc/l >= 0 else 1 + s[i] - arc/l
            #print(f"arc/l = {arc/l}; sA = {s[i]}; s[i] - arc/l == {s[i] - arc/l}; sB = {sB}")
            xB = x_coord(sB, h, r)
            yB = y_coord(sB, h, r)
            if x[i] > xB:
                theta[i] = atan((y[i] - yB) / (x[i] - xB))
            elif x[i] < xB:
                theta[i] = pi + atan((y[i] - yB) / (x[i] - xB))
            elif x[i] == xB and y[i] > yB:
                theta[i] = pi/2
            elif x[i] == xB and y[i] < yB:
                theta[i] = 3*pi/2                
            else :
                raise ValueError('no solution with h=0 for B point')
            yBtab[i] = yB
            xBtab[i] = xB
        elif  s[i]<cbis/l :
            #A on vertical downward translation and B still on turn :
            xA=x[i]
            yA=y[i]
            xO=0
            yO=0
            
            (x3,y3,x4,y4)= coord_intersec_2circles(xO,yO,r,xA,yA,cbis)
            
            #only the solution with y>0 is the correct B point
            #print('y3='+str(y3))
            #print('y4='+str(y4))
            if y3>=0:
                yB=y3
                xB=x3
            elif y4>=0:
                yB=y4
                xB=x4
            elif y3>=0 and y4>=0 and y3!=y4:
                raise ValueError('2 differents intersection points in the half plane y>0 :'
                        +' can not define which is the good one...!')
            else :
                raise ValueError('no solution with y>0 for B point')
    

            theta[i]=-pi+atan(abs(yA-yB)/abs(xA-xB))
            yBtab[i]=yB
            xBtab[i]=xB

        elif s[i]>=cbis/l and s[i]<h/l :
            #downward vertical translation part :
            theta[i]=gamma[i]
            yBtab[i]=y[i]+cbis
            xBtab[i]=x[i]

        elif s[i]>=h/l and s[i]<(h+arc)/l :
            #CofR (A) is on turn part, B is stil on downward vertical translation part.
            #print('A in turn, B is downward vertical')

            xA=x[i]
            yA=y[i]
            (x3,y3,x4,y4)= coord_intersec_circle_line(xA,yA,cbis,-r,0,-r,-h)
            #the solution with y>yO is the correct B point:
            if y3>=-h:
                yB=y3
                xB=x3
            elif y4>=-h:
                yB=y4
                xB=x4
            elif y3>=-h and y4>=-h and y3!=y4:
                raise ValueError('2 differents intersection points in the half plane y>-h :'
                        +' can not define which is the good one...!')
            else :
                raise ValueError('no solution with y>-h for B point')
            #print('x3='+str(x3))
            #print('y3='+str(y3))
            #print('x4='+str(x4))
            #print('y4='+str(y4))
            theta[i]=-pi/2+atan(abs(xA-xB)/abs(yA-yB))
            yBtab[i]=yB
            xBtab[i]=xB
            
            


        elif s[i]>=(h+arc)/l and s[i]<0.5 :
            #both A and B point are on turn part
            sA=s[i]
            sB=sA-arc/l
            xB=x_coord(sB,h,r)
            yB=y_coord(sB,h,r)
            xA=x[i]
            yA=y[i]
            if yB>yA:
                theta[i]=-pi/2+atan(abs(xA-xB)/abs(yA-yB))
            else :
                theta[i]=atan(abs(yA-yB)/abs(xA-xB))
            yBtab[i]=yB
            xBtab[i]=xB

        elif s[i]>=0.5 and s[i]<0.5+cbis/l :
            #A on vertical upward translation an B still on turn :
            xA=x[i]
            yA=y[i]
            xO=0
            yO=-h
            
            (x3,y3,x4,y4)= coord_intersec_2circles(xO,yO,r,xA,yA,cbis)
            
            #only the solution with y<yO is the correct B point
            #print('y3='+str(y3))
            #print('y4='+str(y4))
            if y3<=yO:
                yB=y3
                xB=x3
            elif y4<=yO:
                yB=y4
                xB=x4
            elif y3<=yO and y4<=yO and y3!=y4:
                raise ValueError('2 differents intersection points in the half plane y>0 :'
                        +' can not define which is the good one...!')
            else :
                raise ValueError('no solution with y>0 for B point')
    

            theta[i]=atan(abs(yA-yB)/abs(xA-xB))
            yBtab[i]=yB
            xBtab[i]=xB

        elif s[i]>=0.5+cbis/l and s[i]<0.5+h/l :
            #upward vertical translation part :
            theta[i]=gamma[i]
            yBtab[i]=y[i]-cbis
            xBtab[i]=x[i]

        elif s[i]>=0.5+h/l and s[i]<0.5+(h+arc)/l :
            #CofR (A) is on turn part, B is stil on upward vertical translation part :

            xA=x[i]
            yA=y[i]
            (x3,y3,x4,y4)= coord_intersec_circle_line(xA,yA,cbis,r,0,r,-h)
            #the solution with y<0 is the correct B point:
            if y3<=0:
                yB=y3
                xB=x3
            elif y4<=0:
                yB=y4
                xB=x4
            elif y3<=0 and y4<=0 and y3!=y4:
                raise ValueError('2 differents intersection points in the half plane y>-h :'
                        +' can not define which is the good one...!')
            else :
                raise ValueError('no solution with y>-h for B point')
            #print('x3='+str(x3))
            #print('y3='+str(y3))
            #print('x4='+str(x4))
            #print('y4='+str(y4))
            theta[i]=pi-atan(abs(yA-yB)/abs(xA-xB))
            yBtab[i]=yB
            xBtab[i]=xB
            

        elif s[i]>=0.5+(h+arc)/l :
            #both A and B point are on turn part  :
            sA=s[i]
            sB=sA-arc/l
            #print('sA='+str(sA))
            #print('sB='+str(sB))
            #print('arc='+str(arc))
            #print('l='+str(l))
            xB=x_coord(sB,h,r)
            yB=y_coord(sB,h,r)
            xA=x[i]
            yA=y[i]
            if yA>yB:
                theta[i]=pi-atan(abs(yA-yB)/abs(xA-xB))
            else :
                theta[i]=-pi+atan(abs(yA-yB)/abs(xA-xB))
            yBtab[i]=yB
            xBtab[i]=xB
    #print(f"len(xBtab) = {len(xBtab)}; len(yBtab) = {len(yBtab)}")
    return (theta,xBtab,yBtab)




def add_settingAngle(data):
    """ read data containing 's', 'gamma' and 'theta'keys.
    return a dictionnary with same keys than data plus the key 'beta'
        which is a 1D array which gives the angle between the chord and the 
        tangent to the trajectory (called setting angle).
        beta=theta-gamma

    Args:
        data: {}\n

    Returns:
        dataUp : {}\n

    A way you might use me is:\n
        data = add_settingAngle(data)

    """
    if 's' not in data:
            raise ValueError(' data does not have "s" key')
    if 'gamma' not in data:
            raise ValueError(' data does not have "gamma" key')
    if 'theta' not in data:
            raise ValueError(' data does not have "theta" key')

    dataUp = copy.deepcopy(data)
    size=len(dataUp['s'])

    dataUp['beta']=np.zeros(size)

    for i in range(size):
        dataUp['beta'][i]=dataUp['thetaTot'][i]-dataUp['gammaTot'][i]    

    return (dataUp)



    

def add_Urel(data,h,r,w,Uref,Lambda):
    """ read data containing 's' key and geometric parameters :
        h : length of vertical translation\n
        r : radius of rotation \n
        w : lenght of horizontal translation\n

        Uref is the velocity impose far upstream of the foil. 
        Uref is the component along x direction (the y-component is nul)

        Lambda is the dimensionless velocity magnitude of the center of rotation 
            Lambda=V/Uref (V being the velocity magnitude of the center of rotation )

    return a dictionnary with same keys than data plus the key 'relU'
        which is a sub-dictionnary containing 2 keys 'Ux' and 'Uy'.
        Ux' and Uy' are 1D array which give the relative velocity of the 
        fluid with respect to the foil center of rotation respectively 
        in the x and y direction.
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation.

    Args:
        data: {}\n
        h: float\n
        r: float\n
        w: float\n

    Returns:
        dataUp : {}\n

    A way you might use me is:\n
        data = add_relU(data,w,h,r)

    """
    if 's' not in data:
            raise ValueError(' data does not have "s" key')

    dataUp = copy.deepcopy(data)
    size=len(dataUp['s'])

    dataUp['Urel']={'Ux':np.zeros(size),'Uy':np.zeros(size),'Umag':np.zeros(size)}

    for i in range(size):
        ( dataUp['Urel']['Ux'][i] , dataUp['Urel']['Uy'][i] , dataUp['Urel']['Umag'][i])=Urel_coord(dataUp['s'][i],h,r,w,Uref,Lambda)

    return (dataUp)



def Urel_coord(s,h,r,w,Uref,Lambda):
    """ read several parameters :
        s : position parameter 0<=s<=1\n
        h : length of vertical translation\n
        r : radius of rotation \n
        w : lenght of horizontal translation\n
        Uref is the velocity impose far upstream of the foil. 
        Uref is the component along x direction (the y-component is nul)

        Lambda is the dimensionless velocity magnitude of the center of rotation 
            Lambda=V/Uref (V being the velocity magnitude of the center of rotation )
    y-axe is the vertical ascending direction
    x-axe is  pointing in downstream direction 
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation. it coordonate are (x,y)=(-r,0).
    Returns the Ux, Uy and Umag.

    Args:
        s : float\n
        h : float\n
        r : float\n
        w : float\n

    Returns:
        (Ux,Uy,Umag) : (float,float,float)\n

    A way you might use me is:\n
        (Ux,Uy,Umag)=Urel_coord(s,h,r,w)

    """
    if s<0 or s>1 :
            raise ValueError(' s parameter is not 0<s<1')


    #distance traveled by the center of rotation in one period :
    L=2*(h+w+pi*r)

    #the velocity magnitude of the center of rotation  :
    V=Lambda*Uref 



    if s<h/L :
        #downward vertical translation part :
        Ux=Uref
        Uy=V

    elif s>=h/L and s<(h+0.5*pi*r)/L :
        #first rotation of pi/2 :
        alpha=(L/r)*s-h/r  #alpha between 0 and pi/2
        Ux=Uref-V*sin(alpha)
        Uy=V*cos(alpha)

    elif s>=(h+0.5*pi*r)/L and s<(h+0.5*pi*r+w)/L :
        #x increasing horizontaltranslation :
        Ux=Uref-V
        Uy=0

    elif s>=(h+0.5*pi*r+w)/L  and s<0.5 :
        #second rotation of pi/2  :
        a=-pi*L/(2*h+2*w+pi*r-L)
        alpha=a*(s-0.5)+0.5*pi #alpha between  0 and pi/2
        Ux=Uref-V*cos(alpha)
        Uy=-V*sin(alpha)

    elif s>=0.5 and s<0.5+h/L :
        #vertical upward translation :
        y=-h+s*L-(h+w+pi*r)
        Ux=Uref
        Uy=-V

    elif s>=0.5+h/L and s<0.5+(h+0.5*pi*r)/L :
        #third rotation of pi/2  :
        alpha=(L/r)*s-h/r-0.5*L/r  #alpha between 0 and pi/2
        Ux=Uref+V*sin(alpha)
        Uy=-V*cos(alpha)

    elif s>=0.5+(h+0.5*pi*r)/L and s<0.5+(h+0.5*pi*r+w)/L  :
        #x decreasing horizontaltranslation :
        Ux=Uref+V
        Uy=0

    elif s>=0.5+(h+0.5*pi*r+w)/L  :
        #s>=1-0.5pi*r/L and s<1
        #fourth rotation of pi/2  :
        alpha=(L/r)*(s-1)+pi/2 #alpha between 0 and pi/2
        Ux=Uref+V*cos(alpha)
        Uy=V*sin(alpha)

    else : 
        raise ValueError('s='+str(s)+' does not belong to any defined trajectory sections')
    Umag=sqrt(Ux**2+Uy**2)

    return (Ux,Uy,Umag)




    
def add_totRot(data,h,r,w,beta0):
    """ read data containing 's','beta' keys and geometric parameters :
        h : length of vertical translation\n
        r : radius of rotation \n
        w : lenght of horizontal translation\n
        beta0 : setting angle on vertical translation parts\n
    return a dictionnary with same keys than data plus the key 'totRot'
        which is a 1D array which gives the z-component of rotation vector 
        between the chord in a position s and the s=0 in degree.
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation.

    Args:
        data: {}\n
        h: float\n
        r: float\n
        w: float\n
        beta0: float\n

    Returns:
        dataUp : {}\n

    A way you might use me is:\n
        data = add_totRot(data,h,r,w,beta0)

    """
    if 's' not in data:
            raise ValueError(' data does not have "s" key')

    dataUp = copy.deepcopy(data)
    size=len(dataUp['s'])

    dataUp['totRot']=np.zeros(size)

    for i in range(size):
        dataUp['totRot'][i]=totRot(dataUp['s'][i],h,r,w,dataUp['beta'][i],beta0)
        

    return (dataUp)



def totRot(s,h,r,w,beta,beta0):
    """ read several parameters :
        s : position parameter 0<=s<=1\n
        h : length of vertical translation\n
        r : radius of rotation \n
        w : lenght of horizontal translation\n
        beta0 : setting angle on vertical translation parts\n
        beta  : setting angle of point with s positon parameters\n
    y-axe is the vertical ascending direction
    x-axe is  pointing in downstream direction 
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation. it coordonate are (x,y)=(-r,0).
    Returns the z-component of rotation vector 
        between the chord in the position s and the s=0 in degree.

    Args:
        s : float\n
        h : float\n
        r : float\n
        w : float\n
        beta0: float\n

    Returns:
        rot : float\n

    A way you might use me is:\n
        rot=totRot(s,h,r,w)

    """
    if s<0 or s>1 :
            raise ValueError(' s parameter is not 0<s<1')


    #distance traveled by the center of rotation in one period :
    L=2*(h+w+pi*r)

    '''rotTang=angle between the trajectory tangent at s 
        and the trajectory tangente at s=0 [degree]    '''

    if s<h/L :
        #downward vertical translation part :
        rotTang=0

    elif s>=h/L and s<(h+0.5*pi*r)/L :
        #first rotation of pi/2 :
        alpha=(L/r)*s-h/r  #alpha between 0 and pi/2
        rotTang=RadToDeg(alpha)

    elif s>=(h+0.5*pi*r)/L and s<(h+0.5*pi*r+w)/L :
        #x increasing horizontaltranslation :
        rotTang=90

    elif s>=(h+0.5*pi*r+w)/L and s<0.5 :
        #second rotation of pi/2  :
        a=-pi*L/(2*h+2*w+pi*r-L)
        alpha=a*(s-0.5)+0.5*pi #alpha between 0 and pi/2
        rotTang=90+RadToDeg(alpha)

    elif s>=0.5 and s<0.5+h/L :
        #vertical upward translation :
        rotTang=180

    elif s>=0.5+h/L and s<0.5+(h+0.5*pi*r)/L :
        #third rotation of pi/2  :
        alpha=(L/r)*s-h/r-0.5*L/r  #alpha between 0 and pi/2
        rotTang=180+RadToDeg(alpha)

    elif s>=0.5+(h+0.5*pi*r)/L and s<0.5+(h+0.5*pi*r+w)/L  :
        #x decreasing horizontaltranslation :
        rotTang=270

    elif s>=0.5+(h+0.5*pi*r+w)/L :
        #s>=1-0.5pi*r/L and s<1
        #fourth rotation of pi/2  :
        alpha=(L/r)*(s-1)+pi/2 #alpha between 0 and pi/2
        rotTang=270+RadToDeg(alpha)

    else : 
        raise ValueError('s='+str(s)+' does not belong to any defined trajectory sections')

    return (beta-beta0+rotTang)



    
def add_totRot_cumul(data):
    """ read data containing 'totRot' key.
    return a dictionnary with same keys than data plus a modify dotRot key.
        which is a 1D array which gives the z-component of rotation vector 
        between the chord in a position s and the s=0 in degree.
        The value of totRot is not any more inferior to 360° but has a continue value.
        During the second blade revolution, 360°<totRot<=720°
    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation.

    Args:
        data: {}\n


    Returns:
        dataUp : {}\n

    A way you might use me is:\n
        data = add_totRot_cumul(data)

    """
    if 'totRot' not in data:
            raise ValueError(' data does not have "totRot" key')

    dataUp = copy.deepcopy(data)
    size=len(dataUp['totRot'])


    for i in range(size):
        dataUp['totRot'][i]=totRot(dataUp['s'][i],h,r,w,dataUp['beta'][i],beta0)
        

    return (dataUp)



def add_zeroRot(data):
    """ read data containing 's' key.
    return a dictionnary with same keys than data plus the key 'totRot'
        which is a 1D array which gives the z-component of rotation vector 
        between the chord in a position s and the s=0 in degree. 
        This angle is fixed, egal to 0.

    The point caracterize by s=0 is taken at the start of the downward 
        vertical translation.

    Args:
        data: {}\n

    Returns:
        dataUp : {}\n

    A way you might use me is:\n
        data = add_zeroRot(data)

    """
    if 's' not in data:
            raise ValueError(' data does not have "s" key')

    dataUp = copy.deepcopy(data)
    size=len(dataUp['s'])

    dataUp['totRot']=np.zeros(size)

    for i in range(size):
        dataUp['totRot'][i]=0

    return (dataUp)




def add_Time(data,T,N,Nstep):
    """ read the dictionnary data and 3 parameters :
        T : period of revolution of one blade (s)\n
        N : number of period to be simulated\n
        Nstep : number of time step per period\n

    return a dictionnary with same keys than data plus the key 'Time'
        which is a 1D array. 

    Args:
        data : dict\n
        T : float\n
        N : int\n
        Nstep : int\n

    Returns:
        data : dict\n

    A way you might use me is:\n
        dataTot=add_Time(data,T,N)

    """
    dataTot = copy.deepcopy(data)
    thetaTot=copy.deepcopy(dataTot['theta'])
    theta0=copy.deepcopy(dataTot['theta'])

    dataTot=concatDict(dataTot,N,popLast=True) 


    Time0=np.linspace(0, T, num=Nstep+1)


    Time=np.linspace(0, T, num=Nstep+1)
    size=len(Time)
    if N>1:
        for i in range(N-1):
            Time=np.concatenate((Time[:-1],Time0))

    Time=make_it_continuous(Time,0,T)
    dataTot['Time']=Time


    dataTot['thetaTot']=make_it_continuous(dataTot['theta'],-pi,pi)

    return (dataTot)


def add_Time_fixedBlade(data,T,N,Nstep):
    """ read the dictionnary data and 2 parameters :
        T : period of revolution of one foil (s)\n
        N : number of period to be simulated\n
        Nstep : number of time step per period\n

    return a dictionnary with same keys than data plus the key 'Time'
        which is a 1D array. 

    Args:
        data : dict\n
        T : float\n
        N : int\n

    Returns:
        data : dict\n

    A way you might use me is:\n
        dataTot=add_Time_fixedBlade(data,T,N)

    """
    dataTot = copy.deepcopy(data)

    dataTot=concatDict(dataTot,N,popLast=True) 

    Time0=np.linspace(0, T, num=Nstep+1)


    Time=np.linspace(0, T, num=Nstep+1)
    size=len(Time)
    if N>1:
        for i in range(N-1):
            Time=np.concatenate((Time[:-1],Time0+(i+1)*T*np.ones(size)))

    dataTot['Time']=Time

    return (dataTot)


def add_multiBlades(data,Nb):
    """ Read data dict and Nb=number of blades. 

    Return a dictionnary with Nb subdict, eatch of those subdict
    contained same keys than data.
    Subdict are named 'blade1', blade2',...
    If Nb=3 :
        blade1 is identical to data. (blade0['s'] start at 0)
        blade2['s'] start at 1/3.
        blade3['s'] start at 2/3.
         
    
    Args:
        data : dict\n
        Nb : int\n

    Returns:
        dataBld : dict\n

    A way you might use me is:\n
        dataBld=add_multiBlades(data,Nb)

    """
    dataBld={}

    for i in range(Nb):
        datatmp = copy.deepcopy(data)
        dataBld['blade'+str(i+1)]=shift_data(datatmp,'s',i/Nb)
                            

    return (dataBld)


def init_multiBlades_data(s,Nb):
    """ Read data s=1D array and Nb=number of blades. 

    Return a dictionnary with Nb subdict named 'blade1', blade2',...
    Itch subdict contains key 's', a 1D array define as follow :
    If Nb=3 :
        blade1['s'] is identical to s (s[0]=s[-1]=0)
        blade2['s'] start at 1/3.
        blade3['s'] start at 2/3.
         
    
    Args:
        s : array [float]\n
        Nb : int\n

    Returns:
        data : dict\n

    A way you might use me is:\n
        data=init_multiBlades_data(s,Nb)

    """
    data={}

    for i in range(Nb):
        s_tmp = copy.deepcopy(s)
        key='blade'+str(i+1)
        data[key]={}
        data[key]['s']=shift_array(s,i/Nb)
                            

    return (data)



def write_motionData(dataBld,path_output):
    """ dataBld is a dictionary of subdictionaries.
        Eatch of it sub dict has following keys :
            ['Time']
            ['CofR']['x']
            ['CofR']['y']
            ['thetaTot']
            they are 1D array of them size.
        It will write one motion file for every subdict in dataBld.
        
            
    Args:
        dataBld: dict{array}\n
        output_file: str\n

    A way you might use me is:\n
        write_motionData(dataBld,'./constant')
    """

    for sub in dataBld :

        if 'Time' not in dataBld[sub]:
            raise ValueError('Time not in dataBld['+sub+'] keys')
        if 'CofR' not in dataBld[sub]:
            raise ValueError('CofR not in dataBld['+sub+'] keys')
        if 'thetaTot' not in dataBld[sub]:
            raise ValueError('thetaTot not in dataBld['+sub+'] keys')


        size=len(dataBld[sub]['Time'])

        x0=dataBld[sub]['CofR']['x'][0]
        y0=dataBld[sub]['CofR']['y'][0]
        theta0=dataBld[sub]['thetaTot'][0]

        t=dataBld[sub]['Time']
        x=dataBld[sub]['CofR']['x']-np.ones(size)*x0
        y=dataBld[sub]['CofR']['y']-np.ones(size)*y0

        theta=dataBld[sub]['thetaTot']-np.ones(size)*theta0
        theta=[RadToDeg(a) for a in theta]

        txt = '''{size}
        (
        '''.format(size=size)

        for i in range(size):
            txt+='''({t}        ( ({x} {y} 0.0)     (0 0 {theta} ) ) )\n'''.format(t=t[i],x=x[i],y=y[i],theta=theta[i])

        txt+='''
    )'''

    
        path = os.path.join(path_output,'motion_'+sub+'.dat')
        print('Save motion_'+sub+'.dat in path ' + path)

        with open(path, 'w') as f:
            f.write(txt)



def write_dynamicMeshDict(dataBld,path_output,path_wd):
    """ this function write the constant/dynamicMeshDict file
        for a multi blade simulation.
        dataBld is dictionary of subdict. Eatch subdic contain data 
        of 1 blade motion.
        path_output is the relative path to constant directory
        path_wd is the absolute path of the working directory
        
            
    Args:
        dataBld : {}\n
        path_output: str\n
        path_wd: str\n

    A way you might use me is:\n
        write_dynamicMeshDict(dataBld,'./constant',path_wd)
    """

    txt = '''/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2012                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    class       dictionary;
    object      dynamicMeshDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

dynamicFvMesh       dynamicOversetFvMesh;

solver              multiSolidBodyMotionSolver;

multiSolidBodyMotionSolverCoeffs
{'''
    for sub in dataBld :
        
        x0 = dataBld[sub]['CofR']['x'][0]
        y0 = dataBld[sub]['CofR']['y'][0]
        absPathMotionData = os.path.join(path_wd,path_output, f'motion_{sub}.dat')

        txt += f'''\n
    set_{sub}
    {{
       solidBodyMotionFunction tabulated6DoFMotion;
       tabulated6DoFMotionCoeffs
       {{
           CofG            ( {x0} {y0} 0.0 );
           timeDataFileName    "{absPathMotionData}";

        }}
    }}'''

    txt += '''

}

// ************************************************************************* //'''
    path = os.path.join(path_output,'dynamicMeshDict')
    print('Save dynamicMeshDict in path ' + path)

    with open(path, 'w') as f:
        f.write(txt)
    


def pop_UselessTime(data,r,h):
    """ data is a dictionary. It has following keys :
            ['CofR']['x']
            ['CofR']['y']
            they are 1D array of them size.
        solidBodyMotionFunction tabulated6DoFMotion use 
        interpolated values of motion.dat file. 
        Vertical translation parts with constant setting angle        
        therefor do not recuire more than 2 points to describe them 
        (starting and ending points).
        Here we delete this useless information
        
                        
            
    Args:
        data: dict{array}\n
        r: float\n
        h: float\n

    A way you might use me is:\n
        data=pop_UselessTime(data,r,h)
    """

    dataMin = copy.deepcopy(data)
    size=len(data['CofR']['x'])
    delList=[]

    for i in range(size):
        delete=False
        if dataMin['CofR']['x'][i]==-r or dataMin['CofR']['x'][i]==r: delete=True
        if dataMin['CofR']['y'][i]<-0.95*h or dataMin['CofR']['y'][i]>-0.05*h: delete=False
        if delete :
            delList.append(i) 

    for key in dataMin:
        if type(dataMin[key])==dict:
            for subkey in dataMin[key]:
                dataMin[key][subkey]=np.delete(dataMin[key][subkey],delList)
        else:
            dataMin[key]=np.delete(dataMin[key],delList)

    return (dataMin)



def write_init_pos(data,dir_path,file_name):
    """ this function write position (x,y,theta) of every blade
        at t=0.
        The text file generated will be use by template_Allrun_airfoilMesh.py
        
            
    Args:
        data : {}\n
        dir_path: str\n
        file_name: str\n

    A way you might use me is:\n
        write_init_pos(data,'../constant','blades_init_pos.dat')
    """
    
    data0={}
    data0['x']=[]
    data0['y']=[]
    data0['theta']=[]

    for key in data :
        data0['x'].append(data[key]['CofR']['x'][0])
        data0['y'].append(data[key]['CofR']['y'][0])
        data0['theta'].append(data[key]['theta'][0])

    write_data(data0,dir_path,file_name,decimal=10)

    return None






























