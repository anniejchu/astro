'''
from visual import *
r = vector(-3,4,0)

circle1 = sphere(pos = vector(1,1,1), radius = 0.25, color = color.cyan)
circle2 = sphere(pos = r, radius = 0.15, color = color.red)
pointer = arrow(pos = circle1.pos, axis=circle2.pos-circle1.pos, color = color.blue)


while circle2.x < 10:
	#sphere(pos=r, radius=0.5, color = color.cyan)
	rate(10)
	#circle2.pos = r
	circle2.x = circle2.x +1

'''

from visual import *
import numpy as np 

#setup
scene.width=800
scene.height=600

#constants
c = 3.00*(10**8) #m/s
G = 6.67408*(10**-11) #m^3/kgs^2

#variables
mL = 5 #mass of lens in solar masses
xS0 = vector(-5,0,0) #initial location of the star
dL = 4000 #distance to the lens in pc
dS = 8000 #distance to the star in pc
     #proper motion: maybe input the the radial velocity and tangetial velocity
#muS = pass #diff_angle(v1,v2) to find the angle between 2 vectors
#muL = pass 

'''

#calculations for angular einsten ring
mL = mL*(7.3477*10**22) #solar mass to kg
dL = dL*(3.0857*10**16)
dS = dL*(3.0857*10**16)

DIdist = (dL**-1)-(dS**-1)
pt1 = (4*G*mL)/(c**2)
thetaE = np.sqrt(DIdist*pt1)
'''

#-----------------------------------------
LENS = sphere(pos = vector(0,0, dL), radius = 30, color=color.green)
