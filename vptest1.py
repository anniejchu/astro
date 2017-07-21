'''
from visual import *
import numpy as np 
from astropy import constants as const 
from astropy import units
import pdb
from astropy.time import Time 
import time

#setup
scene.width=800
scene.height=600

#variables
mL = 5 #mass of lens in solar masses
xS0 = vector(-5,0,0) #initial location of the star
dL = 4000 #distance to the lens in pc
dS = 8000 #distance to the star in pc
     #proper motion: maybe input the the radial velocity and tangetial velocity
#muS = pass #diff_angle(v1,v2) to find the angle between 2 vectors
#muL = pass 

#conversion factors and constants
kappa = 4.0 * const/G * units.rad / const.c**2*units.au
kappa = kappa.to(units.mas / units.solMass)

#calculating the einstein ring
inv_dist_diff = (1.0 / (dL * units.pc)) - (1.0 / (dS * units.pc))
thetaE = units.rad * np.sqrt((4.0 * const.G * mL * units.M_sun / const.c**2) * inv_dist_diff)


#source positioning 
sPos = vector(-1000.0, 200.0, -dS)
STAR = sphere(pos=sPos, radius = 100, color = color.yellow)
STAR.velocity = vector(250, 0, 0)
deltaT = 0.005
t = -10

#lens positioning
lPos = vector(0, 0, -dL)
LENS = sphere(pos = lPos, radius = 100, color=color.blue)
#ER = ring(pos = lPos, radius = -- ) #BRING IN EQUATION WITH DR. LU'S CODE


def movingsource():
	while t < -t: 
		STAR.pos.x = STAR.pos.x + STAR.velocity*deltaT
		t = t+deltaT
		rate(100)

print(thetaE)

'''

from visual import *
import numpy as np

#setup
scene.width = 800
scene.heigh = 600

#variables 
imL = 5 #solar mass MUST CONVERT
idL = 4000 #pc
idS = 8000 #pc

# conversion
mL = imL * (1.99 * (10 ** 30))
dL = idL * (30856775714409184)
dS = idS * (30856775714409184)
G = 6.67 * 10**-11
c = 3.0 * 10**8
radtomas = 206264806.247

#calculating einstain radius in radians
inv_dist_diff = (1.0/dL)-(1.0/dS)
thetaE = (np.sqrt((4.0*G*mL/c**2) * inv_dist_diff))*radtomas


print(thetaE)