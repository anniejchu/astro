
from visual import *
import numpy as np

#random variables
origin = vector(0,0,0)
#setup
scene.width = 200
scene.height = 100

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
mtopc = 30856776000000000


#calculating einstein radius 
inv_dist_diff = (1.0/dL)-(1.0/dS)
thetaE = (np.sqrt((4.0*G*mL/c**2) * inv_dist_diff)) # radians
thetaE1 = thetaE * radtomas #in mas

#finding einstein radius distance from center
degreeE = (thetaE * 180)/np.pi # radians to degrees
eradiusm = np.tan(degreeE)*dL #meters
eradiuspc = eradiusm/mtopc
eradiuspc_adjusted = eradiuspc*500000

#lens positioning
lPos = vector(0, 0, -idL)
LENS = sphere(pos = lPos, radius = 100, color=color.blue)
ER = ring(pos = lPos, radius = eradiuspc_adjusted, axis = (0,0,1), thickness=50, color= color.white)

#source positioning 
sPos = vector(-3000.0, 200.0, -idS)
STAR = sphere(pos=sPos, radius = 100, color = color.yellow)
STAR.velocity = vector(250, 0, 0)

#deltaT = speed , xlimit = 4000 
'''
def movingsource(deltaT, xlimit):
	while STAR.pos.x < xlimit: 
		STAR.pos = STAR.pos + STAR.velocity*deltaT
		rate(100)
'''
def movingsource(t, deltaT, timelimit):
	#FIND EINSTEIN CROSSING TIME, MAKE IT RELEVANT TO THIS FUNCTION
	#REMEMBER THAT THE EINSTEIN CROSSING TIME IS JUST HALF OF THE GRAPH
	while t < timelimit:
		STAR.pos = STAR.pos + STAR.velocity*deltaT
		t = t+deltaT
		rate(100)

movingsource(-20, 0.05, 50)

#light curve info ASK NIJAID WHAT THIS ANGLE IS
	#radians
	#BE ABLE TO USE THIS ANGLE TO CREATE A VECTOR USE TRIG????
thetaS = diff_angle(STAR.pos, origin) 
thetaPOS = (thetaS + np.sqrt((thetaS**2)+(4*thetaE)))/2
thetaNEG = (thetaS - np.sqrt((thetaS**2)+(4*thetaE)))/2

	#mas
thetaS1 = thetaS * radtomas
thetaPOS1 = (thetaS1 + np.sqrt((thetaS1**2)+(4*thetaE1)))/2
thetaNEG1= (thetaS1 - np.sqrt((thetaS1**2)+(4*thetaE1)))/2


#light magnification 
#USE THIS TO EITHER CHANGE THE COLOR OR INTENSITY(OPAQUENESS)
	#radians
u = diff_angle(STAR.pos, LENS.pos)
amp = ((u**2)+2)/(u*(np.sqrt((u**2)+4)))
	#mas
u1 = u *radtomas
amp1 = ((u1**2)+2)/(u1*(np.sqrt((u1**2)+4)))


